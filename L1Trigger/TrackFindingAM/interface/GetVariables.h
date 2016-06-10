#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <math.h>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include "../interface/L1TrackTriggerTree.h"
#include "../interface/GetTrackParameters.h"
#include "../interface/StubsCombination.h"


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta
template <class T>
double extrapolateR(const double & R, const double & z, const int layer, const double & tgTheta,
                    const std::vector<int> & uniqueLayers, const std::vector<double> & originalR, const T & originalZ)
{
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation.
template <class T>
double extrapolateRSecondOrder(const double & R, const double & z, const int layer, const double & tgTheta,
                               const double & chargeOverTwoRho,
                               const std::vector<int> & uniqueLayers,
                               const std::vector<double> & originalR, const T & originalZ)
{
  return R;
}


template <class T>
double extrapolateRExact(const double & R, const double & z, const int layer, const double & tgTheta,
                         const double & chargeOverTwoRho,
                         const std::vector<int> & uniqueLayers,
                         const std::vector<double> & originalR, const T & originalZ)
{
  return R;
}


double correctPhiForNonRadialStrips(const double & phi, const double & stripPitch, const float & stripIndex,
                                    const double & extrapolatedR, const double & R, const int layer);


class CorrectPhiForNonRadialStripsLookup {
 public:
  CorrectPhiForNonRadialStripsLookup() {}
  double correctPhiForNonRadialStrips(const double &phi, const double &stripPitch,
                                      const double &extrapolatedR, const double &R,
                                      const double & z,
                                      const int layer) const
  {
    return phi;
  }
};


// Abstract base class
class GetTreeVariable
{
public:
  GetTreeVariable(const std::string & name, const std::set<int> & layers) : name_(name), layers_(layers) {}
  virtual ~GetTreeVariable() {}
  virtual double at(const int k) = 0;
  // virtual double at(const int k) = 0;
  virtual bool layer(const int layer) { return (layers_.count(layer) != 0); }
  unsigned int layersNum() { return layers_.size(); }
  const std::set<int> * const layers() { return &layers_; }
  std::string name() const { return name_; }
  // void resetSeed() { generator_.seed(0); }

protected:
  std::string name_;
  // std::unordered_set<int> layers_;
  std::set<int> layers_;
  std::default_random_engine generator_;
};


class ParametrizedMagneticField
{
 public:
  ParametrizedMagneticField():
      c1(3.8114),
      b0(-3.94991e-06),
      b1(7.53701e-06)
      // , a (2.43878e-11)
  {}
  inline double B0Z(const double z) const {
    return b0*z*z + b1*z + c1;
  }

 private:
  double c1;
  double b0;
  double b1;
  // double a;
};


// Phi variable of the stubs
class GetVarPhi : public GetTreeVariable
{
public:
  GetVarPhi(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarPhi() {}
  virtual double at(const int k) {return std::atan2(var_y->at(k), var_x->at(k));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


// Z variable of the stubs
class GetVarZ : public GetTreeVariable
{
public:
  GetVarZ(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_z(tree->m_stub_z) {}
  virtual ~GetVarZ() {}
  virtual double at(const int k) {return var_z->at(k);}
private:
  std::vector<float> * var_z;
};


// R variable of the stubs
class GetVarR : public GetTreeVariable
{
public:
  GetVarR(const std::string & name, std::shared_ptr<L1TrackTriggerTree> tree, const std::set<int> & layers) :
      GetTreeVariable(name, layers), var_x(tree->m_stub_x), var_y(tree->m_stub_y) {}
  virtual ~GetVarR() {}
  virtual double at(const int k) {return std::sqrt(std::pow(var_x->at(k), 2) + std::pow(var_y->at(k), 2));}
private:
  std::vector<float> * var_x;
  std::vector<float> * var_y;
};


class EstimatorSimple
{
 public:
  EstimatorSimple(const TString & inputFileName) {
    // open matrix file and read V and D arrays
    // std::cout << "opening " + inputFileName + " for reading" << std::endl;

    std::ifstream inputFile;
    inputFile.open(inputFileName);
    if (!inputFile) {
      // std::cout << "EstimatorSimple: Error opening " + inputFileName << std::endl;
      return;
      // throw;
    }

    // Read number of variables and number of track parameters
    int nVars = 0;
    inputFile >> nVars;

    // Skip required layers
    int l;
    for (int v = 0; v < nVars; ++v) {
      inputFile >> l;
    }

    // Read mean values
    double x;
    for (int i = 0; i < nVars; ++i) {
      inputFile >> x;
      means_.push_back(x);
    }
    // Read parameter mean value
    inputFile >> parameterMean_;
    // std::cout << "parameterMean_ = " << parameterMean_ << std::endl;

    // Read coefficients
    for (int i = 0; i < nVars; ++i) {
      inputFile >> x;
      coeff_.push_back(x);
    }
  }

  template <class T>
  double estimate(const T & var)
  {
    double estimatedParameter = 0.;
    for (unsigned int i=0; i<var.size(); ++i) {
      estimatedParameter += (var[i]-means_[i])*coeff_[i];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }

  template <class T, class U>
  double estimate(const T & var1, const U & var2)
  {
    double estimatedParameter = 0.;
    for (unsigned int i=0; i<var1.size(); ++i) {
      estimatedParameter += (var1[i]-means_[i*2])*coeff_[i*2];
      estimatedParameter += (var2[i]-means_[i*2+1])*coeff_[i*2+1];
    }
    // When it is estimated the mean value is subtracted. We add it back.
    return (estimatedParameter + parameterMean_);
  }
 private:
  std::vector<double> means_;
  std::vector<double> coeff_;
  double parameterMean_;
};


class TransformBase
{
 public:
  TransformBase(const std::string & name) :
      name_(name)
  {}
  TransformBase(const std::string & name, const std::string & preEstimateFileName,
                const std::vector<double> & meanRadius) :
      name_(name), estimator_(std::make_shared<EstimatorSimple>(preEstimateFileName)), meanRadius_(meanRadius)
  {}
  TransformBase(const std::string & name, const std::vector<double> & meanRadius) :
      name_(name), meanRadius_(meanRadius)
  {}
  TransformBase(const std::string & name,
                const std::vector<double> & meanRadius, const std::vector<double> & meanZ) :
      name_(name), meanRadius_(meanRadius), meanZ_(meanZ)
  {}
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const = 0;
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const = 0;
  std::string getName() { return name_; }
 protected:
  std::string name_;
  std::shared_ptr<EstimatorSimple> estimator_;
  std::vector<double> meanRadius_;
  std::vector<double> meanZ_;
};


class TransformPropagatePhi : public TransformBase
{
 public:
  TransformPropagatePhi(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagatePhi() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    return stubsCombination.phi(index);
  }
};


class TransformPropagateR : public TransformBase
{
 public:
  TransformPropagateR(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagateR() {}
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    return stubsCombination.R(index);
//    return vars.at(index*3+1);
  }
};


class TransformPropagateZ : public TransformBase
{
 public:
  TransformPropagateZ(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagateZ() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    return stubsCombination.z(index);
//    return vars.at(index*3+2);
  }
};


class TransformCorrectedPhiFirstOrder : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrder(const std::string & name, const std::string & preEstimateChargeOverPtFileName,
                                  const std::vector<double> & meanRadius) :
      TransformBase(name, preEstimateChargeOverPtFileName, meanRadius)
  {}
  virtual ~TransformCorrectedPhiFirstOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    std::vector<double> originalPhi;
//    for (const std::vector<Stub>::const_iterator it = stubsCombination.begin(); it != stubsCombination.end(); ++it) {
//      originalPhi.push_back(it->phi());
//    }
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double DeltaR = R - meanRadius_[index];
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2.);
  }
};


class TransformCorrectedPhiSecondOrder : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrder(const std::string & name, const std::string & preEstimateChargeOverPtFileName,
                                   const std::vector<double> & meanRadius) :
  TransformBase(name, preEstimateChargeOverPtFileName, meanRadius)
  {}
  virtual ~TransformCorrectedPhiSecondOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    std::vector<double> originalPhi;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    // double deltaPhi = estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.;
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderExtrapolatedR : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedR(const std::string & name,
                                                const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                const std::string & firstOrderTgThetaCoefficientsFileName,
                                                const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedR() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index*3+1);
    double z = stubsCombination.z(index*3+2);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size()/3; ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
//    double cotTheta = estimatorTgTheta_->estimate(originalR, originalZ);
//    double genChargeOverTwoRho = genChargeOverPt*3.8114*0.003/2.;

    // If this is a 2S module in the disks
    R = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
//    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
//    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2.);
//    return (phi + genChargeOverPt*DeltaR*3.8114*0.003/2.);
//    return (phi + genChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(genChargeOverPt*3.8114*0.003/2., 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder(const std::string & name,
                                                           const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                           const std::string & firstOrderTgThetaCoefficientsFileName,
                                                           const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
//    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;

    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);

    // If this is a 2S module in the disks
    R = extrapolateRSecondOrder(R, z, stubsCombination.layer(index), tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection(const std::string & name,
                                                           const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                           const std::string & firstOrderTgThetaCoefficientsFileName,
                                                           const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrection() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
//    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);

    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStrips(phi, 0.009, stubsCombination.stub(index).strip(), extrapolatedR, R, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup(const std::string & name,
                                                                                         const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                         const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                         const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
//    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;

    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);

    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    // phi = correctPhiForNonRadialStrips(phi, 0.009, stubsCombination.stub(index).strip(), extrapolatedR, R, layer);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_genTheta : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_genTheta(const std::string & name,
                                                                                                  const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                                  const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                                  const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_genTheta() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
//    double estimatedChargeOverPt = estimator_->estimate(originalPhi);
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;

//    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    double tgTheta = 1./stubsCombination.genCotTheta();

    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    // phi = correctPhiForNonRadialStrips(phi, 0.009, stubsCombination.stub(index).strip(), extrapolatedR, R, layer);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN(const std::string & name,
                                                                                             const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                             const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                             const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, 1./stubsCombination.genCotTheta(), chargeOverTwoRho,
                                                   stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
//    if (stubsCombination.layer(index) > 10 && R > 61.) {
//      double exactExtrapolatedR = sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
//      std::cout << "extrapolatedR = " << extrapolatedR << std::endl;
//      std::cout << "exactExtrapolatedR = " << exactExtrapolatedR << std::endl;
//      std::cout << std::endl;
//    }
    return phi;
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN_ExactExtrapolation : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN_ExactExtrapolation(const std::string & name,
                                                                                                                const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                                                const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                                                const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup_GEN_ExactExtrapolation() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = R;
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      extrapolatedR = sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
    }
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedZFirstOrder : public TransformBase
{
 public:
  TransformCorrectedZFirstOrder(const std::string & name,
                                const std::string & firstOrderCotThetaCoefficientsFileName,
                                const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius)
  {}
  virtual ~TransformCorrectedZFirstOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    return (z - DeltaR*cotTheta);
  }
 private:
};


class TransformCorrectedZSecondOrder : public TransformBase
{
 public:
  TransformCorrectedZSecondOrder(const std::string & name,
                                 const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                 const std::string & firstOrderCotThetaCoefficientsFileName,
                                 const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPt_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPtCoefficientsFileName))
  {}
  virtual ~TransformCorrectedZSecondOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double oneOverRho = (3.8114*0.003)*estimatorChargeOverPt_->estimate(originalPhi);
    return (z - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*cotTheta);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


class TransformCorrectedPhiFirstOrderPz : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrderPz(const std::string & name,
                                    const std::string & firstOrderChargeOverPzCoefficientsFileName,
                                    const std::string & firstOrderCotThetaCoefficientsFileName,
                                    const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPz_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPzCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiFirstOrderPz() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double oneOverRhoz = (3.8114*0.003)*estimatorChargeOverPz_->estimate(originalPhi);
    return (phi + oneOverRhoz*cotTheta*DeltaR*3.8114*0.003/2.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPz_;
};


// Using generator-level c/pT
class TransformCorrectedPhiSecondOrderGen : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGen() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + stubsCombination.genChargeOverPt()*DeltaR*3.8114*0.003/2. + RCube*std::pow(stubsCombination.genChargeOverPt()*3.8114*0.003/2., 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderGenDeltaZ : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGenDeltaZ(const std::string & name,
                                            const std::vector<double> & meanRadius,
                                            const std::vector<double> & meanZ) :
      TransformBase(name, meanRadius, meanZ)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGenDeltaZ() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    // Correct the R for the z variation. This is needed only in the 2S modules of the disks
    // where the R resolution is too low to see the DeltaZ variation within the disks.
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      // std::cout << "R = " << R << std::endl;
      // R must increase when moving from smaller z to bigger z. If z < meanZ the sign of the added them must be positive.
      R += (meanZ_[index] - z)/stubsCombination.genCotTheta();
      // std::cout << "adjusted R = " << R << std::endl;
      // std::cout << "gen R = " << sin((z - genZ0)/genCotTheta*chargeOverTwoRho)/chargeOverTwoRho << std::endl;
    }
    double DeltaR = R - meanRadius_[index];
    return (phi + DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderGenExactR : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGenExactR(const std::string & name,
                                            const std::vector<double> & meanRadius,
                                            const std::vector<double> & meanZ) :
      TransformBase(name, meanRadius, meanZ)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGenExactR() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      R = sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
    }
    // double meanR = meanRadius_[index];
    double DeltaR = R - meanRadius_[index];
    // double deltaPhi = DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.;
    return (phi + DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.);
  }
};


class TransformCorrectedPhiSecondOrderGenExactRNonRadialStripCorrection : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderGenExactRNonRadialStripCorrection(const std::string & name,
                                                                    const std::vector<double> & meanRadius,
                                                                    const std::vector<double> & meanZ) :
      TransformBase(name, meanRadius, meanZ)
  {}
  virtual ~TransformCorrectedPhiSecondOrderGenExactRNonRadialStripCorrection() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      double extrapolatedR =
          sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
      // Correct the phi coordinate using the extrapolated R at the given z
      // strip pitch = 0.009 cm
      // phi = phi - 0.009*(stubsCombination.stub(index).strip() - 512)*(extrapolatedR - R)/(R*extrapolatedR);
      phi = correctPhiForNonRadialStrips(phi, 0.009, stubsCombination.stub(index).strip(), extrapolatedR, R,
                                         stubsCombination.layer(index));
      // phi = phi - 0.009*(stubsCombination.stub(index).strip() - 512)*(extrapolatedR - R)/(R*R);
      // Replace the R for later use
      R = extrapolatedR;
    }
    // double meanR = meanRadius_[index];
    double DeltaR = R - meanRadius_[index];
    // double deltaPhi = DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.;
    return (phi + DeltaR*chargeOverTwoRho + std::pow(R*chargeOverTwoRho, 3)/6.);
  }
};


class TransformCorrectedPhiExactGen : public TransformBase
{
 public:
  TransformCorrectedPhiExactGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiExactGen() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
  }
};


class TransformCorrectedPhiExactGenExactR : public TransformBase
{
 public:
  TransformCorrectedPhiExactGenExactR(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedPhiExactGenExactR() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    // Compute R using z
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      R = sin((z - stubsCombination.genZ0()) / stubsCombination.genCotTheta() * chargeOverTwoRho) / chargeOverTwoRho;
    }
    return (phi + asin(R*chargeOverTwoRho) - asin(meanRadius_[index]*chargeOverTwoRho));
  }
};


// Using generator-level c/pT and cot(theta)
class TransformCorrectedZSecondOrderGen : public TransformBase
{
 public:
  TransformCorrectedZSecondOrderGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedZSecondOrderGen() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double DeltaR = R - meanRadius_[index];
    double oneOverRho = (3.8114*0.003)*stubsCombination.genChargeOverPt();
    return (z - (DeltaR + 1/24.*std::pow(R, 3)*(oneOverRho*oneOverRho))*stubsCombination.genCotTheta());
  }
};


class TransformCorrectedZExactGen : public TransformBase
{
 public:
  TransformCorrectedZExactGen(const std::string & name, const std::vector<double> & meanRadius) :
      TransformBase(name, meanRadius)
  {}
  virtual ~TransformCorrectedZExactGen() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double rho = 1./((3.8114*0.003)*stubsCombination.genChargeOverPt());
    return (z - 2*rho*stubsCombination.genCotTheta()*asin(R/(2*rho)) + 2*rho*stubsCombination.genCotTheta()*asin(meanRadius_[index]/(2*rho)));
  }
};


class TransformExtrapolatedR : public TransformBase
{
 public:
  TransformExtrapolatedR(const std::string & name,
                         const std::string & firstOrderTgThetaCoefficientsFileName,
                         const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderTgThetaCoefficientsFileName, meanRadius)
  {}
  virtual ~TransformExtrapolatedR() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double tgTheta = estimator_->estimate(originalR, originalZ);
//    double tgTheta = 1./genCotTheta;
//    double extrapolatedR = extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
//    return extrapolatedR;
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
    return extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ, genZ0, genChargeOverPt);

//    double chargeOverTwoRho = (3.8114*0.003)*genChargeOverPt/2.;
//    return sin((z - genZ0) * chargeOverTwoRho / genCotTheta) / chargeOverTwoRho;

  }
 private:
};


class TransformExtrapolatedRSecondOrder : public TransformBase
{
 public:
  TransformExtrapolatedRSecondOrder(const std::string & name,
                                    const std::string & firstOrderTgThetaCoefficientsFileName,
                                    const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                    const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderTgThetaCoefficientsFileName, meanRadius),
      estimatorChargeOverPt_(std::make_shared<EstimatorSimple>(firstOrderChargeOverPtCoefficientsFileName))
  {}
  virtual ~TransformExtrapolatedRSecondOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double tgTheta = estimator_->estimate(originalR, originalZ);
//    double tgTheta = 1./estimator_->estimate(originalR, originalZ);
//    double chargeOverPt = estimatorChargeOverPt_->estimate(originalPhi);
    double chargeOverTwoRho = estimatorChargeOverPt_->estimate(originalPhi)*3.8114*0.003/2.;

    //    double tgTheta = 1./genCotTheta;
//    double extrapolatedR = extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
//    return extrapolatedR;
//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ);
    double extrapolatedR = extrapolateRSecondOrder(R, z, stubsCombination.layer(index), tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
//    double genDistanceR = stubsCombination.genTrackDistanceLongitudinalR(index);
//    StubsCombination transformedStubsCombination(stubsCombination);
//    transformedStubsCombination.setR(index, extrapolatedR);
//    double genDistanceExtrapolatedR = transformedStubsCombination.genTrackDistanceLongitudinalR(index);
//    std::cout << "R = " << R << "genTrackDistanceR = " << genDistanceR << std::endl;
//    std::cout << "extrapolatedR = " << extrapolatedR << ", genTrackDistanceExtrapolatedR = " << genDistanceExtrapolatedR << std::endl;
    return extrapolatedR;
//    return extrapolateRSecondOrder(R, z, stubsCombination.layer(index), tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);

//    return extrapolateR(R, z, uniqueLayers[index], tgTheta, uniqueLayers, originalR, originalZ, genZ0, genChargeOverPt);

//    double chargeOverTwoRho = (3.8114*0.003)*genChargeOverPt/2.;
//    return sin((z - genZ0) * chargeOverTwoRho / genCotTheta) / chargeOverTwoRho;

  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


class TransformExtrapolatedRExact : public TransformBase
{
 public:
  TransformExtrapolatedRExact(const std::string & name) :
      TransformBase(name)
  {}
  virtual ~TransformExtrapolatedRExact() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
//  virtual double operator()(const int index, const std::vector<double> & vars, const std::vector<int> & uniqueLayers,
//                            const double & genChargeOverPt, const double & genCotTheta, const double & genZ0) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i = 0; i < stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double chargeOverTwoRho = (3.8114 * 0.003) * stubsCombination.genChargeOverPt() / 2.;
//    double exactExtrapolatedR = extrapolateRExact(R, z, stubsCombination.layer(index), 1./stubsCombination.genCotTheta(), chargeOverTwoRho,
//                             stubsCombination.layers(), originalR, originalZ);
//    double extrapolatedR = extrapolateRSecondOrder(R, z, stubsCombination.layer(index), 1./stubsCombination.genCotTheta(), chargeOverTwoRho,
//                                                   stubsCombination.layers(), originalR, originalZ);
//    return exactExtrapolatedR;
    return extrapolateRExact(R, z, stubsCombination.layer(index), 1./stubsCombination.genCotTheta(), chargeOverTwoRho,
                             stubsCombination.layers(), originalR, originalZ);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


#endif // GETVARIABLES_H
