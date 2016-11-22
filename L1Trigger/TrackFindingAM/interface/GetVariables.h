#ifndef GETVARIABLES_H
#define GETVARIABLES_H

#include <memory>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <stdexcept>
#include "L1Trigger/TrackFindingAM/interface/L1TrackTriggerTree.h"
#include "L1Trigger/TrackFindingAM/interface/GetTrackParameters.h"
#include "L1Trigger/TrackFindingAM/interface/StubsCombination.h"


double computeRotationFactor(const std::vector<double> & vars);


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta
template <class T>
double extrapolateR(const double & R, const double & z, const int layer, const double & tgTheta,
                    const std::vector<int> & uniqueLayers, const std::vector<double> & originalR, const T & originalZ)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        return (originalR[i] + (z - originalZ[i]) * tgTheta);
      }
    }
  }
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation with only two terms.
template <class T>
double extrapolateRSecondOrderFirstTwoTermsOnly(const double & R, const double & z, const int layer, const double & tgTheta,
                               const double & chargeOverTwoRho,
                               const std::vector<int> & uniqueLayers,
                               const std::vector<double> & originalR, const T & originalZ,
                               std::vector<double> & firstOrderTerm,
                               std::vector<double> & secondOrderTerm1,
                               std::vector<double> & secondOrderTerm2)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i]) * tgTheta;
        firstOrderTerm.push_back(deltaZTgTheta);
        double term1 = -deltaZTgTheta * std::pow(originalR[i], 2);
        secondOrderTerm1.push_back(term1 * std::pow(chargeOverTwoRho, 2) / 2.);
        double term2 = -originalR[i] * std::pow(deltaZTgTheta, 2);
        secondOrderTerm2.push_back(term2 * std::pow(chargeOverTwoRho, 2) / 2.);
        double secondOrderTerm = (term1 + term2) * std::pow(chargeOverTwoRho, 2) / 2.;
        return (originalR[i] + deltaZTgTheta + secondOrderTerm);
      }
    }
  }
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation with only one term.
template <class T>
double extrapolateRSecondOrderFirstTermOnly(const double & R, const double & z, const int layer, const double & tgTheta,
                                            const double & chargeOverTwoRho,
                                            const std::vector<int> & uniqueLayers,
                                            const std::vector<double> & originalR, const T & originalZ,
                                            std::vector<double> & firstOrderTerm,
                                            std::vector<double> & secondOrderTerm1)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i]) * tgTheta;
        firstOrderTerm.push_back(deltaZTgTheta);
        double term1 = -deltaZTgTheta * std::pow(originalR[i], 2);
        secondOrderTerm1.push_back(term1 * std::pow(chargeOverTwoRho, 2) / 2.);
        double secondOrderTerm = term1 * std::pow(chargeOverTwoRho, 2) / 2.;
        return (originalR[i] + deltaZTgTheta + secondOrderTerm);
      }
    }
  }
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a first order approximation.
template <class T>
double extrapolateRFirstOrder(const double & R, const double & z, const int layer, const double & tgTheta,
                              const double & chargeOverTwoRho,
                              const std::vector<int> & uniqueLayers,
                              const std::vector<double> & originalR, const T & originalZ,
                              std::vector<double> & firstOrderTerm)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i]) * tgTheta;
        firstOrderTerm.push_back(deltaZTgTheta);
        return (originalR[i] + deltaZTgTheta);
      }
    }
  }
  return R;
}


// Extrapolate R to from the outermost PS module to the given z position using the given tgTheta.
// Use a second order approximation.
template <class T>
double extrapolateRSecondOrder(const double & R, const double & z, const int layer, const double & tgTheta,
                               const double & chargeOverTwoRho,
                               const std::vector<int> & uniqueLayers,
                               const std::vector<double> & originalR, const T & originalZ,
                               std::vector<double> & firstOrderTerm,
                               std::vector<double> & secondOrderTerm1,
                               std::vector<double> & secondOrderTerm2,
                               std::vector<double> & secondOrderTerm3)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i])*tgTheta;
        firstOrderTerm.push_back(deltaZTgTheta);
        double term1 = -deltaZTgTheta*std::pow(originalR[i], 2);
        secondOrderTerm1.push_back(term1*std::pow(chargeOverTwoRho, 2)/2.);
        double term2 = -originalR[i]*std::pow(deltaZTgTheta, 2);
        secondOrderTerm2.push_back(term2*std::pow(chargeOverTwoRho, 2)/2.);
        double term3 = -std::pow(deltaZTgTheta, 3)/3.;
        secondOrderTerm3.push_back(term3*std::pow(chargeOverTwoRho, 2)/2.);
        double secondOrderTerm = (term1+term2+term3)*std::pow(chargeOverTwoRho, 2)/2.;
        return (originalR[i] + deltaZTgTheta + secondOrderTerm);
      }
    }
  }
  return R;
}


template <class T>
double extrapolateRSecondOrder(const double & R, const double & z, const int layer, const double & tgTheta,
                               const double & chargeOverTwoRho,
                               const std::vector<int> & uniqueLayers,
                               const std::vector<double> & originalR, const T & originalZ)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i])*tgTheta;
        double term1 = -deltaZTgTheta*std::pow(originalR[i], 2);
        double term2 = -originalR[i]*std::pow(deltaZTgTheta, 2);
        double term3 = -std::pow(deltaZTgTheta, 3)/3.;
        double secondOrderTerm = (term1+term2+term3)*std::pow(chargeOverTwoRho, 2)/2.;
        return (originalR[i] + deltaZTgTheta + secondOrderTerm);
      }
    }
  }
  return R;
}


template <class T>
double extrapolateRExact(const double & R, const double & z, const int layer, const double & tgTheta,
                         const double & chargeOverTwoRho,
                         const std::vector<int> & uniqueLayers,
                         const std::vector<double> & originalR, const T & originalZ)
{
  if (layer > 10 && R > 61.) {
    // Extrapolate R from the outermost PS module in this stubs combination
    for (int i = uniqueLayers.size() - 1; i >= 0; --i) {
      if (uniqueLayers[i] < 8 || (uniqueLayers[i] > 10 && originalR[i] < 61.)) {
        double deltaZTgTheta = (z - originalZ[i])*tgTheta;
        return sin(deltaZTgTheta*chargeOverTwoRho + asin(originalR[i]*chargeOverTwoRho))/chargeOverTwoRho;
      }
    }
  }
  return R;
}


double correctPhiForNonRadialStrips(const double & phi, const double & stripPitch, const float & stripIndex,
                                    const double & extrapolatedR, const double & R, const int layer);


class CorrectPhiForNonRadialStripsLookup
{
 public:
  CorrectPhiForNonRadialStripsLookup();
//  ~CorrectPhiForNonRadialStripsLookup()
//  {
//    lookupOneOverRSquaredOutput_.close();
//  }
  double correctPhiForNonRadialStrips(const double &phi, const double &stripPitch,
//                                      const float & inputStripIndex,
                                      const double &extrapolatedR, const double &R,
                                      const double & z,
                                      const int layer) const;
 protected:
  std::vector<std::vector<double> > rns_;

  // Half length of the module perpendicular to the strip direction = strip_pitch*number_of_strips/2
  // We use number_of_strips/2+1 as it works better with the floors, roundings and precision.
  double d_;
  // Multiplicative factor used to scale the R of the stub before flooring it. It allows to gain
  // enough precision to distinguish all rings.
  int factor_;

  // Lookup table providing all the numbers necessary to perform the radial strip corrections in the 2S
  // modules of the disks. It also includes the 1/R^2, event though it is not used in this correction.
  // This number is used for the extrapolation of R we store it here for convenience.
  // The limited precision value of int(R*4) is sufficient to distinguish all the R values. Two values
  // are provided and they are valid for different z. The lookup_z map allows to select the appropriate
  // value.The shift is roughly 2pi/N, however we prefer to use the phi of the central strip in module 0.
  std::unordered_map<int, std::vector<double> > lookupR_;

  // Look-up table to determine the shift to take for the given z
  std::unordered_map<int, int> lookupZ_;

//  mutable std::ofstream lookupOneOverRSquaredOutput_;
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

protected:
  std::string name_;
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
  EstimatorSimple(const TString & inputFileName, const double & scaleFactor = 1.)
  {
    std::ifstream inputFile;
    inputFile.open(inputFileName);
    if (!inputFile) {
      // std::cout << "EstimatorSimple: Error opening " + inputFileName << std::endl;
      // throw std::runtime_error("EstimatorSimple: Error opening " + inputFileName);
      throw std::runtime_error("EstimatorSimple: Error opening " + std::string(inputFileName.Data()));
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
    parameterMean_ *= scaleFactor;

    // Read coefficients
    for (int i = 0; i < nVars; ++i) {
      inputFile >> x;
      coeff_.push_back(x*scaleFactor);
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

  /// Only needed for common call with other methods.
  void write() {};

 protected:
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
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    return stubsCombination.R(index);
  }
};


class TransformPropagateZ : public TransformBase
{
 public:
  TransformPropagateZ(const std::string & name) : TransformBase(name)
  {}
  virtual ~TransformPropagateZ() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    return stubsCombination.z(index);
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
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    std::vector<double> originalPhi;
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

    // If this is a 2S module in the disks
    R = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + estimatedChargeOverPt*DeltaR*3.8114*0.003/2. + RCube*std::pow(estimatedChargeOverPt*3.8114*0.003/2., 3)/6.);
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


class TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup(const std::string & name,
                                                                                        const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                        const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                        const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup() {}
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
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


class TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections : public TransformBase
{
 public:
  TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections(const std::string & name,
                                                                                                                      const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                                                      const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                                                      const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiSecondOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections() {}
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    // R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR + RCube*std::pow(chargeOverTwoRho, 3)/6.);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections(const std::string & name,
                                                                                                                      const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                                                      const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                                                      const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookupNoExtRInCoordinatesCorrections() {}
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    // R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    // double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup(const std::string & name,
                                                                                       const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                       const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                       const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiFirstOrderExtrapolatedRFirstOrderNonRadialStripCorrectionLookup() {}
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    // If this is a 2S module in the disks
    double extrapolatedR = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    // double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
};


class TransformCorrectedPhiFirstOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup : public TransformBase
{
 public:
  TransformCorrectedPhiFirstOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup(const std::string & name,
                                                                                        const std::string & firstOrderChargeOverPtCoefficientsFileName,
                                                                                        const std::string & firstOrderTgThetaCoefficientsFileName,
                                                                                        const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderChargeOverPtCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedPhiFirstOrderExtrapolatedRSecondOrderNonRadialStripCorrectionLookup() {}
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    phi = correctPhiForNonRadialStripsLookup_.correctPhiForNonRadialStrips(phi, 0.009, // stubsCombination.stub(index).strip(),
                                                                           extrapolatedR, R, z, layer);
    R = extrapolatedR;

    double DeltaR = R - meanRadius_[index];
    // double RCube = R*R*R;
    return (phi + chargeOverTwoRho*DeltaR);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  CorrectPhiForNonRadialStripsLookup correctPhiForNonRadialStripsLookup_;
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
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
    double chargeOverTwoRho = estimator_->estimate(originalPhi)*3.8114*0.003/2.;

    double tgTheta = 1./stubsCombination.genCotTheta();

    // If this is a 2S module in the disks
    int layer = stubsCombination.layer(index);
    double extrapolatedR = extrapolateRSecondOrder(R, z, layer, tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
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


class TransformCorrectedZFirstOrderExtrapolatedR : public TransformBase
{
 public:
  TransformCorrectedZFirstOrderExtrapolatedR(const std::string & name,
                                             const std::string & firstOrderCotThetaCoefficientsFileName,
                                             const std::string & firstOrderTgThetaCoefficientsFileName,
                                             const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName))
  {}
  virtual ~TransformCorrectedZFirstOrderExtrapolatedR() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    // double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    // If this is a 2S module in the disks
    // int layer = stubsCombination.layer(index);
    // If this is a 2S module in the disks
    double extrapolatedR = extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
    R = extrapolatedR;
    double DeltaR = R - meanRadius_[index];
    return (z - DeltaR*cotTheta);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
};


class TransformCorrectedZFirstOrderExtrapolatedRSecondOrder : public TransformBase
{
 public:
  TransformCorrectedZFirstOrderExtrapolatedRSecondOrder(const std::string & name,
                                                        const std::string & firstOrderCotThetaCoefficientsFileName,
                                                        const std::string & firstOrderTgThetaCoefficientsFileName,
                                                        const std::string & firstOrderCOverTwoRhoCoefficientsFileName,
                                                        const std::vector<double> & meanRadius) :
      TransformBase(name, firstOrderCotThetaCoefficientsFileName, meanRadius),
      estimatorTgTheta_(std::make_shared<EstimatorSimple>(firstOrderTgThetaCoefficientsFileName)),
      estimatorCOverTwoRho_(std::make_shared<EstimatorSimple>(firstOrderCOverTwoRhoCoefficientsFileName))
  {}
  virtual ~TransformCorrectedZFirstOrderExtrapolatedRSecondOrder() {}
  virtual double operator()(const StubsCombination & stubsCombination, const int index) const
  {
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    // double DeltaR = R - meanRadius_[index];
    std::vector<double> originalPhi;
    std::vector<double> originalR;
    std::vector<double> originalZ;
    for (size_t i=0; i<stubsCombination.size(); ++i) {
      originalPhi.push_back(stubsCombination.phi(i));
      originalR.push_back(stubsCombination.R(i));
      originalZ.push_back(stubsCombination.z(i));
    }
    double cotTheta = estimator_->estimate(originalR, originalZ);
    double tgTheta = estimatorTgTheta_->estimate(originalR, originalZ);
    double chargeOverTwoRho = estimatorCOverTwoRho_->estimate(originalPhi);
    // If this is a 2S module in the disks
    // int layer = stubsCombination.layer(index);
    // If this is a 2S module in the disks
    double extrapolatedR = extrapolateRSecondOrder(R, z, stubsCombination.layer(index), tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    R = extrapolatedR;
    double DeltaR = R - meanRadius_[index];
    return (z - DeltaR*cotTheta);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorTgTheta_;
  std::shared_ptr<EstimatorSimple> estimatorCOverTwoRho_;
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
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    // Correct the R for the z variation. This is needed only in the 2S modules of the disks
    // where the R resolution is too low to see the DeltaZ variation within the disks.
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      // R must increase when moving from smaller z to bigger z. If z < meanZ the sign of the added them must be positive.
      R += (meanZ_[index] - z)/stubsCombination.genCotTheta();
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
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      R = sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
    }
    double DeltaR = R - meanRadius_[index];
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
  {
    double phi = stubsCombination.phi(index);
    double R = stubsCombination.R(index);
    double z = stubsCombination.z(index);
    double chargeOverTwoRho = stubsCombination.genChargeOverPt()*3.8114*0.003/2.;
    if (stubsCombination.layer(index) > 10 && R > 61.) {
      double extrapolatedR =
          sin((z - stubsCombination.genZ0()) * chargeOverTwoRho / stubsCombination.genCotTheta()) / chargeOverTwoRho;
      // Correct the phi coordinate using the extrapolated R at the given z
      phi = correctPhiForNonRadialStrips(phi, 0.009, stubsCombination.stub(index).strip(), extrapolatedR, R,
                                         stubsCombination.layer(index));
      // Replace the R for later use
      R = extrapolatedR;
    }
    double DeltaR = R - meanRadius_[index];
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
    return extrapolateR(R, z, stubsCombination.layer(index), tgTheta, stubsCombination.layers(), originalR, originalZ);
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
    double chargeOverTwoRho = estimatorChargeOverPt_->estimate(originalPhi)*3.8114*0.003/2.;
    double extrapolatedR = extrapolateRSecondOrder(R, z, stubsCombination.layer(index), tgTheta, chargeOverTwoRho, stubsCombination.layers(), originalR, originalZ);
    return extrapolatedR;
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
    return extrapolateRExact(R, z, stubsCombination.layer(index), 1./stubsCombination.genCotTheta(), chargeOverTwoRho,
                             stubsCombination.layers(), originalR, originalZ);
  }
 private:
  std::shared_ptr<EstimatorSimple> estimatorChargeOverPt_;
};


#endif // GETVARIABLES_H
