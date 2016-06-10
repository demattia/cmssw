//
// Created by Marco De Mattia on 7/30/15.
//

#include "../interface/StubsCombination.h"

void StubsCombination::clear()
{
  genChargeOverPt_ = 0.;
  genPhi0_ = 0.;
  genD0_ = 0.;
  genCotTheta_ = 0.;
  genZ0_ = 0.;
  stubs_.clear();
  combinationIndex_ = 0;
}

void StubsCombination::pushStub(const double & phi, const double & R, const double & z,
                                const int layer, const float strip)
{
  stubs_.push_back(Stub(phi, R, z, layer, strip));
}


void StubsCombination::setGenTrack(const double & genChargeOverPt, const double & genPhi0, const double & genD0,
                                   const double & genCotTheta, const double & genZ0)
{
  genChargeOverPt_ = genChargeOverPt;
  genPhi0_ = genPhi0;
  genD0_ = genD0;
  genCotTheta_ = genCotTheta;
  genZ0_ = genZ0;
}


void StubsCombination::build(const StubsCombination & stubsCombination, const std::vector<int> & combination)
{
  setGenTrack(stubsCombination.genChargeOverPt(), stubsCombination.genPhi0(), stubsCombination.genD0(),
              stubsCombination.genCotTheta(), stubsCombination.genZ0());
  stubs_.clear();
  for (auto i : combination) stubs_.push_back(stubsCombination.stub(i));
  setCombinationIndex();
}


std::vector<double> StubsCombination::phiVector() const
{
  std::vector<double> phiVec;
  for (const Stub & s : stubs_) phiVec.push_back(s.phi());
  return phiVec;
}


std::vector<double> StubsCombination::RVector() const
{
  std::vector<double> RVec;
  for (const Stub & s : stubs_) RVec.push_back(s.R());
  return RVec;
}


std::vector<double> StubsCombination::zVector() const
{
  std::vector<double> zVec;
  for (const Stub & s : stubs_) zVec.push_back(s.z());
  return zVec;
}


std::vector<int> StubsCombination::layers() const
{
  std::vector<int> layersVec;
  for (const Stub & s : stubs_) layersVec.push_back(s.layer());
  return layersVec;
}


std::vector<double> StubsCombination::variables() const
{
  std::vector<double> variablesVec;
  for (const Stub & s : stubs_) {
    variablesVec.push_back(s.phi());
    variablesVec.push_back(s.R());
    variablesVec.push_back(s.z());
  }
  return variablesVec;
}


double StubsCombination::genTrackDistanceTransverse(const int index) const
{
  const Stub & stub = stubs_.at(index);
  double rho = genChargeOverPt_ != 0 ? (1./genChargeOverPt_)/(3.8114*0.003) : 10000.;
  double R = stub.R();
  double phiGen = genPhi0_ - asin((genD0_*genD0_ + 2*genD0_*rho + R*R)/(2*R*(rho+genD0_)));
  double deltaPhi = (stub.phi() - phiGen);
  if (deltaPhi > M_PI) deltaPhi -= M_PI;
  else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  return deltaPhi;
}


double StubsCombination::genTrackDistanceTransverseFromZ(const int index) const
{
  const Stub & stub = stubs_.at(index);
  double rho = genChargeOverPt_ != 0 ? (1./genChargeOverPt_)/(3.8114*0.003) : 10000.;
  double phiGen = genPhi0_ - (stub.z()-genZ0_)/(2*rho*genCotTheta_);
  double deltaPhi = (stub.phi() - phiGen);
  if (deltaPhi > M_PI) deltaPhi -= M_PI;
  else if (deltaPhi < -M_PI) deltaPhi += M_PI;
  return deltaPhi;
}


double StubsCombination::genTrackDistanceLongitudinal(const int index) const
{
  const Stub & stub = stubs_.at(index);
  double rho = genChargeOverPt_ != 0 ? (1./genChargeOverPt_)/(3.8114*0.003) : 10000.;
  double R = stub.R();
  double zGen = genZ0_ + 2*rho*genCotTheta_*asin((genD0_*genD0_ + 2*genD0_*rho + R*R)/(2*R*(rho+genD0_)));
  return (stub.z() - zGen);
}


double StubsCombination::genTrackDistanceLongitudinalR(const int index) const
{
  const Stub & stub = stubs_.at(index);
  double rho = genChargeOverPt_ != 0 ? (1./genChargeOverPt_)/(3.8114*0.003) : 10000.;
  double z = stub.z();
  double RGen = 2*rho*sin((z - genZ0_)/(2*rho*genCotTheta_));
  return (stub.R() - RGen);
}
