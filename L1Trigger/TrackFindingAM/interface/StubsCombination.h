//
// Created by Marco De Mattia on 7/30/15.
//

#ifndef REMOTEPROJECTS_STUBSCOMBINATION_H
#define REMOTEPROJECTS_STUBSCOMBINATION_H

#include <vector>
#include <cmath>
#include "Stub.h"
#include "CombinationIndex.h"

class StubsCombination
{
 public:
  StubsCombination() :
      genChargeOverPt_(0.), genPhi0_(0.), genD0_(0.), genCotTheta_(0.), genZ0_(0.), combinationIndex_(0)
  {}
  void clear();
  void pushStub(const double & phi, const double & R, const double & z, const int layer, const float strip);
  void setGenTrack(const double & genChargeOverPt, const double & genPhi0, const double & genD0,
                   const double & genCotTheta, const double & genZ0);
  void build(const StubsCombination & stubsCombination, const std::vector<int> & combination);
  double genChargeOverPt() const { return genChargeOverPt_; }
  double genPhi0() const { return genPhi0_; }
  double genD0() const { return genD0_; }
  double genCotTheta() const { return genCotTheta_; }
  double genZ0() const { return genZ0_; }
  const Stub & stub(const int i) const { return stubs_.at(i); }
  void setCombinationIndex() { combinationIndex_ = combinationIndex(stubs_); }
  unsigned long getCombinationIndex() const { return combinationIndex_; }
  size_t size() const { return stubs_.size(); }
  double phi(const int index) const { return stubs_.at(index).phi(); }
  double R(const int index) const { return stubs_.at(index).R(); }
  double z(const int index) const { return stubs_.at(index).z(); }
  int layer(const int index) const { return stubs_.at(index).layer(); }
  void setPhi(const int index, const double & phi) { stubs_.at(index).setPhi(phi); }
  void setR(const int index, const double & R) { stubs_.at(index).setR(R); }
  void setZ(const int index, const double & z) { stubs_.at(index).setZ(z); }
  std::vector<double> phiVector() const;
  std::vector<double> RVector() const;
  std::vector<double> zVector() const;
  std::vector<double> variables() const;
  std::vector<int> layers() const;
  std::vector<Stub>::const_iterator begin() const { return stubs_.begin(); }
  std::vector<Stub>::const_iterator end() const { return stubs_.end(); }
  double genTrackDistanceTransverse(const int index) const;
  double genTrackDistanceTransverseFromZ(const int index) const;
  double genTrackDistanceLongitudinal(const int index) const;
  double genTrackDistanceLongitudinalR(const int index) const;

 private:
  double genChargeOverPt_;
  double genPhi0_;
  double genD0_;
  double genCotTheta_;
  double genZ0_;
  std::vector<Stub> stubs_;
  unsigned long combinationIndex_;
};


#endif //REMOTEPROJECTS_STUBSCOMBINATION_H
