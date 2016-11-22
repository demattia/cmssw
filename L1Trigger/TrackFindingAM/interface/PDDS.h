//
// Created by Denis Rathjens on 11/17/16.
//

#ifndef _PDDS_H
#define _PDDS_H

#include <vector>
#include "L1Trigger/TrackFindingAM/interface/StubsCombination.h"
#include "L1Trigger/TrackFindingAM/interface/Road.h"

/**
 * Build all possible combinations of 6/6 and 5/5 pairs of stubs. Also build all possible combinations of 5 stubs out of
 * sets of 6 layers/disks without repetitions.
 */

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct PairAssignment //structure for saving pairs for the pairwise Delta Delta S combination cleanup
{
  unsigned Stub1Address;
  unsigned Stub2Address;
  bool SurvivingPair;
  PairAssignment(unsigned Stub1Address_, unsigned Stub2Address_, bool SurvivingPair_){
    Stub1Address=Stub1Address_;
    Stub2Address=Stub2Address_;
    SurvivingPair=SurvivingPair_;
  }
};

struct IndexInfo //contains compact information about the target kind of combination
{
  unsigned int target;
  int layers[6];
};

class PDDS 
{
 public:
  PDDS();
  std::vector<StubsCombination> combine(Road & road);
  void SetFive(const bool SetFive) { SetFive_=SetFive;}

 private:

  void initialize(std::vector<unsigned> layer, std::vector<float> radii);
  unsigned int Index;
  unsigned int CombinationKind(std::vector<unsigned> layer_, std::vector<float> radii_);
  unsigned int IndexInterpreter(unsigned int Index_);
  std::vector<PairAssignment> Stage1Pair(std::vector<std::pair<unsigned, float> > first, std::vector<std::pair<unsigned, float> > second, float cut);
  std::vector<std::vector<unsigned> > CombinationBuilder(std::vector<PairAssignment> Layer10, std::vector<PairAssignment> Layer32, std::vector<PairAssignment> Layer54);
  float SebCoarsener(float DS_, int layer_, float z_, unsigned ring_);
  std::vector<std::vector<unsigned> > CombinationAll(std::vector<Stub> stubs_);
  IndexInfo IndexDefiner(unsigned int target_);
  float cuts_[3];
  bool SetFive_;
  unsigned int DummyAllowance_;
  bool verbose_;

};

#endif //PDDS_H
