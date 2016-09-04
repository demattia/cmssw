//
// Created by Marco De Mattia on 7/5/16.
//

#ifndef LINEARIZEDTRACKFIT_ADVANCEDCOMBINATIONBUILDER_H
#define LINEARIZEDTRACKFIT_ADVANCEDCOMBINATIONBUILDER_H

#include "L1Trigger/TrackFindingAM/interface/CombinationBuilderBase.h"

/**
 * Build all possible combinations of 6/6 and 5/5 stubs. Also build all possible combinations of 5 stubs out of
 * sets of 6 layers/disks without repetitions.
 */

class AdvancedCombinationBuilder : public CombinationBuilderBase
{
 public:
  AdvancedCombinationBuilder();
  virtual void initialize(const Road & road);
  virtual StubsCombination nextCombination();

 private:

  struct CombUnit
  {
    CombUnit(const unsigned int index, const int range) :
        index_(index), counter_(range), range_(range), store0_(0), store1_(0)
    {}
    int comb() const { return counter_; }

    void setRange(const int range)
    {
      range_ = range;
      counter_ = range;
    }
    int decrement(const int value, std::vector<bool> & zeroEnable)
    {
      if (value == 0) return 0;
      store0_ = (counter_ == 0);
      store1_ = (counter_ == 1);
      counter_ -= value;
      bool carryOver = (zeroEnable.at(index_+1) && store1_) || store0_;
      if (carryOver) counter_ = range_;
      zeroEnable[index_] = store0_;
      return (carryOver ? 1 : 0);
    }
    void propagateZeroEnable(std::vector<bool> & zeroEnable)
    {
      zeroEnable[index_] = (counter_ == 0) || zeroEnable[index_+1];
    }

   private:
    unsigned int index_;
    int counter_;
    int range_;
    bool store0_;
    bool store1_;
  };

//  int isZero(const size_t & value) const;
  bool rangesOr() const;

  std::vector<CombUnit> combUnits_;
  std::vector<bool> zeroEnable_;
};

#endif //LINEARIZEDTRACKFIT_ADVANCEDCOMBINATIONBUILDER_H
