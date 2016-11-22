//
// Created by Marco De Mattia on 2/28/16.
//

#ifndef LINEARIZEDTRACKFIT_SIMPLECOMBINATIONBUILDER_H
#define LINEARIZEDTRACKFIT_SIMPLECOMBINATIONBUILDER_H

#include "L1Trigger/TrackFindingAM/interface/CombinationBuilderBase.h"

/**
 * Builds all possible combinations from the stubs in the 6 layers in a road.
 * The current hardware design allows at most 4 stubs per layer. If there are more than 4 only the last 4 are retained.
 */

class SimpleCombinationBuilder : public CombinationBuilderBase
{
 public:
  SimpleCombinationBuilder();
  virtual void initialize(const Road & road);
  virtual StubsCombination nextCombination();

 private:

  struct CombUnit
  {
    CombUnit(const int range) :
        counter_(range), range_(range)
    {}
    int comb() const { return counter_; }
    void setRange(const int range)
    {
      range_ = range;
      counter_ = range;
    }
    int decrement(const int value)
    {
      if (value == 0) return 0;
      int carryOver = counter_ == 0 ? 1 : 0;
      counter_ -= value;
      if (carryOver != 0) counter_ = range_;
      return carryOver;
    }
   private:
    int counter_;
    int range_;
  };

  std::vector<CombUnit> combUnits_;
};

#endif //LINEARIZEDTRACKFIT_SIMPLECOMBINATIONBUILDER_H
