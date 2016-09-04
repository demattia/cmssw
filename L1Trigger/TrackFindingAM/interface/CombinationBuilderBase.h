//
// Created by Marco De Mattia on 7/5/16.
//

#ifndef LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H
#define LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H

#include "L1Trigger/TrackFindingAM/interface/StubsCombination.h"
#include "L1Trigger/TrackFindingAM/interface/Road.h"

/**
 * Defines the interface of a combination builder and implements some of the common data members and methods.
 */

class CombinationBuilderBase
{
 public:
  CombinationBuilderBase();
  virtual void initialize(const Road & road) = 0;
  virtual StubsCombination nextCombination() = 0;
  int totalCombinations() const { return totalCombinations_; }

 protected:
  void fillGenInfo(StubsCombination & stubsCombination);
  std::vector<int> carryOver_;
  const Road * road_;
  int totalCombinations_;
};

#endif //LINEARIZEDTRACKFIT_COMBINATIONBUILDERBASE_H_H
