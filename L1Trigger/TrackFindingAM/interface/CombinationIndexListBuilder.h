//
// Created by Marco De Mattia on 7/21/15.
//

#ifndef REMOTEPROJECTS_COMBINATIONINDEXLISTBUILDER_H
#define REMOTEPROJECTS_COMBINATIONINDEXLISTBUILDER_H

#include <algorithm>
#include "CombinationIndex.h"
#include "CombinationsGenerator.h"

class CombinationIndexListBuilder
{
 public:
  /// Returns the list of combination indexes for the 6/6 and 5/6 of a given combination of 6 layers/disks and radii
  void allCombinationIndexes(const std::vector<int> &layers, const std::vector<double> &radius,
                             std::vector<unsigned long> &combinationIndexList, const bool fiveOutOfSix);

  void fillDefaultIndexList(std::vector<unsigned long> & combinationIndexList, const bool fiveOutOfSix);

 private:
  CombinationsGenerator combinationsGenerator_;
};



#endif //REMOTEPROJECTS_COMBINATIONINDEXLISTBUILDER_H
