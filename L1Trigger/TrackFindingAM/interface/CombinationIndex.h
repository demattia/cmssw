//
// Created by Marco De Mattia on 7/11/15.
//

#ifndef REMOTEPROJECTS_COMBINATIONINDEX_H
#define REMOTEPROJECTS_COMBINATIONINDEX_H

#include <vector>
#include <bitset>
#include "Stub.h"


inline void combinationIndex(const std::vector<int> & layers, std::bitset<32> & bits) { for (auto l : layers) bits.set(l, 1); }


unsigned long combinationIndex(const std::vector<int> & layers, const int region);


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius);


unsigned long combinationIndex(const std::vector<Stub> & stubs);


void allCombinationIndexes(const std::vector<int> & layers, const std::vector<double> & radius,
                           std::vector<unsigned long> & combinationIndexList, const bool fiveOutOfSix);

#endif //REMOTEPROJECTS_COMBINATIONINDEX_H
