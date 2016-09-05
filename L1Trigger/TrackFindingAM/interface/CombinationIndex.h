//
// Created by Marco De Mattia on 7/11/15.
//

#ifndef REMOTEPROJECTS_COMBINATIONINDEX_H
#define REMOTEPROJECTS_COMBINATIONINDEX_H

#include <vector>
#include <bitset>
#include "L1Trigger/TrackFindingAM/interface/Stub.h"


inline void combinationIndex(const std::vector<int> & layers, std::bitset<32> & bits) { for (auto l : layers) bits.set(l, 1); }


unsigned long combinationIndex(const std::vector<int> & layers, const int region);


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius, const int regionsNumber);


unsigned long combinationIndex(const std::vector<Stub> & stubs, const int regionsNumber);


void allCombinationIndexes(const std::vector<int> & layers, const std::vector<double> & radius,
                           std::vector<unsigned long> & combinationIndexList, const bool fiveOutOfSix);


template <class T>
inline void setLayerRadiusBits(const int layer, const double & radius, T & bits, const int regionsNumber)
{
  if (layer > 10) {
    if (radius < 61.) {
      bits.set(layer + 5, 1);
    }
    // Split the 2S modules part of the disks
    else if (regionsNumber == 14 &&
        (((layer == 11 || layer == 12) && radius < 82.5) ||
        ((layer == 13) && radius < 82.5) ||
        ((layer == 14) && radius < 77.) ||
        ((layer == 15) && radius < 77.))) {
      bits.set(layer + 10, 1);
    }
//    else if (((layer == 11 || layer == 12) && radius < 82.5) ||
//             ((layer == 13) && radius < 77.) ||
//             ((layer == 14) && radius < 72.) ||
//             ((layer == 15) && radius < 77.)) {
//      bits.set(layer + 10, 1);
//    }
  }
}



#endif //REMOTEPROJECTS_COMBINATIONINDEX_H
