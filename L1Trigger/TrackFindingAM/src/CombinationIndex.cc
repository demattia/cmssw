//
// Created by Marco De Mattia on 7/11/15.
//

#include "L1Trigger/TrackFindingAM/interface/CombinationIndex.h"
#include <iostream>


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius,
                               const int regionsNumber)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);

  // Set bits to determine the type of modules in the disks (2S vs PS) and upper or lower 2S radius.
  // Their positions in the bitset are the disk index + 5, since there
  // are 5 disks per side and they never appear together.
  for (unsigned int i=0; i<layers.size(); ++i) {
    setLayerRadiusBits(layers[i], radius[i], bits, regionsNumber);
  }
  return bits.to_ulong();
}


unsigned long combinationIndex(const std::vector<Stub> & stubs, const int regionsNumber)
{
  std::bitset<32> bits;
  for (const Stub & stub : stubs) {
    int layer = stub.layer();
    bits.set(layer, 1);
    // Set bits to determine the type of modules in the disks (2S vs PS).
    // Their positions in the bitset are the disk index + 10, since there are 10 disks in total.
    setLayerRadiusBits(layer, stub.R(), bits, regionsNumber);
  }
  return bits.to_ulong();
}
