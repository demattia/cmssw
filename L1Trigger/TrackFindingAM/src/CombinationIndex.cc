//
// Created by Marco De Mattia on 7/11/15.
//

#include "../interface/CombinationIndex.h"
#include <iostream>

unsigned long combinationIndex(const std::vector<int> & layers, const int region)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);
  int lastDisk = 10;
  if (region == 9) lastDisk = 15;
  if (region == 8) lastDisk = 14;
  if (region == 7) lastDisk = 13;
  if (region == 6) lastDisk = 12;
  if (region == 5) lastDisk = 11;

  for (int disk=11; disk<=lastDisk; ++disk) {
    // Only set the bit of the disk is there
    if (bits[disk]) bits.set(disk+10, 1);
  }

  return bits.to_ulong();
}


unsigned long combinationIndex(const std::vector<int> & layers, const std::vector<double> & radius)
{
  std::bitset<32> bits;
  // Set the bits for the layers
  combinationIndex(layers, bits);

  // Set bits to determine the type of modules in the disks (2S vs PS).
  // Their positions in the bitset are the disk index + 10, since there are 10 disks in total.
  for (unsigned int i=0; i<layers.size(); ++i) {
    if (layers[i] > 10 && radius[i] < 61.) bits.set(layers[i]+10, 1);
  }
  return bits.to_ulong();
}


unsigned long combinationIndex(const std::vector<Stub> & stubs)
{
  std::bitset<32> bits;
  for (const Stub & stub : stubs) {
    int layer = stub.layer();
    bits.set(layer, 1);
    // Set bits to determine the type of modules in the disks (2S vs PS).
    // Their positions in the bitset are the disk index + 10, since there are 10 disks in total.
    if (layer > 10 && stub.R() < 61.) bits.set(layer+10, 1);
  }
  return bits.to_ulong();
}
