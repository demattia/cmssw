//
// Created by Marco De Mattia on 7/5/16.
//

#include "L1Trigger/TrackFindingAM/interface/AdvancedCombinationBuilder.h"


AdvancedCombinationBuilder::AdvancedCombinationBuilder() :
    zeroEnable_(std::vector<bool>(7, false))
{
  for (unsigned int layer = 0; layer < 6; ++layer) {
    combUnits_.push_back(CombUnit(layer, 0));
  }
}


void AdvancedCombinationBuilder::initialize(const Road & road)
{
  // Set the ranges for the counters to this combination
  totalCombinations_ = 1;
  int emptyLayers = 0;
  size_t layer = 0;
  for (auto l = road.beginLayer(); l != road.endLayer(); ++l, ++layer) {
    int range = l->size();
    combUnits_.at(layer).setRange(range);
    totalCombinations_ *= range;
    if (range == 0) ++emptyLayers;
  }
  // Count also all possible 5/6 without repetitions
  layer = 0;
  for (auto l1 = road.beginLayer(); l1 != road.endLayer(); ++l1, ++layer) {
    int combs = 1;
    size_t layer2 = 0;
    for (auto l2 = road.beginLayer(); l2 != road.endLayer(); ++l2, ++layer2) {
      if (layer != layer2) combs *= l2->size();
    }
    totalCombinations_ += combs;
  }
  // If there is more than one missing stub compute the combinations
  // as in the simple combination builder
  if (emptyLayers > 1) {
    int combs = 1;
    for (auto l = road.beginLayer(); l != road.endLayer(); ++l) {
      if (l->size() > 0) {
        combs *= l->size();
      }
    }
    // We do equal since there are no 6/6 and 5/6 (or 5/5) combinations in this case
    totalCombinations_ = combs;
    // std::cout << "TotalCombinations for a road with less than 5 layers = " << totalCombinations_ << std::endl;
  }
  carryOver_ = std::vector<int>(7, 0);
  carryOver_[0] = 1;

  road_ = &road;

  // Initialize the zero-enable. This must be after the pointer to the road is stored.
  zeroEnable_[6] = rangesOr();
  // std::cout << "Empty layers = " << emptyLayers << std::endl;
}


StubsCombination AdvancedCombinationBuilder::nextCombination()
{
  std::vector<int> stubIndexes;
  StubsCombination stubsCombination;
  if (carryOver_[6] == 0) {
    // The zero-enable propagates backwards
    for (int layer = 5; layer >= 0; --layer) {
      combUnits_[layer].propagateZeroEnable(zeroEnable_);
    }
    for (size_t layer = 0; layer < 6; ++layer) {
      if (road_->stubsNum(layer) > 0 && combUnits_[layer].comb() != 0) {
        // The -1 is because the comb() returns 0 for a dummy stub and starts counting real stubs from 1
        stubsCombination.pushStub(road_->getStub(layer, combUnits_[layer].comb()-1));
        stubIndexes.push_back(combUnits_[layer].comb());
      }
      else {
        stubIndexes.push_back(0);
      }
    }
    for (size_t layer = 0; layer < 6; ++layer) {
      carryOver_[layer + 1] = combUnits_[layer].decrement(carryOver_[layer], zeroEnable_);
    }
  }

//  if (totalCombinations_ > 50) {
//    std::cout << "stubIndexes = ";
//    for (auto s : stubIndexes) std::cout << s << " ";
//    std::cout << std::endl;
//  }

  // Fill generator-level information
  fillGenInfo(stubsCombination);

  return stubsCombination;
}


//int AdvancedCombinationBuilder::isZero(const size_t & value) const
//{
//  // Returns 0 if the number is zero or 1 otherwise.
//  // In an FPGA this could be implemented for instance as an OR of the bits of the number.
//  if (value == 0) return 0;
//  return 1;
//}


bool AdvancedCombinationBuilder::rangesOr() const
{
  // Returns 1 if at least one range is zero or 0 otherwise.
  for (auto l = road_->beginLayer(); l != road_->endLayer(); ++l) {
    if (l->size() == 0) return true;
  }
  return false;
}

