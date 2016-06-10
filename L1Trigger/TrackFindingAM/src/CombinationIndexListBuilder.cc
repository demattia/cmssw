//
// Created by Marco De Mattia on 7/21/15.
//

#include "../interface/CombinationIndexListBuilder.h"

void CombinationIndexListBuilder::allCombinationIndexes(const std::vector<int> &layers,
                                                        const std::vector<double> &radius,
                                                        std::vector<unsigned long> &combinationIndexList,
                                                        const bool fiveOutOfSix)
{
  combinationIndexList.push_back(combinationIndex(layers, radius));
  // std::cout << "6/6 combination index = " << combinationIndexList.back() << std::endl;
  // Generate names for the six 5/6 cases in this region
  if (fiveOutOfSix) {
    for (int i = 0; i < 6; ++i) {
      std::vector<int> indexes(combinationsGenerator_.combination(i, 6));
      std::vector<int> layers_5_6;
      std::vector<double> radius_5_6;
      for (auto index : indexes) {
        layers_5_6.push_back(layers.at(index));
        radius_5_6.push_back(radius.at(index));
      }
      combinationIndexList.push_back(combinationIndex(layers_5_6, radius_5_6));
    }
  }
}


void CombinationIndexListBuilder::fillDefaultIndexList(std::vector<unsigned long> & combinationIndexList,
                                                       const bool fiveOutOfSix)
{
  // Region 1
  std::vector<int> layers = {5, 6, 7, 8, 9, 10};
  std::vector<double> radius = {0., 0., 0., 0., 0., 0.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 2
  layers = std::vector<int>{5, 6, 7, 8, 9, 11};
  radius = std::vector<double>{0., 0., 0., 0., 0., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 3
  layers = std::vector<int>{5, 6, 7, 8, 11, 12};
  radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 4
  layers = std::vector<int>{5, 6, 7, 11, 12, 13};
  radius = std::vector<double>{0., 0., 0., 100., 100., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 5
  layers = std::vector<int>{5, 6, 11, 12, 13, 14};
  radius = std::vector<double>{0., 0., 0., 100., 100., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 6
  layers = std::vector<int>{5, 6, 11, 12, 14, 15};
  radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 7
  layers = std::vector<int>{5, 11, 12, 13, 14, 15};
  radius = std::vector<double>{0., 0., 0., 0., 100., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 8
  radius = std::vector<double>{0., 0., 0., 0., 0., 100.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);
  // Region 9
  radius = std::vector<double>{0., 0., 0., 0., 0., 0.};
  allCombinationIndexes(layers, radius, combinationIndexList, fiveOutOfSix);

  // Remove duplicates
  std::sort(combinationIndexList.begin(), combinationIndexList.end());
  combinationIndexList.erase(std::unique(combinationIndexList.begin(), combinationIndexList.end()),
                             combinationIndexList.end());
}
