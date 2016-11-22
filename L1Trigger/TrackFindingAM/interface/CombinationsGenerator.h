//
// Created by Marco De Mattia on 7/9/15.
//

#ifndef REMOTEPROJECTS_COMBINATIONSGENERATOR_H
#define REMOTEPROJECTS_COMBINATIONSGENERATOR_H

#include <vector>
#include <map>
#include <iostream>
#include <numeric>

class CombinationsGenerator
{
 public:
  CombinationsGenerator()
  {
    // Generate all combinations for 6/8 and for 5/8 when there are 8 and store them.
    std::vector<int> indexes(8);
    std::iota(indexes.begin(), indexes.end(), 0);
    combinations(indexes, 7);
    combinations(indexes, 6);
    combinations(indexes, 5);

    // Generate all combinations for 6/7 and for 5/7 when there are 7 and store them.
    indexes = std::vector<int>(7);
    std::iota(indexes.begin(), indexes.end(), 0);
    // combinations(indexes, 7);
    combinations(indexes, 6);
    combinations(indexes, 5);

    // Generate all default combinations for 5/6 when there are 6 and store them.
    indexes = std::vector<int>(6);
    std::iota(indexes.begin(), indexes.end(), 0);
    // combinations(indexes, 6);
    combinations(indexes, 5);
  }

  int combinationsSize(const int combinationsIndex)
  {
    return combs_[combinationsIndex].size();
  }

  std::vector<int> combination(const int index, const int combinationsIndex)
  {
    return combs_[combinationsIndex][index];
  }

 private:
  std::vector<std::vector<int> > combinations(const std::vector<int> & indexes, const int combLength)
  {
    combLength_ = combLength;
    // The index of this set of combinations is the length in the input vector.
    combinationsIndex_ = indexes.size();
    auto c = combs_.find(combinationsIndex_);
    if (c == combs_.end()) combs_.insert(std::make_pair(combinationsIndex_, std::vector<std::vector<int> >()));
    if (combLength > combinationsIndex_) {
      std::cout << "Error: combination length = " << combLength << " larger than number of layers/disks = " << combinationsIndex_ << std::endl;
      throw;
    }
    if (combLength == combinationsIndex_) {
      combs_[combinationsIndex_].push_back(indexes);
      return combs_[combinationsIndex_];
    }
    combination_ = std::vector<int>(combLength_);
    generateCombinations(indexes, 0, 0);
    return combs_[combinationsIndex_];
  }

  // Note: it is fine to pass the values for a small number of recursions. If the number grows the temporary copies
  // can require too much memory. It is better in that case to increase and decrease values in data members.
  void generateCombinations(const std::vector<int> & indexes, const int start, const int positionInCombination)
  {
    if (positionInCombination == combLength_) {
      combs_[combinationsIndex_].push_back(combination_);
    }
    else {
      for (int i=start; i<=combinationsIndex_ - combLength_ + positionInCombination; ++i) {
        combination_[positionInCombination] = indexes[i];
        generateCombinations(indexes, i+1, positionInCombination+1);
      }
    }
  }

  std::map<int, std::vector<std::vector<int> > > combs_;
  int combLength_;
  int combinationsIndex_;
  std::vector<int> combination_;
};

#endif //REMOTEPROJECTS_COMBINATIONSGENERATOR_H
