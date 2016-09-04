//
// Created by Marco De Mattia on 2/28/16.
//

#ifndef LINEARIZEDTRACKFIT_ROAD_H
#define LINEARIZEDTRACKFIT_ROAD_H

#include <vector>
#include <map>
#include "L1Trigger/TrackFindingAM/interface/Stub.h"
#include "L1Trigger/TrackFindingAM/interface/Hit.h"

/**
 * A Road contains a vector of 6 vectors, one per layer. Each vector contains the stubs in that layer.
 */

class Road
{
 public:
  // Road(const std::shared_ptr<RoadsTree> & tree, const size_t & roadIndex);
  // Road(const std::shared_ptr<RoadsTree> & tree, const std::map<unsigned long, unsigned long> & stubsIndexes);
  Road(const std::vector<Hit*> & stubs);

  std::vector<std::vector<Stub> >::const_iterator beginLayer() const { return stubs_.begin(); }
  std::vector<std::vector<Stub> >::const_iterator endLayer() const { return stubs_.end(); }
  std::size_t stubsNum(const std::size_t & layer) const { return stubs_.at(layer).size(); }
  Stub getStub(const std::size_t & layer, const int count) const { return stubs_.at(layer).at(count); }

 private:
  int maxStubsPerLayer_;
  std::vector<std::vector<Stub> > stubs_;
};

#endif //LINEARIZEDTRACKFIT_ROAD_H
