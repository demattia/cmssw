//
// Created by Marco De Mattia on 7/5/16.
//

#include "L1Trigger/TrackFindingAM/interface/CombinationBuilderBase.h"


CombinationBuilderBase::CombinationBuilderBase() :
    road_(nullptr), totalCombinations_(1)
{}

void CombinationBuilderBase::fillGenInfo(StubsCombination & stubsCombination)
{
  // bool first = true;
  // int tpId = -1;
  // for (auto stub : stubsCombination) {
  //   if (first) {
  //     tpId = road_->tpId(stub.stubRef());
  //     first = false;
  //   }
  //   else if (tpId != road_->tpId(stub.stubRef())) {
  //     tpId = -1;
  //     break;
  //   }
  // }
  // if (tpId != -1) {
  //   unsigned int stubRef = stubsCombination.begin()->stubRef();
  //   stubsCombination.setGenTrack(road_->genChargeOverPt(stubRef), road_->genPhi0(stubRef), road_->genD0(stubRef),
  //                                road_->genCotTheta(stubRef), road_->genZ0(stubRef));
  // }
}
