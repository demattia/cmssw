//
// Created by Marco De Mattia on 7/25/16.
//

#ifndef SEBSFRAMEWORK_PARAMETERDUPLICATEREMOVAL_H
#define SEBSFRAMEWORK_PARAMETERDUPLICATEREMOVAL_H

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include <vector>
#include <TVector2.h>
#include <assert.h>
#include <cmath>

namespace slhcl1tt{
  class ParameterDuplicateRemoval
  {
   public:
    ParameterDuplicateRemoval() {}
    ~ParameterDuplicateRemoval() {}

    std::vector<TTTrack<Ref_PixelDigi_> > ReduceTracks(const std::vector<TTTrack<Ref_PixelDigi_> > & Tracks);
  };

//structure containing necessary original Track information
  struct TrackCloud{
    std::vector<unsigned> TrackIterators;
    unsigned  BestTrack;
    bool BestRoadCategory; //false = 5/x, true = 6/6
    float BestChi2;
    float BestEta;
    float BestPhi;
    float BestPt;
    float BestCharge;

    template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
    }

    bool Add(unsigned TrackIterator_, float Chi2_, float Eta_, float Phi_, float Pt_, float Charge_){ //test whether a track fits into a cloud or not
      //define binned pT barriers
      std::vector<std::pair<float,float> > PtEdgeAndDelta;
      PtEdgeAndDelta.push_back(std::make_pair(3.,5.));
      PtEdgeAndDelta.push_back(std::make_pair(5.,8.));//3.4));
      PtEdgeAndDelta.push_back(std::make_pair(8.,3.5));//2.8));
      PtEdgeAndDelta.push_back(std::make_pair(13.,15.0));//7.8));
      PtEdgeAndDelta.push_back(std::make_pair(22.,17.9));//18.1));
      PtEdgeAndDelta.push_back(std::make_pair(37.,60.));//33.5));
      PtEdgeAndDelta.push_back(std::make_pair(60.,51.));//52.4));
      PtEdgeAndDelta.push_back(std::make_pair(99.,164.));
      PtEdgeAndDelta.push_back(std::make_pair(164.,270.));
      PtEdgeAndDelta.push_back(std::make_pair(270.,445.));
      PtEdgeAndDelta.push_back(std::make_pair(445.,734.));//447.9));
      PtEdgeAndDelta.push_back(std::make_pair(734.,1210.));
      PtEdgeAndDelta.push_back(std::make_pair(1210.,1995.));
      PtEdgeAndDelta.push_back(std::make_pair(1995.,3290.));
      Charge_=sgn(Charge_);

      //find the appropriate pT difference for a given maximum pT between track and cloud
      float DeltaPt=10.;
      for(int i=PtEdgeAndDelta.size()-1; i>=0; --i){
        if(BestPt > PtEdgeAndDelta[i].first || Pt_ > PtEdgeAndDelta[i].first){
          DeltaPt=PtEdgeAndDelta[i].second;
          break;
        }
      }
      //std::cout<<DeltaPt<<" for "<<BestPt<<"-"<<Pt_<<std::endl;

      if(TrackIterators.size()==0){
        TrackIterators.push_back(TrackIterator_);
        BestChi2=Chi2_;
        BestEta=Eta_;
        BestPhi=Phi_;
        BestPt=Pt_;
        BestCharge=Charge_;
        BestTrack=TrackIterator_;
        return true;
      }//previously 0.0601 && 0.00576
      else if(fabs(Eta_-BestEta)<0.0301 && fabs(TVector2::Phi_mpi_pi(Phi_-BestPhi))<0.00306 && fabs(BestPt-Pt_)<DeltaPt*0.9 && BestCharge==Charge_){ //explicit cloud inclusion test
        if(Chi2_<BestChi2){
          BestChi2=Chi2_;
          BestEta=Eta_;
          BestPhi=Phi_;
          BestPt=Pt_;
          BestCharge=Charge_;
          BestTrack=TrackIterator_;
        }
        TrackIterators.push_back(TrackIterator_);
        return true;
      }
      else return false;
    }
  };
}


#endif //SEBSFRAMEWORK_PARAMETERDUPLICATEREMOVAL_H
