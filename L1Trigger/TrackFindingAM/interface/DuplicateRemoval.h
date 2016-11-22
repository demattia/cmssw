//
// Created by Marco De Mattia on 9/2/16.
//

#ifndef SEBSFRAMEWORK_DUPLICATEREMOVAL_H
#define SEBSFRAMEWORK_DUPLICATEREMOVAL_H

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include <vector>
#include <algorithm>


namespace slhcl1tt {

  class DuplicateRemoval
  {
   public:
    // Constructor
    DuplicateRemoval() {}

    // Destructor
    ~DuplicateRemoval() {}

    // Return flags categorizing as duplicate (1) or not (0)
    // void CheckTracks(std::vector<TTTrack2> &full_am_track_list, int dupRm);

    // Return flags categorizing as duplicate (1) or not (0)
    std::vector<TTTrack<Ref_PixelDigi_> > CheckTracks(const std::vector<TTTrack<Ref_PixelDigi_> >& full_am_track_list, const int dupRm = 2);
  };

}

#endif //SEBSFRAMEWORK_DUPLICATEREMOVAL_H
