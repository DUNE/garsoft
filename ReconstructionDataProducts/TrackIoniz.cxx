//
//  TrackIoniz.cxx
//  garsoft-mrb
//
//  Created by Leo Bellantoni on 22 Feb, 2018

#include "ReconstructionDataProducts/TrackIoniz.h"

namespace gar {
  namespace rec {



    //--------------------------------------------------------------------------
    TrackIoniz::TrackIoniz() {}



    //--------------------------------------------------------------------------
    void TrackIoniz::setData(std::vector<std::pair<float,float>> dSigdX_FWD,
                 std::vector<std::pair<float,float>> dSigdX_BAK) {
        fFWD_dSigdXs = dSigdX_FWD;
        fBAK_dSigdXs = dSigdX_BAK;
    }



  } // rec
} // gar
