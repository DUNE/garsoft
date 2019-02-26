//
//  TrackIoniz.cxx
//  garsoft-mrb
//
//  Created by Leo Bellantoni on 22 Feb, 2018

#include "ReconstructionDataProducts/TrackIoniz.h"

namespace gar {
  namespace rec {
     
    //--------------------------------------------------------------------------
    TrackIoniz::TrackIoniz()
      : fProcessed(false) {}



    //--------------------------------------------------------------------------
    // Two push methods

    void TrackIoniz::push_dSigdX(float dSigvalue, float dXvalue) {
      std::pair pushme = std::make_pair(dSigvalue,dXvalue);
      fdSigdXs.push_back( pushme );
    }

    void TrackIoniz::push_dE_X  (float dEvalue,   float Xvalue) {
      std::pair pushme = std::make_pair(dEvalue,Xvalue);
      fdE_Xs.push_back( pushme );
    }


  } // rec
} // gar
