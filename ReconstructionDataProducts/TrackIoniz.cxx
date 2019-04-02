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
    // Two push, one pull methods

    void TrackIoniz::push_dSigdX(float dSigvalue, float dXvalue) {
      std::pair pushme = std::make_pair(dSigvalue,dXvalue);
      fdSigdXs.push_back( pushme );
    }

    void TrackIoniz::pull_dSigdX() {
      fdSigdXs.pop_back();
    }

    void TrackIoniz::push_dE_X  (float dEvalue,   float Xvalue) {
      std::pair pushme = std::make_pair(dEvalue,Xvalue);
      fdE_Xs.push_back( pushme );
    }


    //--------------------------------------------------------------------------
    // Two sort methods and their predicates

    bool lessThan_dE(std::pair<float,float> a, std::pair<float,float> b)
      {return a.first < b.first;}
    bool lessThan_X(std::pair<float,float> a, std::pair<float,float> b)
      {return a.second < b.second;}

    void TrackIoniz::sort_dE_X_by_dE() {
      std::sort(fdE_Xs.begin(),fdE_Xs.end(), lessThan_dE);
    }

    void TrackIoniz::sort_dE_X_by_X() {
      std::sort(fdE_Xs.begin(),fdE_Xs.end(), lessThan_X);
    }



  } // rec
} // gar
