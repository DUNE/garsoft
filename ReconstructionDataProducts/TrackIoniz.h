//
//  TrackIoniz.h
//
//  Created by Leo Bellantoni on 22 Feb 2019.
//  Track ionization data -- use associations with tracks to find them.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_TrackIoniz_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_TrackIoniz_h

#include "RtypesCore.h"
#include <iostream>
#include <vector>
#include <utility>

namespace gar {
  namespace rec {

    class TrackIoniz {

    public:
      TrackIoniz();

      // let the compiler provide the dtor

    private:

      std::vector<std::pair<float,float>> fdSigdXs;      ///< Segment values ordered by forward fit
      std::vector<std::pair<float,float>> fdE_Xs;        ///< Is dE/dX vs length along trackk
      bool fProcessed;                                   ///< Is calibrated & valid dEdX available?


#ifndef __GCCXML__

    public:

      void push_dSigdX(float dSigvalue, float dXvalue);
      void push_dE_X  (float dEvalue,   float  Xvalue);

      const std::vector<std::pair<float,float>> dSigdXvalues();
      const std::vector<std::pair<float,float>> dE_Xvalues();

      const bool Processed();
      void setProcessedFlag();  

#endif

    };

    inline const std::vector<std::pair<float,float>> TrackIoniz::dSigdXvalues() {return fdSigdXs;}
    inline const std::vector<std::pair<float,float>> TrackIoniz::dE_Xvalues()   {return fdE_Xs;}
    inline const bool TrackIoniz::Processed()        {return fProcessed;}
    inline void       TrackIoniz::setProcessedFlag() {fProcessed = true; return;}

//  This method is not linkable from gar::processIonizationInfo(gar::rec::TrackIoniz)
//  if it is in the ---.cxx file, although push_dSigdX can be found from 
//  gar::rec::tpctrackfit2::KalmanFit.   Why? Why? Why?
/* WTF
    void TrackIoniz::push_dE_X  (float dEvalue,   float Xvalue) {
      std::pair pushme = std::make_pair(dEvalue,Xvalue);
      fdE_Xs.push_back  ( pushme );
    }
FTW*/


  } // rec
} // gar


#endif
