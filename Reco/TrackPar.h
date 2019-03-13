//
// TrackPar.h
//
// Created by Tom Junk, July 2, 2018
// 
// Reasoning to keep this out of gar::rec::Track -- the art data product design guide discourages
// the mixture of i/o and algorithms.  This class contains algorithms needed to manipulate TPCClusters and tracks.
//

#ifndef GARTPCTRACKPAR_H
#define GARTPCTRACKPAR_H

#include "ReconstructionDataProducts/Track.h"
#include "TVector3.h"

namespace gar
{
  namespace rec
  {

    class Track;

    class TrackPar
    {
    public:

      TrackPar(gar::rec::Track const &t);  // constructor from a Track data product

      TrackPar(const float lengthforwards,  // constructor from parameters
	       const float lengthbackwards,
	       const size_t nTPCClusters,
	       const float xbeg,           // x location at beginning of track in cm
	       const float *trackparbeg,   // y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)
	       const float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
	       const float chisqforward,   // chisquared of forwards fit
	       const float xend,           // x location at end of track
	       const float *trackparend,   // y, z, curvature, phi, lambda  -- 5-parameter track (cm, cm, cm-1, radians, lambda)
	       const float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
	       const float chisqbackward,  // chisquared of backwards fit
	       const ULong64_t time);      // timestamp

      TrackPar() {};   // empty constructor

      const float *getTrackParametersBegin();
      const float *getTrackParametersEnd();
      const float *getCovMatBeg();
      const float *getCovMatEnd();
      size_t getNTPCClusters();   
      float getLengthForwards();
      float getLengthBackwards();
      float getChisqForwards();
      float getChisqBackwards();
      float getXBeg();
      float getXEnd();
      int getChargeBeg();   // just returns +1 or -1 depending on the sign of the curvature at the track beginning point
      int getChargeEnd();  // just returns +1 or -1 depending on the sign of the curvature at the track ending point
      ULong64_t getTime();
      TVector3 getXYZBeg();
      TVector3 getXYZEnd();
      TVector3 getPosAtX(const float x, bool usebegpar=true);  // returns 3d track position if x is provided, using track begin track parameters  
      // if usebegpar is true, otherwise use the end parameters

      void setNTPCClusters(const size_t nTPCClusters);
      void setTrackParametersBegin(const float *tparbeg);
      void setTrackParametersEnd(const float *tparend);
      void setCovMatBeg(const float *covmatbeg);
      void setCovMatEnd(const float *covmatend);

      void setLengthForwards(const float lengthforwards);
      void setLengthBackwards(const float lengthbackwards);
      void setChisqForwards(const float chisqforwards);
      void setChisqBackwards(const float chisqbackwards);
      void setXBeg(const float xbeg);
      void setXEnd(const float xend);
      void setTime(const ULong64_t time);

      float DistXYZ(const float *xyz);  // Distance from a point to this track (minimum using forwards and backwards track parameters).
      gar::rec::Track CreateTrack();    // Make a Track data product from this TrackPar instance
      void FitAnotherTrack(TrackPar &othertrack, float &chisquared, float *xyz, float *covmat); // find the best-fit vertex with another track

      // extrapolate both ends of a track to a radius from the center of the detector and provide 3x3 covariance matrices
      void ExtrapolateToRadius(float radius, bool &extrapsucceeded1, float *xyz1, float *covmat1, bool &extrapsucceeded2, float *xyz2, float *covmat2);

      // extrapolate one end of a track to a radius from the center of the detector and provide 3x3 covariance matrices
      // returns true if extrapolation succeeded, false otherwise
      bool ExtrapolateOneEndToRadius(float radius, float *trackpar, float xpar, float *xyz, float *covmat);

    private:
      size_t fNTPCClusters;
      float fXBeg;  // X at which the beginning track parameters are quoted
      float fTrackParametersBegin[5];  // y, z, curvature, phi, lambda
      float fXEnd;
      float fTrackParametersEnd[5]; // y, z, curvature, phi, lambda
      float fChisquaredForwards;
      float fChisquaredBackwards;
      float fLengthForwards;
      float fLengthBackwards;
      float fCovMatBeg[25];
      float fCovMatEnd[25];
      float fYCentBeg;      // center of helix using beginning track parameters
      float fZCentBeg;      
      float fYCentEnd;      // center of helix using end track parameters
      float fZCentEnd;
      bool fBegCentValid;   // flags to indicate if circle center calc is valid (really just curvature != 0)
      bool fEndCentValid;
      ULong64_t fTime;      // timestamp

      void CalcCenter();    // so as not to duplicate code in the constructors
      float DistXYZ_Aux(const float *xyz, bool useBegin);
    };

  } // rec namespace
} // gar namespace

#endif
