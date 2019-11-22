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

      TrackPar(gar::rec::Track const &t, bool reversed=false);  // constructor from a Track data product

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

      const float *getTrackParametersBegin() const;
      const float *getTrackParametersEnd() const;
      const float *getCovMatBeg() const;
      const float *getCovMatEnd() const;

      const float getYCentBeg() const; // center of helix using beginning track parameters
      const float getZCentBeg() const;
      const float getYCentEnd() const; // center of helix using end track parameters
      const float getZCentEnd() const;
      const bool  getBegCentValid() const;  // is circle center calc valid (really just curvature != 0)
      const bool  getEndCentValid() const;

      size_t getNTPCClusters() const;
      float getLengthForwards() const;
      float getLengthBackwards() const;
      float getChisqForwards() const;
      float getChisqBackwards() const;
      float getXBeg() const;
      float getXEnd() const;
      //int getChargeBeg() const;  // just returns +1 or -1 depending on the sign of the curvature at the track beginning point
      //int getChargeEnd() const;  // just returns +1 or -1 depending on the sign of the curvature at the track ending point
      ULong64_t getTime() const;
      TVector3 getXYZBeg() const;
      TVector3 getXYZEnd() const;

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

      gar::rec::Track CreateTrack();    // Make a Track data product from this TrackPar instance
      void FitAnotherTrack(TrackPar &othertrack, float &chisquared, float *xyz, float *covmat); // find the best-fit vertex with another track

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
    };

  } // rec namespace
} // gar namespace

#endif
