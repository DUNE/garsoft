//
//  Track.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//  Modifications and additions by Tom Junk, 2018
//

#include "ReconstructionDataProducts/Track.h"
#include "nutools/MagneticField/MagneticField.h"
#include "TMath.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Track::Track()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    // Track constructor with no errors -- to be called by the Track Cheater

    Track::Track(float  length,
                 float  momentum_beg,
		 float  momentum_end,
                 float *vtx,
                 float *end,
                 float *vtxDir,
                 float *endDir,
		 size_t nhits,
		 int   charge)
      : fLength  (length  )
      , fMomentum_beg(momentum_beg)
      , fMomentum_end(momentum_end)
      , fChisqForward(0)
      , fChisqBackward(0)
      , fNHits(nhits)
    {

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      fVertex[0] = vtx[0];
      fVertex[1] = vtx[1];
      fVertex[2] = vtx[2];
      
      fVtxDir[0] = vtxDir[0];
      fVtxDir[1] = vtxDir[1];
      fVtxDir[2] = vtxDir[2];

      fEnd[0] = end[0];
      fEnd[1] = end[1];
      fEnd[2] = end[2];
      
      fEndDir[0] = endDir[0];
      fEndDir[1] = endDir[1];
      fEndDir[2] = endDir[2];

      for (size_t i=0; i<15; ++i)
	{
	  fCovMatBeg[i] = 0;
	  fCovMatEnd[i] = 0;
	}

      // calculate track parameters from information we have

      fTrackParBeg[0] = vtx[1];
      fTrackParBeg[1] = vtx[2];
      if (momentum_beg > 0)
	{
          fTrackParBeg[2] = 0.3*magfield[0]/momentum_beg;
	}
      else
	{
          fTrackParBeg[2] = 0;
	}
      fTrackParBeg[3] = TMath::ATan2(vtxDir[2],vtxDir[1]);

      if (vtxDir[0] != 0)
	{
	  fTrackParBeg[4] = TMath::Sqrt(vtxDir[0]*vtxDir[0] + vtxDir[1]*vtxDir[1])/vtxDir[0];
	}
      else
	{
	  fTrackParBeg[4] = -999.0;  // a better approximation than zero
	}

      fTrackParEnd[0] = end[1];
      fTrackParEnd[1] = end[2];
      if (momentum_beg > 0)
	{
          fTrackParEnd[2] = 0.3*magfield[0]/momentum_beg;
	}
      else
	{
          fTrackParEnd[2] = 0;
	}
      fTrackParEnd[3] = TMath::ATan2(endDir[2],endDir[1]);

      if (endDir[0] != 0)
	{
	  fTrackParEnd[4] = TMath::Sqrt(endDir[0]*endDir[0] + endDir[1]*endDir[1])/endDir[0];
	}
      else
	{
	  fTrackParEnd[4] = -999.0;  // a better approximation than zero
	}

      return;
    }


    Track::Track(float length,
		 size_t nhits,
		 float xbeg,           // x location at beginning of track in cm
		 float *trackparbeg,   // y, z, curvature, phi, slope  -- 5-parameter track  (cm, cm, cm-1, radians, dy,z/dx)
		 float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
		 float chisqforward,   // chisquared of forwards fit
		 float xend,           // x location at end of track
		 float *trackparend,   // y, z, curvature, phi, slope  -- 5-parameter track (cm, cm, cm-1, radians, dy,z/dx)
		 float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
		 float chisqbackward) // chisquared of forwards fit
      : fLength(length)
      , fChisqForward(chisqforward)
      , fChisqBackward(chisqbackward)
      , fNHits(nhits)
    {
      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      size_t icov = 0;
      for (size_t i=0; i< 5; ++i)
	{
	  fTrackParBeg[i] = trackparbeg[i];
	  fTrackParEnd[i] = trackparend[i];
	  for (size_t j=i; j<5; ++j)
	    {
	      fCovMatBeg[icov] = covmatbeg[i+5*j];
	      fCovMatEnd[icov] = covmatend[i+5*j];
	      ++icov;
	    }
	}

      fVertex[0] = xbeg;
      fVertex[1] = trackparbeg[0];
      fVertex[2] = trackparbeg[1];
      if (trackparbeg[2] != 0)
	{
          fMomentum_beg =  TMath::Abs(0.3*magfield[0]/trackparbeg[2]);  // check constant (kG or T?)
	}
      else
	{
	  fMomentum_beg = 0;
	}
      float hbeg = TMath::Sqrt( 1.0 + TMath::Sq(trackparbeg[4]) );
      float sbeg = TMath::Sin(trackparbeg[3]);
      float cbeg = TMath::Cos(trackparbeg[3]);
      fVtxDir[0] = trackparbeg[4]/hbeg;
      fVtxDir[1] = cbeg/hbeg;
      fVtxDir[2] = sbeg/hbeg;

      fEnd[0] = xend;
      fEnd[1] = trackparend[0];
      fEnd[2] = trackparend[1];
      if (trackparend[2] != 0)
	{
          fMomentum_end =  TMath::Abs(0.3*magfield[0]/trackparend[2]);  // check constant (kG or T?)
	}
      else
	{
	  fMomentum_end = 0;
	}
      float hend = TMath::Sqrt( 1.0 + TMath::Sq(trackparend[4]) );
      float send = TMath::Sin(trackparend[3]);
      float cend = TMath::Cos(trackparend[3]);
      fEndDir[0] = trackparend[4]/hend;
      fEndDir[1] = cend/hend;
      fEndDir[2] = send/hend;
    }


    
    // todo -- check consistency between charge at the beginning and end of track -- probably should
    // fix such a problem in the track fitter if the track appears to change charge sign along the way,
    // indicating that breaking the track in two is probably a good idea
    // n.b. our idea of what the charge is depends on which way the track is going
    // n.b. this has no concept of charge +-2 or other values than abs(charge)=1
 
    int Track::ChargeBeg()
    {
      int icharge=0;
      if (fTrackParBeg[2] > 0)
	{
	  icharge = 1;
	}
      else
	{
	  icharge = -1;
	}
      return icharge;
    }

    int Track::ChargeEnd()
    {
      int icharge=0;
      if (fTrackParEnd[2] > 0)
	{
	  icharge = 1;
	}
      else
	{
	  icharge = -1;
	}
      return icharge;
    }

    // recover symmetric covariance matrices.  Assume the caller owns the memory to a 5x5 array

    void Track::CovMatBegSymmetric(float *cmb) const
    {
      size_t ic=0;
      for (size_t i=0; i<5; ++i)
	{
	  for (size_t j=i; j<5; ++j)
	    {
	      cmb[i + 5*j] = fCovMatBeg[ic];
	      cmb[j + 5*i] = fCovMatBeg[ic];
	      ++ic;
	    }
	}
    }

    void Track::CovMatEndSymmetric(float *cmb) const
    {
      size_t ic=0;
      for (size_t i=0; i<5; ++i)
	{
	  for (size_t j=i; j<5; ++j)
	    {
	      cmb[i + 5*j] = fCovMatEnd[ic];
	      cmb[j + 5*i] = fCovMatEnd[ic];
	      ++ic;
	    }
	}
    }

  } // rec
} // gar

