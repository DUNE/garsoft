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

        Track::Track(const float  length,
        const float  momentum_beg,
        const float  momentum_end,
        const float *vtx,
        const float *end,
        const float *vtxDir,
        const float *endDir,
        const size_t nhits,
        const int   charge,
        const ULong64_t time)
        : fLengthforwards  (length  )
        , fLengthbackwards  (length  )
        , fMomentum_beg(momentum_beg)
        , fMomentum_end(momentum_end)
        , fChisqForward(0)
        , fChisqBackward(0)
        , fNHits(nhits)
        , fTime(time)
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
            float Pt_beg = momentum_beg * TMath::Sqrt(std::max(1E-6,1-TMath::Sq(vtxDir[0])));
            if (momentum_beg > 0)
            {
                fTrackParBeg[2] = (1.0E-2)*0.3*magfield[0]/Pt_beg;
            }
            else
            {
                fTrackParBeg[2] = 0;
            }
            fTrackParBeg[3] = TMath::ATan2(vtxDir[2],vtxDir[1]);

            if (vtxDir[1] != 0 || vtxDir[2] != 0)
            {
                fTrackParBeg[4] = TMath::ATan(vtxDir[0]/TMath::Sqrt(vtxDir[1]*vtxDir[1] + vtxDir[2]*vtxDir[2]));
            }
            else
            {
                fTrackParBeg[4] = -TMath::Pi()/2.0;
            }

            fTrackParEnd[0] = end[1];
            fTrackParEnd[1] = end[2];
            float Pt_end = momentum_end * TMath::Sqrt(std::max(1E-6,1-TMath::Sq(endDir[0])));
            if (momentum_end > 0)
            {
                fTrackParEnd[2] = (1.0E-2)*0.3*magfield[0]/Pt_end;
            }
            else
            {
                fTrackParEnd[2] = 0;
            }
            fTrackParEnd[3] = TMath::ATan2(endDir[2],endDir[1]);

            if (endDir[1] != 0 || endDir[2] != 0)
            {
                fTrackParEnd[4] = TMath::ATan(endDir[0]/TMath::Sqrt(endDir[1]*endDir[1] + endDir[2]*endDir[2]));
            }
            else
            {
                fTrackParEnd[4] = -TMath::Pi()/2.0;
            }

            return;
        }


        Track::Track(const float lengthforwards,
        const float lengthbackwards,
        const size_t nhits,
        const float xbeg,           // x location at beginning of track in cm
        const float *trackparbeg,   // y, z, curvature, phi, slope  -- 5-parameter track  (cm, cm, cm-1, radians, dy,z/dx)
        const float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
        const float chisqforward,   // chisquared of forwards fit
        const float xend,           // x location at end of track
        const float *trackparend,   // y, z, curvature, phi, slope  -- 5-parameter track (cm, cm, cm-1, radians, dy,z/dx)
        const float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
        const float chisqbackward,  // chisquared of forwards fit
        const ULong64_t time)       // timestamp
        : fLengthforwards(lengthforwards),
        fLengthbackwards(lengthbackwards)
<<<<<<< HEAD
      , fChisqForward(chisqforward)
      , fChisqBackward(chisqbackward)
      , fNHits(nhits)
      , fTime(time)
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

      FindDirectionFromTrackParameters(trackparbeg, fVtxDir);

      if (trackparbeg[2] != 0)
	{
	  float poverpt = 1.0/TMath::Sqrt(std::max(1E-6,1.0-TMath::Sq(fVtxDir[0])));
	  fMomentum_beg =  poverpt * TMath::Abs((1.0E-2)*0.3*magfield[0]/trackparbeg[2]);  // check constant (kG or T?)
	}
      else
	{
	  fMomentum_beg = 0;
	}

      fEnd[0] = xend;
      fEnd[1] = trackparend[0];
      fEnd[2] = trackparend[1];

      FindDirectionFromTrackParameters(trackparend, fEndDir);

      if (trackparend[2] != 0)
	{
	  float poverpt = 1.0/TMath::Sqrt(std::max(1E-6,1.0-TMath::Sq(fEndDir[0])));
	  fMomentum_end =  poverpt * TMath::Abs((1.0E-2)*0.3*magfield[0]/trackparend[2]);  // check constant (kG or T?)
	}
      else
	{
	  fMomentum_end = 0;
	}

    }



    // todo -- check consistency between charge at the beginning and end of track -- probably should
    // fix such a problem in the track fitter if the track appears to change charge sign along the way,
    // indicating that breaking the track in two is probably a good idea
    // n.b. our idea of what the charge is depends on which way the track is going
    // n.b. this has no concept of charge +-2 or other values than abs(charge)=1

    int Track::ChargeBeg() const
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

    int Track::ChargeEnd() const
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


    void FindDirectionFromTrackParameters(const float *tparms, float *dir)
    {
      dir[0] = -TMath::Tan(tparms[4]);
      dir[1] = -TMath::Sin(tparms[3]);
      dir[2] = -TMath::Cos(tparms[3]);
      float norm = TMath::Sqrt( 1.0 + dir[0]*dir[0]);
      dir[0] /= norm;
      dir[1] /= norm;
      dir[2] /= norm;
    }


  } // rec namespace
=======
        , fChisqForward(chisqforward)
        , fChisqBackward(chisqbackward)
        , fNHits(nhits)
        , fTime(time)
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

            FindDirectionFromTrackParameters(trackparbeg, fVtxDir);

            if (trackparbeg[2] != 0)
            {
                float poverpt = 1.0/TMath::Sqrt(std::max(1E-6,1.0-TMath::Sq(fVtxDir[0])));
                fMomentum_beg =  poverpt * TMath::Abs((1.0E-2)*0.3*magfield[0]/trackparbeg[2]);  // check constant (kG or T?)
            }
            else
            {
                fMomentum_beg = 0;
            }

            fEnd[0] = xend;
            fEnd[1] = trackparend[0];
            fEnd[2] = trackparend[1];

            FindDirectionFromTrackParameters(trackparend, fEndDir);

            if (trackparend[2] != 0)
            {
                float poverpt = 1.0/TMath::Sqrt(std::max(1E-6,1.0-TMath::Sq(fEndDir[0])));
                fMomentum_end =  poverpt * TMath::Abs((1.0E-2)*0.3*magfield[0]/trackparend[2]);  // check constant (kG or T?)
            }
            else
            {
                fMomentum_end = 0;
            }

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


        void FindDirectionFromTrackParameters(const float *tparms, float *dir)
        {
            dir[0] = -TMath::Tan(tparms[4]);
            dir[1] = -TMath::Sin(tparms[3]);
            dir[2] = -TMath::Cos(tparms[3]);
            float norm = TMath::Sqrt( 1.0 + dir[0]*dir[0]);
            dir[0] /= norm;
            dir[1] /= norm;
            dir[2] /= norm;
        }


    } // rec namespace
>>>>>>> 99080c065a325b43fd285adee0dc284c051764ce
} // gar namespace
