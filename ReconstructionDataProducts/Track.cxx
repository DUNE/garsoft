//
//  Track.cpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright 2016 Brian Rebel. All rights reserved.  Unless he was working
//  for FRA in 2016 :) !
//  Modifications and additions by Tom Junk, 2018; Leo Bellantoni, 2019; etc.
//

#include "DetectorInfo/GArMagneticField.h"
#include "ReconstructionDataProducts/Track.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"
#include "CoreUtils/ServiceUtil.h"

#include "TMath.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        Track::Track()
        {
            // The default constructor is used e.g. by art::DataViewImpl::getHandle
            // Make sure all Track objects are numbered, lest art deep-copy uninitialized
            // Track instances and then operator==() evaluates meaninglessly true.
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();
            return;
        }



        bool Track::operator==(const Track& rhs) const {
            return (this->fIDnumero == rhs.fIDnumero);
        }

        bool Track::operator!=(const Track& rhs) const {
            return (this->fIDnumero != rhs.fIDnumero);
        }

        gar::rec::IDNumber Track::getIDNumber() const {return fIDnumero;}



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
	const int   , //charge,
        const double time)
        : fLengthforwards  (length  )
        , fLengthbackwards  (length  )
        , fMomentum_beg(momentum_beg)
        , fMomentum_end(momentum_end)
        , fChisqForward(0)
        , fChisqBackward(0)
        , fNHits(nhits)
        , fTime(time)
        {
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();

            auto const* magFieldService = gar::providerFrom<mag::MagneticFieldService>();
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
            fTrackParBeg[3] = TMath::ATan2(vtxDir[1],vtxDir[2]);

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
            fTrackParEnd[3] = TMath::ATan2(endDir[1],endDir[2]);

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
        const double time)       // timestamp
        : fLengthforwards(lengthforwards)
        , fLengthbackwards(lengthbackwards)
        , fChisqForward(chisqforward)
        , fChisqBackward(chisqbackward)
        , fNHits(nhits)
        , fTime(time)
        {
          IDNumberGen::create(FirstNumber);
          fIDnumero = IDNumberGen::create()->getNewOne();

          auto const* magFieldService = gar::providerFrom<mag::MagneticFieldService>();
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

          FindDirectionFromTrackParameters(trackparbeg, xbeg,xend, fVtxDir);

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

          FindDirectionFromTrackParameters(trackparend, xend,xbeg, fEndDir);

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



        // Charge is d(phi)/d(time in flight).  We don't know if it went from fVertex
        // to fEnd or the otherway round; but under either assumption we can find
        // if phi is increasing or decreasing.  All we need is the x component of
        // the cross product of the direction with the lever arm (position of track
        // end in the (z,y) plane relative to (Zc, Yc).
        // N.b. this has no concept of charge +-2 or other values than abs(charge)=1

        int Track::ChargeBeg() const {
            float Zlever = +TMath::Sin(fTrackParBeg[3]) / fTrackParBeg[2];
            float Ylever = -TMath::Cos(fTrackParBeg[3]) / fTrackParBeg[2];
            float Xcross = Ylever * fVtxDir[2] - Zlever * fVtxDir[1];
            return ( Xcross>0 ? -1 : +1 );
        }

        int Track::ChargeEnd() const {
            float Zlever = +TMath::Sin(fTrackParEnd[3]) / fTrackParEnd[2];
            float Ylever = -TMath::Cos(fTrackParEnd[3]) / fTrackParEnd[2];
            float Xcross = Ylever * fEndDir[2] - Zlever * fEndDir[1];
            return ( Xcross>0 ? -1 : +1 );
        }



        // recover symmetric covariance matrices.  Assume the caller owns the memory to a 5x5 array
        void Track::CovMatBegSymmetric(float *cmb) const {
            size_t ic=0;
            for (size_t i=0; i<5; ++i) {
                for (size_t j=i; j<5; ++j) {
                    cmb[i + 5*j] = fCovMatBeg[ic];
                    cmb[j + 5*i] = fCovMatBeg[ic];
                    ++ic;
                 }
            }
         }

        void Track::CovMatEndSymmetric(float *cmb) const {
            size_t ic=0;
            for (size_t i=0; i<5; ++i) {
                for (size_t j=i; j<5; ++j) {
                    cmb[i + 5*j] = fCovMatEnd[ic];
                    cmb[j + 5*i] = fCovMatEnd[ic];
                    ++ic;
                }
             }
         }


        void FindDirectionFromTrackParameters(const float *tparms,
            const float thisXend, const float farXend, float *dir) {
            // Direction is either +1 or -1 times d (position) / d(phi)
            // from the track fit equations; the sign is determined by the
            // x position of the far end.
            dir[0] = TMath::Tan(tparms[4]);
            dir[1] = TMath::Sin(tparms[3]);
            dir[2] = TMath::Cos(tparms[3]);

            int sigh = +1;
            if ( dir[0]>0 && farXend<thisXend ) sigh = -1;
            if ( dir[0]<0 && farXend>thisXend ) sigh = -1;
            float norm = sigh * TMath::Sqrt( 1.0 + dir[0]*dir[0]);

            dir[0] /= norm;
            dir[1] /= norm;
            dir[2] /= norm;
        }



    } // rec namespace
} // gar namespace
