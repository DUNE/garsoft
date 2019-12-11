//
//  Track.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  original Copyright 2016 Brian Rebel. All rights reserved.
//  Additions by Tom Junk, 2018
//  Little cute additions by Leo Bellantoni, 2019

#ifndef Track_hpp
#define Track_hpp

#include <stdio.h>
#include "RtypesCore.h"
#include <stdint.h>
#include "IDNumberGen.h"



namespace gar {
    namespace rec {



        // Sometimes we need to specify which end of the track we are using, as in 
        // vertexing.  If the Start/Vtx/beg end of track i is used, want true for 
        // the TrackEnd; if the far/End end, use false.  However, we can't use a bool
        // because the C++ STL does wierd stuff to make std::vector<bool> efficient...
        // thereby breaking all our code.
        typedef int TrackEnd;
        TrackEnd const TrackEndBeg = 1;      TrackEnd const TrackEndEnd = 0;



        class Track {

        public:
            Track();

            // let the compiler provide the dtor



        private:
            static gar::rec::IDNumber const FirstNumber = 100000;
            gar::rec::IDNumber fIDnumero;

            float fLengthforwards;    ///< length of the track in cm from forwards fit
            float fLengthbackwards;   ///< length of the track in cm from backwards fit
            float fMomentum_beg;  ///< momentum of the track at the vertex in GeV/c
            float fMomentum_end;  ///< momentum of the track at the end in GeV/c
            float fVertex[3]; ///< track vertex position in cm -- == "beginning" of track.  Arbitrary choice made by patrec
            float fEnd[3];    ///< track end    position in cm
            float fVtxDir[3]; ///< track vertex direction
            float fEndDir[3]; ///< track end    direction
            float fChisqForward; ///< chisquared forward fit
            float fChisqBackward; ///< chisquared backward fit
            size_t fNHits;        ///< number of hits
            double fTime;         ///< time in ns from trigger

            // use the x from fVertex and fEnd to specify the independent variable -- no need to store them twice

            float fTrackParBeg[5]; ///< Track parameters at beginning of track y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)
            float fTrackParEnd[5]; ///< Track parameters at end of track y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)

            float fCovMatBeg[15]; ///< covariance matrix at beginning of track -- packed in a 1D array, assuming symmetry
            float fCovMatEnd[15]; ///< covariance matrix at end of track

            #ifndef __GCCXML__

        public:

            // Cartesian coordinate constructor -- to be used with the MC Cheater, where
            // the uncertainties are not specified and are set to zero in the constructor.

            Track(const float  length,
            const float  momentum_beg,
            const float  momentum_end,
            const float *vtx,
            const float *end,
            const float *vtxDir,
            const float *endDir,
            const size_t nhits,
            const int    charge,   // units of e.  Need this to convert momentum into  curvature
            const double time);

            // constructor to fill after the fits -- including covariance matrix

            Track(const float lengthforwards,
            const float lengthbackwards,
            const size_t nhits,
            const float xbeg,           // x location at beginning of track in cm
            const float *trackparbeg,   // y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)
            const float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
            const float chisqforward,   // chisquared of forwards fit
            const float xend,           // x location at end of track
            const float *trackparend,   // y, z, curvature, phi, lambda  -- 5-parameter track (cm, cm, cm-1, radians, radians)
            const float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
            const float chisqbackward,  // chisquared of backwards fit
            const double time);      // timestamp

            //Track(gar::rec::TrackPar &tp);    // constructor using the track parameter class

            bool operator==(const Track& rhs) const;
            bool operator!=(const Track& rhs) const;
            gar::rec::IDNumber getIDNumber() const;

            const float* Vertex()   const;
            const float* End()      const;
            const float* VtxDir()   const;
            const float* EndDir()   const;
            float Length();  // returns average length
            float const& LengthForward()   const;  // returns length from forward fit
            float const& LengthBackward()  const;  // returns length from backward fit
            float const& Momentum_beg() const;
            float const& Momentum_end() const;
            float const& ChisqForward() const;
            float const& ChisqBackward() const;
            size_t const& NHits() const;
            const float* TrackParBeg() const;
            const float* TrackParEnd() const;
            void CovMatBegSymmetric(float *cmb) const; // fill a 5x5 array of floats owned by the caller with the values of the covariance matrix at the beginning of the track
            void CovMatEndSymmetric(float *cme) const; // fill a 5x5 array of floats owned by the caller with the values of the covariance matrix at the beginning of the track
            const float* CovMatBegPacked() const;
            const float* CovMatEndPacked() const;
            int ChargeBeg() const;
            int ChargeEnd() const;
            //   the charge assumes the track started at the beginning or the end above, and is expected to change sign
            // depending on which way the track is hypothesized to go
            double  const&       Time()      const;

            // expose this so we can use it elsewhere, and it's used twice in a constructor
            #endif

        };

        inline const float* Track::Vertex()   const { return fVertex;   }
        inline const float* Track::End()      const { return fEnd;      }
        inline const float* Track::VtxDir()   const { return fVtxDir;   }
        inline const float* Track::EndDir()   const { return fEndDir;   }
        inline float Track::Length() { return 0.5*(fLengthforwards+fLengthbackwards);  }
        inline float const& Track::LengthForward()   const { return fLengthforwards;   }
        inline float const& Track::LengthBackward()   const { return fLengthbackwards;   }
        inline float const& Track::Momentum_beg() const { return fMomentum_beg; }
        inline float const& Track::Momentum_end() const { return fMomentum_end; }
        inline float const& Track::ChisqForward() const { return fChisqForward; }
        inline float const& Track::ChisqBackward() const { return fChisqBackward; }
        inline size_t const& Track::NHits() const { return fNHits; }
        inline const float* Track::TrackParBeg() const { return fTrackParBeg; }
        inline const float* Track::TrackParEnd() const { return fTrackParEnd; }
        inline const float* Track::CovMatBegPacked() const { return fCovMatBeg; }
        inline const float* Track::CovMatEndPacked() const { return fCovMatEnd; }
        inline double const& Track::Time() const { return fTime; }

        // non-class functions

        // finds a unit vector pointing along the track momentum at one end.  If 
        // you want to extrapolate beyond the end, flip the sign of the direction
        void FindDirectionFromTrackParameters(const float *tparms, 
             const float thisXend,const float farXend, float *dir);


    } // rec
} // gar

#endif /* Track_hpp */
