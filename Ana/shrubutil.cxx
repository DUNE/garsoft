////////////////////////////////////////////////////////////////////////
// Class:       no class.
/// \file       shrubutil.cxx
/// \brief Code to compute quantities that will go into the anatree
///
/// \version $Id:  $
/// \author  bellanto@fnal.gov
//
// This is a separate file of functions used by anatree_module.cc
// It exists to keep that file of manageable size.
//
////////////////////////////////////////////////////////////////////////



#include "Ana/shrubutil.h"



namespace gar {



    //==========================================================================
    // Process ionization.  Eventually this moves into the reco code.
    void processIonizationInfo_shrub( rec::TrackIoniz& ion, float ionizeTruncate,
                                float& forwardIonVal, float& backwardIonVal ) {
        // Get the ADC data
        std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();

        // NO CALIBRATION SERVICE FOR NOW

        forwardIonVal = processOneDirection_shrub(SigData, ionizeTruncate);

        SigData = ion.getBAK_dSigdXs();
        backwardIonVal = processOneDirection_shrub(SigData, ionizeTruncate);

        return;
    }



    float processOneDirection_shrub(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

        std::vector<std::pair<float,float>> dEvsX;    // Will be the ionization vs distance along track

        // The first hit on the track never had its ionization info stored.  Not a problem
        // really.  Each pair is a hit and the step along the track that ends at the hit
        // For the last hit, just take the step from the n-1 hit; don't guess some distance to
        // (nonexistant!) n+1 hit.  Using pointer arithmetic because you are a real K&R C nerd!
        float distAlongTrack = 0;
        std::vector<std::pair<float,float>>::iterator littlebit = SigData.begin();
        for (; littlebit<(SigData.end()-1); ++littlebit) {
            float dE =   std::get<0>(*littlebit);
           // tpctrackfit2_module.cc fills the TrackIoniz data product so that 
           // this quantity is really dL > 0 not dX, a coordinate on the drift axis
            float dX  = std::get<1>(*littlebit);
            distAlongTrack += dX;    // But count full step to get hit position on track
            // Take dX to be 1/2 the previous + last segment
            dX += std::get<1>(*(littlebit+1));
            float dEdX = dE/(0.5*dX);

            std::pair pushme = std::make_pair(dEdX,distAlongTrack);
            dEvsX.push_back( pushme );
        }

        // Get the truncated mean; first sort then take mean
        std::sort(dEvsX.begin(),dEvsX.end(), lessThan_byE_shrub);

        // Get the dEdX vs length data, truncated.
        int goUpTo = ionizeTruncate * dEvsX.size() +0.5;
        if (goUpTo > (int)dEvsX.size()) goUpTo = dEvsX.size();
        int i = 1;        float returnvalue = 0;
        littlebit = dEvsX.begin();
        for (; littlebit<dEvsX.end(); ++littlebit) {
          returnvalue += std::get<0>(*littlebit);
          ++i;
          if (i>goUpTo) break;
        }
        returnvalue /= goUpTo;
        return returnvalue;
    }




    //==========================================================================
    // Coherent pion analysis specific code
    float computeT_shrub( simb::MCTruth theMCTruth ) {
        // Warning.  You probably want the absolute value of t, not t.
        int nPart = theMCTruth.NParticles();
        enum { nu, mu, pi};
        float E[3], Px[3], Py[3], Pz[3];
        E[nu] = E[mu] = E[pi] = -1e42;

        for (int i=0; i<3;++i) {
            Px[i] = 0; 
            Py[i] = 0;
            Pz[i] = 0;
            E[i]  = 0;
        }
        // Find t from the MCParticles via the
        for (int iPart=0; iPart<nPart; iPart++) {
            simb::MCParticle Part = theMCTruth.GetParticle(iPart);
            int code = Part.PdgCode();
            int mom  = Part.Mother();

            // get the neutrino
            if ( abs(code) == 12 || abs(code) == 14 || abs(code) == 16 ) {
                if (mom == -1) {
                    E[nu] = Part.E();   Px[nu] = Part.Px();   Py[nu] = Part.Py();   Pz[nu] = Part.Pz();
                }
            }

            // get the lepton
            if ( abs(code) == 11 || abs(code) == 13 || abs(code) == 15 ) {
                if (mom == 0) {
                    E[mu] = Part.E();   Px[mu] = Part.Px();   Py[mu] = Part.Py();   Pz[mu] = Part.Pz();
                }
            }

            // get the pion
            if ( code==111 || abs(code)==211 ) {
                if (mom == 1) {
                    E[pi] = Part.E();   Px[pi] = Part.Px();   Py[pi] = Part.Py();   Pz[pi] = Part.Pz();
                }
            }

            // get outa here
            if ( E[nu]!=0 && E[mu]!=0 && E[pi]!=0) break;

        }

        // Compute t; reuse nu 4-vector to get first q, then t.
        E[nu] -= E[mu];   Px[nu] -= Px[mu];   Py[nu] -= Py[mu];   Pz[nu] -= Pz[mu];
        E[nu] -= E[pi];   Px[nu] -= Px[pi];   Py[nu] -= Py[pi];   Pz[nu] -= Pz[pi];
        float t = E[nu]*E[nu] -Px[nu]*Px[nu] -Py[nu]*Py[nu] -Pz[nu]*Pz[nu];
        return t;
    }



} // end namespace gar
