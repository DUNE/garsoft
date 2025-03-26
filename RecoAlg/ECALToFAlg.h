////////////////////////////////////////////////////////////////////////
// Class:       ECALToFAlg
// File:        ECALToFAlg.h
//
// Generated at Tue Mar 05 16:15 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#ifndef GAR_RECOALG_ECALToFAlg_h
#define GAR_RECOALG_ECALToFAlg_h

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Track.h"

#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/BitFieldCoder.h"
#include "Geometry/GeometryGAr.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TVector3.h"

#include "TFile.h"
#include "TTree.h"

#include <Math/Integrator.h>
#include <Math/RootFinder.h>
#include <Math/WrappedFunction.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

struct CalibratedToFScore {
  float beta_max_f1;
  float calibration_a;
  float calibration_b;
};

namespace fhicl {
  class ParameterSet;
}

namespace gar {
  namespace rec {
    namespace alg {

      class ECALToFAlg {
      public:
        ECALToFAlg(fhicl::ParameterSet const& pset, const gar::geo::GeometryCore* geo);

        virtual ~ECALToFAlg();

        void Configure(fhicl::ParameterSet const& pset);

        void LoadScorePars();

        void ResetVariables();

        bool PrepareAlgo(const rec::Track* track, rec::TrackEnd track_end, float t0);

        bool ComputeEntryPoint();

        void AddHits(std::vector<const rec::CaloHit*> hit_vec);

        void ComputeArrivalTime();

        float GetTime();

        float GetBeta();

        float GetMass();

        float GetToFProtonScore();

      private:
        // Configuration parameters
        int fVerbosity; ///< level of verbosity for printouts
        float fMaxAngle;
        std::string fTimeMethod;
        std::string fToFScoreParsFileName;

        std::map<std::pair<std::string, std::string>, CalibratedToFScore> fScorerMap;

        std::string fECALEncoding;
        gar::geo::BitFieldCoder* fFieldDecoder_ECAL;

        float fDriftVelocity;

        std::vector<float> fTPCCent; ///< position of TPC from geometry service; 1 S Boston Ave.
        float fECALInnerRadius;
        float fECALStartX;

        std::unordered_map<int, std::tuple<float, float, float>> fHitMap;

        rec::TrackEnd fiTrackEnd;
        float fT0;

        float fTrackPar[5];
        float fTrackEnd[3];
        float fTrackLength;
        float fTrackMomentum;

        float fPhiMax = -1.;

        float fEntryPoint[3];
        float fEntryPhi;

        float fTrackLengthExtra;
        float fTrackLengthCorrected;

        bool fPropagateToEndCap;

        float fTime = -1.0;
        float fBeta = -1.0;
        float fMass = -1.0;

        void PropagateHitEndCap(const rec::CaloHit* hit, int Layer);

        void PropagateHitBarrel(const rec::CaloHit* hit, int Layer);

        float helix_circle_intersections(float phi,
                                         float y0,
                                         float z0,
                                         float R,
                                         float phi0,
                                         float r);

        float helix_x_intersections(float phi,
                                    float x0,
                                    float t0,
                                    float R,
                                    float lambda0,
                                    float phi0,
                                    float x);

        void helix(float phi, float* trackPar, float* trackEnd, float* projected, float t0);

        float beta(float length, float time);

        float mass(float momentum, float length, float time);

        float Sigmoid(float x, float a, float b);
      };

    } // namespace alg
  }   // namespace rec
} // namespace gar

#endif /* GAR_RECOALG_ECALToFAlg_h */
