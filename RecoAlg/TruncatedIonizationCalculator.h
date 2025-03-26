////////////////////////////////////////////////////////////////////////
// Class:       TruncatedIonizationCalculator
// File:        TruncatedIonizationCalculator.h
//
// Generated at Wed Feb 28 22:55 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#ifndef GAR_RECOALG_TruncatedIonizationCalculator_h
#define GAR_RECOALG_TruncatedIonizationCalculator_h

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"
#include <TMath.h>

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

struct CalibratedCaloScore {
  float dEdx_max_f1;
  float calibration_a;
  float calibration_b;
};

namespace fhicl {
  class ParameterSet;
}

namespace gar {
  namespace rec {
    namespace alg {

      class TruncatedIonizationCalculator {
      public:
        TruncatedIonizationCalculator(fhicl::ParameterSet const& pset);

        virtual ~TruncatedIonizationCalculator();

        void Configure(fhicl::ParameterSet const& pset);

        void LoadScorePars();

        void ClearLists();

        void PrepareAlgo(const rec::Track* track, const rec::TrackIoniz* ionization);

        void ComputeMeanIonization();

        std::pair<float, float> GetIonization();

        float GetdEdxProtonScore();

      private:
        std::vector<std::pair<float, float>> RegroupTrackClusters(
          std::vector<std::pair<float, float>> IonizationData,
          size_t nGroup);

        float CalculateTruncatedMean(std::vector<std::pair<float, float>> IonizationData,
                                     float percentage);

        float CalibrationFunction(float dQdx);

        std::vector<float> CalibrateIonization(std::vector<float> dQdxData);

        float TotalCaloEnergy(std::vector<std::pair<float, float>> IonizationData);

        // Configuration parameters
        int fNGroupCluster;
        float fTruncatePercent;
        float fIonizationEnergy;
        float fGroupGain;
        float fFitA;
        float fFitB;
        float fFitC;
        float fdQdxMax;
        float fdEdxMax;
        std::string fdEdxScoreParsFileName;

        std::map<std::pair<std::string, std::string>, CalibratedCaloScore> fScorerMap;

        float fTrackMomentum;

        std::vector<std::pair<float, float>> fSigDataFWD;
        std::vector<std::pair<float, float>> fSigDataBAK;

        float fTotalCalo;
        float fTruncatedMean;

        float Sigmoid(float x, float a, float b);
      };

    } // namespace alg
  }   // namespace rec
} // namespace gar

#endif /* GAR_RECOALG_TruncatedIonizationCalculator_h */
