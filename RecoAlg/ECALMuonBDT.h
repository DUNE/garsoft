////////////////////////////////////////////////////////////////////////
// Class:       ECALMuonBDT
// File:        ECALMuonBDT.h
//
// Generated at Mon Mar 04 16:25 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#ifndef GAR_RECOALG_ECALMuonBDT_h
#define GAR_RECOALG_ECALMuonBDT_h

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "TVector3.h"
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"

#include <cstdlib>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric>
#include <cmath>
#include <map>
#include <unordered_map>

#include <filesystem>

// 1D cross product
template <typename T>
std::vector<T> CrossProduct1D(std::vector<T> const &a, std::vector<T> const &b) {
  std::vector<T> r (a.size());  
  r[0] = a[1]*b[2]-a[2]*b[1];
  r[1] = a[2]*b[0]-a[0]*b[2];
  r[2] = a[0]*b[1]-a[1]*b[0];
  return r;
}

struct CalibrationBDT {
  float n_estimators;
  float learning_rate;
  float calibration_a;
  float calibration_b;
};

namespace fhicl{
  class ParameterSet;
}

namespace gar {
  namespace rec {
    namespace alg {

      class ECALMuonBDT {
      public:

        ECALMuonBDT(fhicl::ParameterSet const& pset, const gar::geo::GeometryCore* geo);

        virtual ~ECALMuonBDT();

        void Configure(fhicl::ParameterSet const& pset);

        void LoadClassifiers();

        void ResetVariables();

        void PrepareAlgo(const rec::Track* track);

        void AddECALHits(const rec::Cluster* ecal_cluster, std::vector<const rec::CaloHit*> ecal_hit_vec);

        void AddMuIDHits(const rec::Cluster* muid_cluster, std::vector<const rec::CaloHit*> muid_hit_vec);

        void ComputeFeatures();

        void ApplyClassifier();

        float InverseTransformationTMVA(float x);

        float Sigmoid(float x, float a, float b);

        float ReEvaluateTMVA(float x, float learning_rate, float n_estimators, float a, float b);

        std::pair<float, int> GetECALEnergy();

        std::pair<float, int> GetMuIDEnergy();

        float GetScore();

      private:

        // Configuration parameters
        int         fVerbosity;            ///< level of verbosity for printouts
        std::string fBDTWeightDirectory;   ///< directory containing BDT weight XML files
        std::string fBDTSummaryFileName;   ///< ROOT file with additional BDT information
        float       fMaxMomentumECALOnly;  ///< max momentum value that uses ECAL only, in GeV
        float       fTMVAOutputMax;        ///< truncate TMVA outputs beyound this cutoff (needed to apply calibration)

        std::vector<float> fTPCCent;       ///< position of TPC from geometry service; 1 S Boston Ave.

        float fTrackMomentum;

        // The map uses a pair of momentum values as keys, defining the region where we will apply the classifier [p_min, p_max)
        std::map<std::pair<std::string, std::string>, TMVA::Reader*> fClassifierMap;
        std::map<std::pair<std::string, std::string>, CalibrationBDT> fCalibrationMap;

        // Placeholder variables for classifier
        float _ClusterTotalEnergy;
        float _DistHitClusterMean;
        float _DistHitClusterRMS;
        float _DistHitCenterMax;
        float _TOFVelocity;
        float _NLayers;
        float _NHits;
        float _HitMeanEnergy;
        float _HitStdEnergy;
        float _HitMaxEnergy;
        float _Radius90E;
        float _ClusterMuIDTotalEnergy;
        float _DistHitMuIDMax;
        float _DistHitCenterMuIDMax;
        float _HitMuIDMeanEnergy;
        float _HitMuIDStdEnergy;
        float _HitMuIDMaxEnergy;
        float _NLayersMuID;
        float _NHitsMuID;
        float _TOFMuID;
        float _ClusterTotalEnergyOverRecoMomentumFWD;
        float _ClusterMuIDTotalEnergyOverRecoMomentumFWD;

        // ECAL hit properties
        float fECALTotalEnergy;
        int   fNECALHits;

        std::vector<float> fECALHitEnergy;
        std::vector<float> fECALHitDistCluster;
        std::vector<float> fECALHitDistCentre;
        std::vector<float> fECALHitTime;
        std::vector<int>   fECALHitLayer;

        // MuID hit properties
        int   fNMuIDHits;
        float fMuIDTotalEnergy;

        std::vector<float>              fMuIDHitEnergy;
        std::vector<std::vector<float>> fMuIDHitPos;
        std::vector<float>              fMuIDHitDistCentre;
        std::vector<float>              fMuIDHitTime;
        std::vector<int>                fMuIDHitLayer;

        // ECAL features
        float fECALEnergyRatio;

        float fECALHitDistClusterMean;
        float fECALHitDistClusterRMS;

        float fECALHitDistCentreMin;
        float fECALHitDistCentreMax;

        float fECALHitEnergyMean;
        float fECALHitEnergyVar;
        float fECALHitEnergyStd;
        float fECALHitEnergyMax;

        int fECALHitLayerMin;
        int fECALHitLayerMax;
        int fECALHitNLayers;

        float fECALToFVelocity;
        float fECALRadius90E;

        float fECALHitDistCentroid;
        float fECALHitTimeCentroid;

        // MuID features
        float fMuIDEnergyRatio;

        float fMuIDHitDistCentreMax;
        float fMuIDHitDistMax;

        float fMuIDHitEnergyMean;
        float fMuIDHitEnergyVar;
        float fMuIDHitEnergyStd;
        float fMuIDHitEnergyMax;

        int fMuIDHitLayerMin;
        int fMuIDHitLayerMax;
        int fMuIDHitNLayers;

        float fMuIDHitDistCentroid;
        float fMuIDHitTimeCentroid;

        float fMuIDToFVelocity;

        float fMuonScore;

      };

    } // namespace alg
  } // namespace rec
} // namespace gar

#endif /* GAR_RECOALG_ECALMuonBDT_h */
