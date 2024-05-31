////////////////////////////////////////////////////////////////////////
// Class:       ECALMuonBDT
// File:        ECALMuonBDT.cxx
//
// Generated at Mon Mar 04 16:25 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ECALMuonBDT.h"
#include "RecoAlg/Loader.h"

namespace gar {
  namespace rec {
    namespace alg {

      //----------------------------------------------------------------------------
      ECALMuonBDT::ECALMuonBDT(fhicl::ParameterSet const& pset, const gar::geo::GeometryCore* geo)
      {

        this->Configure(pset);

        fTPCCent.push_back(geo->TPCXCent());
        fTPCCent.push_back(geo->TPCYCent());
        fTPCCent.push_back(geo->TPCZCent());

        this->LoadClassifiers();

        return;
      }

      //----------------------------------------------------------------------------
      ECALMuonBDT::~ECALMuonBDT()
      {
        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::Configure(fhicl::ParameterSet const& pset)
      {

        fVerbosity           = pset.get<int>("Verbosity",              1);
        fBDTWeightDirectory  = pset.get<std::string>("BDTWeightDirectory", "$GARSOFT_DIR/MVAData/weights");
        fBDTSummaryFileName  = pset.get<std::string>("BDTSummaryFileName", "/pnfs/dune/persistent/users/fmlopez/GAr/MVAData/gar_bdt_summary_v00_01_00.root");
        fMaxMomentumECALOnly = pset.get<float>("MaxMomentumECALOnly", 0.8);
        fTMVAOutputMax       = pset.get<float>("TMVAOutputMax", 1.0);

        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::LoadClassifiers()
      {

        // Load parameters from TTree used to assign "muon-ness" score
        // to particles based on BDT output
        WildcardSource loader = WildcardSource(fBDTSummaryFileName);
        TFile *infile = loader.GetNextFile();
        TTree *tree = (TTree*) infile->Get("tree");

        std::vector<std::string>* _p0     = 0;
        std::vector<std::string>* _sigmap = 0;
        Double_t _n_estimators;
        Double_t _learning_rate;
        Double_t _calibration_a;
        Double_t _calibration_b;

        tree->SetBranchAddress("p0_value",           &_p0);
        tree->SetBranchAddress("sigmap_value",       &_sigmap);
        tree->SetBranchAddress("n_estimators",       &_n_estimators);
        tree->SetBranchAddress("learning_rate",      &_learning_rate);
        tree->SetBranchAddress("calibrated_a",       &_calibration_a);
        tree->SetBranchAddress("calibrated_b",       &_calibration_b);

        // Read tree entries and create the map between (p0, sigmap) and calibration structs
        for(int i=0; i<tree->GetEntries(); i++){
          tree->GetEntry(i);

          CalibrationBDT calibration;

          calibration.n_estimators  = (float)_n_estimators;
          calibration.learning_rate = (float)_learning_rate;
          calibration.calibration_a = (float)_calibration_a;
          calibration.calibration_b = (float)_calibration_b;

          fCalibrationMap[std::make_pair(_p0->at(0), _sigmap->at(0))] = calibration;
        }

        // Load BDT from the TMVA xml weight files
        
        std::string delimiter = "_"; // filenames are separated by underscores

        // Check all files in the provided directory
        for (const auto & entry : std::filesystem::directory_iterator(Wildcard(fBDTWeightDirectory).at(0))){

            // Get path and filename (without extension) of the current file
            std::string path     = entry.path().string();
            std::string filename = entry.path().stem();

            std::string p0     = "";
            std::string sigmap = "";

            size_t pos = 0;
            std::string token = ""; // initialise token to empty string

            // Find all instances of the delimiter in the filename
            while ((pos = filename.find(delimiter)) != std::string::npos) {
                if (token == "p0") {
                  // If the previous substring was "p0" the next one is the central momentum value
                  p0 = filename.substr(0, pos);
                } else if (token == "sigmap") {
                  // If the previous substring was "sigmap" the next one is the momentum width
                  sigmap = filename.substr(0, pos);
                }

                // Get substring until next delimiter
                token = filename.substr(0, pos);
                
                // Remove substring from string + delimiter length
                filename.erase(0, pos + delimiter.length());
            }

            // Create Reader object
            TMVA::Reader* reader = new TMVA::Reader("Silent");

            // Add variables to Reader, there must be a better way...
            if (std::stof(p0) >= fMaxMomentumECALOnly) {
              reader->AddVariable("ClusterTotalEnergy",                     &_ClusterTotalEnergy);
              reader->AddVariable("DistHitClusterMean",                     &_DistHitClusterMean);
              reader->AddVariable("DistHitClusterRMS",                      &_DistHitClusterRMS);
              reader->AddVariable("DistHitCenterMax",                       &_DistHitCenterMax);
              reader->AddVariable("TOFVelocity",                            &_TOFVelocity);
              reader->AddVariable("NLayers",                                &_NLayers);
              reader->AddVariable("NHits",                                  &_NHits);
              reader->AddVariable("HitMeanEnergy",                          &_HitMeanEnergy);
              reader->AddVariable("HitStdEnergy",                           &_HitStdEnergy);
              reader->AddVariable("HitMaxEnergy",                           &_HitMaxEnergy);
              reader->AddVariable("Radius90E",                              &_Radius90E);
              reader->AddVariable("ClusterMuIDTotalEnergy",                 &_ClusterMuIDTotalEnergy);
              reader->AddVariable("DistHitMuIDMax",                         &_DistHitMuIDMax);
              reader->AddVariable("DistHitCenterMuIDMax",                   &_DistHitCenterMuIDMax);
              reader->AddVariable("HitMuIDMeanEnergy",                      &_HitMuIDMeanEnergy);
              reader->AddVariable("HitMuIDStdEnergy",                       &_HitMuIDStdEnergy);
              reader->AddVariable("HitMuIDMaxEnergy",                       &_HitMuIDMaxEnergy);
              reader->AddVariable("NLayersMuID",                            &_NLayersMuID);
              reader->AddVariable("NHitsMuID",                              &_NHitsMuID);
              reader->AddVariable("TOFMuID",                                &_TOFMuID);
              reader->AddVariable("ClusterTotalEnergy/RecoMomentumFWD",     &_ClusterTotalEnergyOverRecoMomentumFWD);
              reader->AddVariable("ClusterMuIDTotalEnergy/RecoMomentumFWD", &_ClusterMuIDTotalEnergyOverRecoMomentumFWD);
            } else {
              reader->AddVariable("ClusterTotalEnergy",                     &_ClusterTotalEnergy);
              reader->AddVariable("HitMeanEnergy",                          &_HitMeanEnergy);
              reader->AddVariable("HitStdEnergy",                           &_HitStdEnergy);
              reader->AddVariable("HitMaxEnergy",                           &_HitMaxEnergy);
              reader->AddVariable("NHits",                                  &_NHits);
              reader->AddVariable("NLayers",                                &_NLayers);
              reader->AddVariable("DistHitClusterMean",                     &_DistHitClusterMean);
              reader->AddVariable("DistHitClusterRMS",                      &_DistHitClusterRMS);
              reader->AddVariable("DistHitCenterMax",                       &_DistHitCenterMax);
              reader->AddVariable("Radius90E",                              &_Radius90E);
              reader->AddVariable("TOFVelocity",                            &_TOFVelocity);
              reader->AddVariable("ClusterTotalEnergy/RecoMomentumFWD",     &_ClusterTotalEnergyOverRecoMomentumFWD);

            }

            // Book BDT classifier...
            reader->BookMVA("BDTG", path);
            // ...and add it to the classifer map
            fClassifierMap[std::make_pair(p0, sigmap)] = reader;
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::ResetVariables()
      {
        // Reset input ECAL data
        fNECALHits = 0;
        fECALTotalEnergy = 0.0;
        fECALHitEnergy.clear();
        fECALHitDistCluster.clear();
        fECALHitDistCentre.clear();
        fECALHitTime.clear();
        fECALHitLayer.clear();

        // Reset input MuID data
        fNMuIDHits = 0;
        fMuIDTotalEnergy = 0.0;
        fMuIDHitEnergy.clear();
        fMuIDHitPos.clear();
        fMuIDHitDistCentre.clear();
        fMuIDHitTime.clear();
        fMuIDHitLayer.clear();

        // Reset ECAL features
        fECALEnergyRatio        = 0.0;
        fECALHitDistClusterMean = -9999.;
        fECALHitDistClusterRMS  = -9999.;
        fECALHitDistCentreMin   = -9999.;
        fECALHitDistCentreMax   = -9999.;
        fECALHitEnergyMean      = -9999.;
        fECALHitEnergyVar       = -9999.;
        fECALHitEnergyStd       = -9999.;
        fECALHitEnergyMax       = -9999.;
        fECALHitLayerMin        = -1;
        fECALHitLayerMax        = -1;
        fECALHitNLayers         = 0;
        fECALToFVelocity        = -9999.;
        fECALRadius90E          = -9999.;
        fECALHitDistCentroid    = -9998.;
        fECALHitTimeCentroid    = -9998.;

        // Reset MuID features
        fMuIDEnergyRatio        = 0.0;
        fMuIDHitDistMax         = -9999.;
        fMuIDHitDistCentreMax   = -9999.;
        fMuIDHitEnergyMean      = -9999.;
        fMuIDHitEnergyVar       = -9999.;
        fMuIDHitEnergyStd       = -9999.;
        fMuIDHitEnergyMax       = -9999.;
        fMuIDHitLayerMin        = -1;
        fMuIDHitLayerMax        = -1;
        fMuIDHitNLayers         = 0;
        fMuIDHitDistCentroid    = -9999.;
        fMuIDHitTimeCentroid    = -9999.;
        fMuIDToFVelocity        = -9999.;

        // Reset value of muon score
        fMuonScore = 0.0;

      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::PrepareAlgo(const rec::Track* track)
      {

        ResetVariables();

        fTrackMomentum = 0.5*(track->Momentum_beg()+track->Momentum_end());

        return;

      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::AddECALHits(const rec::Cluster* ecal_cluster, std::vector<const rec::CaloHit*> ecal_hit_vec)
      {

        fECALTotalEnergy += ecal_cluster->Energy();

        float ClusterX = ecal_cluster->Position()[0]; float ClusterY = ecal_cluster->Position()[1]; float ClusterZ = ecal_cluster->Position()[2];
        std::vector<float> ClusterPos{ClusterX, ClusterY, ClusterZ};

        float ClusterMainAxisX = ecal_cluster->EigenVectors()[0]; float ClusterMainAxisY = ecal_cluster->EigenVectors()[1]; float ClusterMainAxisZ = ecal_cluster->EigenVectors()[2];
        std::vector<float> ClusterMainAxis{ClusterX+ClusterMainAxisX, ClusterY+ClusterMainAxisY, ClusterZ+ClusterMainAxisZ};

        for (auto const& ecal_hit: ecal_hit_vec) {
          fECALHitEnergy.push_back(ecal_hit->Energy());

          float HitX = ecal_hit->Position()[0]; float HitY = ecal_hit->Position()[1]; float HitZ = ecal_hit->Position()[2];
          std::vector<float> HitPos{HitX, HitY, HitZ};

          std::vector<float> HitDiff1;
          std::transform(HitPos.begin(), HitPos.end(), ClusterPos.begin(), std::back_inserter(HitDiff1), std::minus<float>());
          std::vector<float> HitDiff2;
          std::transform(HitPos.begin(), HitPos.end(), ClusterMainAxis.begin(), std::back_inserter(HitDiff2), std::minus<float>());

          std::vector<float> HitDistVec = CrossProduct1D(HitDiff1, HitDiff2);
          fECALHitDistCluster.push_back(std::sqrt(std::inner_product(HitDistVec.begin(), HitDistVec.end(), HitDistVec.begin(), 0.0)));

          std::vector<float> HitDistCent;
          std::transform(HitPos.begin(), HitPos.end(), fTPCCent.begin(), std::back_inserter(HitDistCent), std::minus<float>());
          fECALHitDistCentre.push_back(std::sqrt(std::inner_product(HitDistCent.begin(), HitDistCent.end(), HitDistCent.begin(), 0.0)));

          fECALHitTime.push_back(ecal_hit->Time().first);
          fECALHitLayer.push_back(ecal_hit->Layer());

          fNECALHits++;
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::AddMuIDHits(const rec::Cluster* muid_cluster, std::vector<const rec::CaloHit*> muid_hit_vec)
      {

        fMuIDTotalEnergy += muid_cluster->Energy();

        for (auto const& muid_hit: muid_hit_vec) {

          fMuIDHitEnergy.push_back(muid_hit->Energy());

          float HitMuIDX = muid_hit->Position()[0]; float HitMuIDY = muid_hit->Position()[1]; float HitMuIDZ = muid_hit->Position()[2];
          std::vector<float> HitMuIDPos{HitMuIDX, HitMuIDY, HitMuIDZ};

          fMuIDHitPos.push_back(HitMuIDPos);

          std::vector<float> HitDistMuIDCent;
          std::transform(HitMuIDPos.begin(), HitMuIDPos.end(), fTPCCent.begin(), std::back_inserter(HitDistMuIDCent), std::minus<float>());
          fMuIDHitDistCentre.push_back(std::sqrt(std::inner_product(HitDistMuIDCent.begin(), HitDistMuIDCent.end(), HitDistMuIDCent.begin(), 0.0)));

          fMuIDHitTime.push_back(muid_hit->Time().first);
          fMuIDHitLayer.push_back(muid_hit->Layer());
          
          fNMuIDHits++;

        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::ComputeFeatures()
      {

        if (fNECALHits > 0) {

          // Get ratio between energy deposited in ECAL and momentum measured in the TPC
          fECALEnergyRatio = fECALTotalEnergy/fTrackMomentum;

          // Find mean and RMS of the ditance between hits and cluster main axis
          fECALHitDistClusterMean = std::accumulate(fECALHitDistCluster.begin(), fECALHitDistCluster.end(), 0.0) / fECALHitDistCluster.size();
          fECALHitDistClusterRMS = std::sqrt(std::inner_product(fECALHitDistCluster.begin(), fECALHitDistCluster.end(), fECALHitDistCluster.begin(), 0.0) / fECALHitDistCluster.size());

          // Find min and max values of the distance between the hits and the center of the TPC
          fECALHitDistCentreMin = *std::min_element(fECALHitDistCentre.begin(), fECALHitDistCentre.end());
          fECALHitDistCentreMax = *std::max_element(fECALHitDistCentre.begin(), fECALHitDistCentre.end());

          // Get the first and last layer where hits were present
          fECALHitLayerMin = *std::min_element(fECALHitLayer.begin(), fECALHitLayer.end());
          fECALHitLayerMax = *std::max_element(fECALHitLayer.begin(), fECALHitLayer.end());
          fECALHitNLayers = fECALHitLayerMax-fECALHitLayerMin+1;

          // Obtain slope from the fit of hit time and hit distance to the center
          if(fNECALHits > 1){
            TGraph* graph = new TGraph(fECALHitTime.size(), &fECALHitTime[0], &fECALHitDistCentre[0]);
            TF1* linearFit = new TF1("linearFit", "pol1", fECALHitDistCentreMin, fECALHitDistCentreMax);
            graph->Fit(linearFit, "Q");
            fECALToFVelocity = linearFit->GetParameter(1);
          }

          // Find mean, std and max of the hit energy
          fECALHitEnergyMean = std::accumulate(fECALHitEnergy.begin(), fECALHitEnergy.end(), 0.0) / fECALHitEnergy.size();
          fECALHitEnergyVar  = std::accumulate(fECALHitEnergy.begin(), fECALHitEnergy.end(), 0.0, [&](float accumulator, const float& val) {return accumulator + ((val - fECALHitEnergyMean)*(val - fECALHitEnergyMean) / (fNECALHits - 1));});
          fECALHitEnergyStd = std::sqrt(fECALHitEnergyVar);
          fECALHitEnergyMax  = *std::max_element(fECALHitEnergy.begin(), fECALHitEnergy.end());

          // Create index for sorted values of the ditance between hits and cluster main axis
          std::vector<int> indices(fECALHitDistCluster.size());
          std::iota(indices.begin(), indices.end(), 0);
          std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool {return fECALHitDistCluster[A] < fECALHitDistCluster[B];});

          // Get radius containing 90% of the total ECAL energy
          float summed_energy = 0.;
          for (int j=0; j<fNECALHits; j++) {
            summed_energy += fECALHitEnergy[indices[j]];
            if (summed_energy >= 0.9*fECALTotalEnergy) {
              fECALRadius90E = fECALHitDistCluster[indices[j]];
              break;
            }
          }

          // Compute energy-weighted centroids for distance and time of ECAL hits
          fECALHitDistCentroid = std::inner_product(fECALHitDistCentre.begin(), fECALHitDistCentre.end(), fECALHitEnergy.begin(), 0.0)/std::accumulate(fECALHitEnergy.begin(), fECALHitEnergy.end(), 0.0);
          fECALHitTimeCentroid = std::inner_product(fECALHitTime.begin(), fECALHitTime.end(), fECALHitEnergy.begin(), 0.0)/std::accumulate(fECALHitEnergy.begin(), fECALHitEnergy.end(), 0.0);

        }

        if (fNMuIDHits > 0) {

          // Compute the maximum distance between a pair of hits in the MuID (measure of spread)
          for (size_t m=0; m<fMuIDHitPos.size(); m++) {
            for (size_t n=0; n<fMuIDHitPos.size(); n++){
              std::vector<float> HitDistMuID;
              std::transform(fMuIDHitPos[m].begin(), fMuIDHitPos[m].end(), fMuIDHitPos[n].begin(), std::back_inserter(HitDistMuID), std::minus<float>());
              float distHit_muid = std::sqrt(std::inner_product(HitDistMuID.begin(), HitDistMuID.end(), HitDistMuID.begin(), 0.0));
              if (distHit_muid > fMuIDHitDistMax) fMuIDHitDistMax = distHit_muid;
            }
          }

          // Get ratio between energy deposited in MuID and momentum measured in the TPC
          fMuIDEnergyRatio = fMuIDTotalEnergy/fTrackMomentum;

          // Find max values of the distance between the hits and the center of the TPC
          fMuIDHitDistCentreMax = *std::max_element(fMuIDHitDistCentre.begin(), fMuIDHitDistCentre.end());

          // Find mean, std and max of the hit energy
          fMuIDHitEnergyMean = std::accumulate(fMuIDHitEnergy.begin(), fMuIDHitEnergy.end(), 0.0) / fMuIDHitEnergy.size();
          fMuIDHitEnergyVar  = std::accumulate(fMuIDHitEnergy.begin(), fMuIDHitEnergy.end(), 0.0, [&](float accumulator, const float& val) {return accumulator + ((val - fMuIDHitEnergyMean)*(val - fMuIDHitEnergyMean) / (fNMuIDHits - 1));});
          fMuIDHitEnergyStd = std::sqrt(fMuIDHitEnergyVar);
          fMuIDHitEnergyMax  = *std::max_element(fMuIDHitEnergy.begin(), fMuIDHitEnergy.end());

          // Get the first and last layer where hits were present
          fMuIDHitLayerMin = *std::min_element(fMuIDHitLayer.begin(), fMuIDHitLayer.end());
          fMuIDHitLayerMax = *std::max_element(fMuIDHitLayer.begin(), fMuIDHitLayer.end());
          fMuIDHitNLayers = fMuIDHitLayerMax-fMuIDHitLayerMin+1;

          // Compute energy-weighted centroids for distance and time of ECAL hits
          fMuIDHitDistCentroid = std::inner_product(fMuIDHitDistCentre.begin(), fMuIDHitDistCentre.end(), fMuIDHitEnergy.begin(), 0.0)/std::accumulate(fMuIDHitEnergy.begin(), fMuIDHitEnergy.end(), 0.0);
          fMuIDHitTimeCentroid = std::inner_product(fMuIDHitLayer.begin(), fMuIDHitLayer.end(), fMuIDHitEnergy.begin(), 0.0)/std::accumulate(fMuIDHitEnergy.begin(), fMuIDHitEnergy.end(), 0.0);

          // Use centroids as a measure of velocity
          fMuIDToFVelocity = (fMuIDHitDistCentroid-fECALHitDistCentroid)/(fMuIDHitTimeCentroid-fECALHitTimeCentroid);
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALMuonBDT::ApplyClassifier()
      {

        if (fNECALHits == 0) return; // return if the track doesn't make it to the ECAL

        _ClusterTotalEnergy     = fECALTotalEnergy;
        _DistHitClusterMean     = fECALHitDistClusterMean;
        _DistHitClusterRMS      = fECALHitDistClusterRMS;
        _DistHitCenterMax       = fECALHitDistCentreMax;
        _TOFVelocity            = fECALToFVelocity;
        _NLayers                = (float)fECALHitNLayers;
        _NHits                  = (float)fNECALHits;
        _HitMeanEnergy          = fECALHitEnergyMean;
        _HitStdEnergy           = fECALHitEnergyStd;
        _HitMaxEnergy           = fECALHitEnergyMax;
        _Radius90E              = fECALRadius90E;
        _ClusterMuIDTotalEnergy = fMuIDTotalEnergy;
        _DistHitMuIDMax         = fMuIDHitDistMax;
        _DistHitCenterMuIDMax   = fMuIDHitDistCentreMax;
        _HitMuIDMeanEnergy      = fMuIDHitEnergyMean;
        _HitMuIDStdEnergy       = fMuIDHitEnergyStd;
        _HitMuIDMaxEnergy       = fMuIDHitEnergyMax;
        _NLayersMuID            = (float)fMuIDHitNLayers;
        _NHitsMuID              = (float)fNMuIDHits;
        _TOFMuID                = fMuIDToFVelocity;

        _ClusterTotalEnergyOverRecoMomentumFWD     = fECALEnergyRatio;
        _ClusterMuIDTotalEnergyOverRecoMomentumFWD = fMuIDEnergyRatio;

        for (auto& [key, clf]: fClassifierMap) {

          float p_min = std::stof(key.first)-std::stof(key.second);
          float p_max = std::stof(key.first)+std::stof(key.second);

          if ((fTrackMomentum >= p_min)&&(fTrackMomentum < p_max)) {
            fMuonScore = clf->EvaluateMVA("BDTG");

            if (abs(fMuonScore) >= fTMVAOutputMax) {
              // Whoops! Too big...
              fMuonScore = (fMuonScore ? (fMuonScore < 0) ? -1 : 1 : 0)*fTMVAOutputMax; // don't forget the sign
            }

            // Apply corresponding probability calibration
            CalibrationBDT calibration = fCalibrationMap[std::make_pair(key.first, key.second)];
            fMuonScore = ReEvaluateTMVA(fMuonScore, calibration.learning_rate, calibration.n_estimators, calibration.calibration_a, calibration.calibration_b);

            break;
          }
        }
      }

      //----------------------------------------------------------------------------
      float ECALMuonBDT::InverseTransformationTMVA(float x)
      {
        return TMath::Log((1+x)/(1-x))/2.0;
      }

      //----------------------------------------------------------------------------
      float ECALMuonBDT::Sigmoid(float x, float a, float b)
      {
        return 1.0/(1.0+TMath::Exp(a*x+b));
      }

      //----------------------------------------------------------------------------
      float ECALMuonBDT::ReEvaluateTMVA(float x, float learning_rate, float n_estimators, float a, float b)
      {
        float y = InverseTransformationTMVA(x);
        y = y*learning_rate*n_estimators;
        return Sigmoid(y, a, b);
      }

      //----------------------------------------------------------------------------
      std::pair<float, int> ECALMuonBDT::GetECALEnergy()
      {
        return std::make_pair(fECALTotalEnergy, fNECALHits);
      }

      //----------------------------------------------------------------------------
      std::pair<float, int> ECALMuonBDT::GetMuIDEnergy()
      {
        return std::make_pair(fMuIDTotalEnergy, fNMuIDHits);
      }

      //----------------------------------------------------------------------------
      float ECALMuonBDT::GetScore()
      {
        return fMuonScore;
      }

    } // namespace alg
  } // namespace rec
} // namespace gar
