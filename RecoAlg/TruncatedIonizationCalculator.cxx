////////////////////////////////////////////////////////////////////////
// Class:       TruncatedIonizationCalculator
// File:        TruncatedIonizationCalculator.cxx
//
// Generated at Wed Feb 28 22:55 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/TruncatedIonizationCalculator.h"
#include "RecoAlg/Loader.h"

namespace gar {
  namespace rec {
    namespace alg {

      //----------------------------------------------------------------------------
      TruncatedIonizationCalculator::TruncatedIonizationCalculator(fhicl::ParameterSet const& pset)
      {

        this->Configure(pset);

        this->LoadScorePars();

        return;
      }

      //----------------------------------------------------------------------------
      TruncatedIonizationCalculator::~TruncatedIonizationCalculator()
      {
        return;
      }

      //----------------------------------------------------------------------------
      void TruncatedIonizationCalculator::Configure(fhicl::ParameterSet const& pset)
      {

        fNGroupCluster          = pset.get<int>("NGroupCluster",              4);
        fTruncatePercent        = pset.get<float>("TruncatePercent",        0.6);
        fIonizationEnergy       = pset.get<float>("IonizationEnergy",   26.4e-6);
        fGroupGain              = pset.get<float>("GroupGain",             2.24);
        fFitA                   = pset.get<float>("FitA",                 0.883);
        fFitB                   = pset.get<float>("FitB",                  5.60);
        fFitC                   = pset.get<float>("FitC",                  4.94);
        fdQdxMax                = pset.get<float>("dQdxMax",              1.5e5);
        
        fdEdxMax = CalibrationFunction(fdQdxMax);
        
        fdEdxScoreParsFileName  = pset.get<std::string>("dEdxScoreParsFileName", "/pnfs/dune/persistent/users/fmlopez/GAr/MVAData/gar_dEdx_proton_score_v00_01_00.root");

        return;
      }

      //----------------------------------------------------------------------------
      void TruncatedIonizationCalculator::LoadScorePars()
      {
        // Load parameters from TTree used to assign "proton-ness" score
        // to particles based on <dE/dx> and momentum

        WildcardSource loader = WildcardSource(fdEdxScoreParsFileName);
        TFile *infile = loader.GetNextFile();
        TTree *tree = (TTree*) infile->Get("tree");

        std::vector<std::string>* _p_min     = 0;
        std::vector<std::string>* _p_max     = 0;
        Double_t _dEdx_max_f1;
        Double_t _calibration_a;
        Double_t _calibration_b;

        tree->SetBranchAddress("p_min",           &_p_min);
        tree->SetBranchAddress("p_max",           &_p_max);
        tree->SetBranchAddress("dEdx_max_f1",     &_dEdx_max_f1);
        tree->SetBranchAddress("calibrated_a",    &_calibration_a);
        tree->SetBranchAddress("calibrated_b",    &_calibration_b);

        // Read tree entries and create the map between (p0, sigmap) and calibration structs
        for(int i=0; i<tree->GetEntries(); i++){
          
          tree->GetEntry(i);
          CalibratedCaloScore calibration;

          calibration.dEdx_max_f1   = (float)_dEdx_max_f1;
          calibration.calibration_a = (float)_calibration_a;
          calibration.calibration_b = (float)_calibration_b;

          fScorerMap[std::make_pair(_p_min->at(0), _p_max->at(0))] = calibration;
        }

      }

      //----------------------------------------------------------------------------
      void TruncatedIonizationCalculator::ClearLists()
      {
        fTrackMomentum = 0.0;

        fSigDataFWD.clear();
        fSigDataBAK.clear();
      }

      //----------------------------------------------------------------------------
      void TruncatedIonizationCalculator::PrepareAlgo(const rec::Track* track, const rec::TrackIoniz* ionization)
      {

        //Clear the lists
        ClearLists();

        fTrackMomentum = 0.5*(track->Momentum_beg()+track->Momentum_end());

        // Fill corresponding ionization information
        std::vector<std::pair<float,float>> SigDataFWD;
        fSigDataFWD = ionization->getFWD_dSigdXs();
        fSigDataBAK = ionization->getBAK_dSigdXs();

        return;

      }

      //----------------------------------------------------------------------------
      void TruncatedIonizationCalculator::ComputeMeanIonization()
      {

        // Compute the sum of the deposited energy
        float TotalCaloFWD = TotalCaloEnergy(fSigDataFWD);
        float TotalCaloBAK = TotalCaloEnergy(fSigDataBAK);

        fTotalCalo = 0.5*(TotalCaloFWD+TotalCaloBAK);

        // Get the new ionization data regrouping the clusters in groups of fNGroupCluster clusters
        std::vector<std::pair<float,float>> NewSigDataFWD;
        NewSigDataFWD = RegroupTrackClusters(fSigDataFWD, fNGroupCluster);
        float TruncatedMeanFWD = CalculateTruncatedMean(NewSigDataFWD, fTruncatePercent);
        
        std::vector<std::pair<float,float>> NewSigDataBAK;
        NewSigDataBAK = RegroupTrackClusters(fSigDataBAK, fNGroupCluster);
        float TruncatedMeanBAK = CalculateTruncatedMean(NewSigDataBAK, fTruncatePercent);

        fTruncatedMean = 0.5*(TruncatedMeanFWD+TruncatedMeanBAK);

        return;
      }

      //----------------------------------------------------------------------------
      std::pair<float,float> TruncatedIonizationCalculator::GetIonization()
      {
        return std::make_pair(fTotalCalo, fTruncatedMean);
      }

      //----------------------------------------------------------------------------
      std::vector<std::pair<float, float>> TruncatedIonizationCalculator::RegroupTrackClusters(std::vector<std::pair<float, float>> IonizationData, size_t nGroup) {

        // Get the total number of input clusters...
        size_t nClusters = IonizationData.size();
        // ...the number of new clusters...
        size_t nNewClusters = nClusters/nGroup;
        // ...and the number of leftover clusters
        size_t nLeft = nClusters-nNewClusters*nGroup;

        // Create output collection
        std::vector<std::pair<float, float>> NewIonizationData;

        // Simply group the clusters in groups of nGroup
        // adding the energies and step sizes
        for(size_t i=0; i<nNewClusters; ++i){
          float newdE = 0;
          float newdX = 0;
          for(size_t j=0; j<nGroup; ++j){
            newdE += IonizationData[i*nGroup+j].first;
            newdX += IonizationData[i*nGroup+j].second;
          }
          std::pair<float, float> newData = std::make_pair(newdE, newdX);
          NewIonizationData.push_back(newData);
        }

        // Do not forget about the clusters that may be left
        // if the total number is not a multiple of nGroup
        if(nLeft > 0){
          float newdE = 0;
          float newdX = 0;
          for(size_t k=0; k<nLeft; ++k){
            newdE += IonizationData[nNewClusters*nGroup+k].first;
            newdX += IonizationData[nNewClusters*nGroup+k].second;
          }
          std::pair<float, float> newData = std::make_pair(newdE, newdX);
          NewIonizationData.push_back(newData); 
        }

        return NewIonizationData;

      }

      //----------------------------------------------------------------------------
      float TruncatedIonizationCalculator::CalculateTruncatedMean(std::vector<std::pair<float, float>> IonizationData, float percentage) {

        // Create dQdx vector from vector of pairs
        std::vector<float> dQdXvector;
        std::transform(IonizationData.begin(), IonizationData.end(), std::back_inserter(dQdXvector), [](const std::pair<float, float>& pair) -> float { return pair.first/pair.second; });

        // Apply energy calibration
        std::vector<float> dEdXvector;
        dEdXvector = CalibrateIonization(dQdXvector);

        std::sort(dEdXvector.begin(), dEdXvector.end());

        // Truncate, e.g. resize vector
        size_t newSize = dEdXvector.size()*percentage;
        dEdXvector.resize(newSize);

        float TruncatedMean;
        TruncatedMean = std::accumulate(dEdXvector.begin(), dEdXvector.end(), 0.0) / dEdXvector.size();

        return TruncatedMean;

      }

      //----------------------------------------------------------------------------
      float TruncatedIonizationCalculator::CalibrationFunction(float dQdx){
        return (std::exp(dQdx*fFitB*fIonizationEnergy/(fFitC*fGroupGain))-fFitA)/(fFitB)*1e3;
      }

      std::vector<float> TruncatedIonizationCalculator::CalibrateIonization(std::vector<float> dQdxData) {

        std::vector<float> dEdxData;
        std::transform(dQdxData.begin(), dQdxData.end(), std::back_inserter(dEdxData), [this](const float& dQdx) -> float { return CalibrationFunction(dQdx); });
        return dEdxData;

      }

      //----------------------------------------------------------------------------
      float TruncatedIonizationCalculator::TotalCaloEnergy(std::vector<std::pair<float, float>> IonizationData) {

        float TotalEnergy = 0.0;
        for (size_t i = 0; i<IonizationData.size(); ++i) {
          float deltadEdx = CalibrationFunction(IonizationData[i].first/IonizationData[i].second);
          
          if (deltadEdx <= fdEdxMax) {
            TotalEnergy += deltadEdx*IonizationData[i].second;
          } else {  
            TotalEnergy += fdEdxMax*IonizationData[i].second;
          }
        }

        return TotalEnergy;

      }

      //----------------------------------------------------------------------------
      float TruncatedIonizationCalculator::GetdEdxProtonScore()
      {

      float ProtonScore = 0.0;

      for (auto& [key, clf]: fScorerMap) {

          float p_min = std::stof(key.first);
          float p_max = std::stof(key.second);

          if ((fTrackMomentum >= p_min)&&(fTrackMomentum < p_max)) {
            // Apply corresponding probability calibration
            CalibratedCaloScore calibration = fScorerMap[std::make_pair(key.first, key.second)];
            ProtonScore = Sigmoid(fTruncatedMean-calibration.dEdx_max_f1, calibration.calibration_a, calibration.calibration_b);
            break;
          } else if ((fTrackMomentum >= 0.10)&&(fTrackMomentum < 0.30)) {
            // for momenta in the range 100-300 MeV we apply a simple cut at 50 keV/cm
            if (fTruncatedMean >= 50.0) ProtonScore = 1.0;
            break;
          } else if (fTrackMomentum < 0.10) {
            // for momenta lower than 100 MeV we apply a simple cut at 150 keV/cm
            if (fTruncatedMean >= 150.0) ProtonScore = 1.0;
            break;
          }

        }

        return ProtonScore;

      }

      //----------------------------------------------------------------------------
      float TruncatedIonizationCalculator::Sigmoid(float x, float a, float b)
      {
        return 1.0/(1.0+TMath::Exp(a*x+b));
      }

    } // namespace alg
  } // namespace rec
} // namespace gar
