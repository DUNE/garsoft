////////////////////////////////////////////////////////////////////////
// Class:       ECALToFAlg
// File:        ECALToFAlg.cxx
//
// Generated at Tue Mar 05 16:15 BST 2024 by Francisco Martinez Lopez
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ECALToFAlg.h"
#include "RecoAlg/Loader.h"

namespace gar {
  namespace rec {
    namespace alg {

      //----------------------------------------------------------------------------
      ECALToFAlg::ECALToFAlg(fhicl::ParameterSet const& pset, const gar::geo::GeometryCore* geo)
      {

        this->Configure(pset);

        this->LoadScorePars();

        auto detProp   = gar::providerFrom<detinfo::DetectorPropertiesService>();
        fDriftVelocity = detProp->DriftVelocity(detProp->Efield(), detProp->Temperature());

        fECALEncoding = geo->GetECALCellIDEncoding();
        fFieldDecoder_ECAL = new gar::geo::BitFieldCoder(fECALEncoding);

        fTPCCent.push_back(geo->TPCXCent());
        fTPCCent.push_back(geo->TPCYCent());
        fTPCCent.push_back(geo->TPCZCent());

        fECALInnerRadius = geo->GetECALInnerBarrelRadius();
        fECALStartX      = geo->GetECALEndcapStartX();

        return;
      }

      //----------------------------------------------------------------------------
      ECALToFAlg::~ECALToFAlg()
      {
        return;
      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::Configure(fhicl::ParameterSet const& pset)
      {

        fVerbosity             = pset.get<int>("Verbosity",       1);
        fMaxAngle              = pset.get<float>("MaxAngle",      1.0);
        fTimeMethod            = pset.get<std::string>("TimeMethod",      "Average");
        fToFScoreParsFileName  = pset.get<std::string>("ToFScoreParsFileName", "/pnfs/dune/persistent/users/fmlopez/GAr/MVAData/gar_tof_proton_score_v00_01_00.root");

        if (fTimeMethod.compare("Earliest") == 0) {
          MF_LOG_DEBUG("ECALToFAlg") << "Using time of earliest hit as arrival time";
        } else if (fTimeMethod.compare("Average") == 0) {
          MF_LOG_DEBUG("ECALToFAlg") << "Using average to compute arrival times";
        } else if (fTimeMethod.compare("Fit") == 0) {
          MF_LOG_DEBUG("ECALToFAlg") << "Using intercept of fit as arrival time";
        } else {
          throw cet::exception("ECALToFAlg") << "Unable to determine which algorithm to use, bail";
        }

        return;
      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::LoadScorePars()
      {
        // Load parameters from TTree used to assign "proton-ness" score
        // to particles based on time-of-flight and momentum
        //WildcardSource loader = WildcardSource(fToFScoreParsFileName);
        std::vector<std::string> FileVector = { fToFScoreParsFileName };
        FileListSource loader = FileListSource(FileVector);

        TFile *infile = loader.GetNextFile();
        TTree *tree = (TTree*) infile->Get("tree");

        std::vector<std::string>* _p_min     = 0;
        std::vector<std::string>* _p_max     = 0;
        Double_t _beta_max_f1;
        Double_t _calibration_a;
        Double_t _calibration_b;

        tree->SetBranchAddress("p_min",           &_p_min);
        tree->SetBranchAddress("p_max",           &_p_max);
        tree->SetBranchAddress("beta_max_f1",     &_beta_max_f1);
        tree->SetBranchAddress("calibrated_a",    &_calibration_a);
        tree->SetBranchAddress("calibrated_b",    &_calibration_b);

        // Read tree entries and create the map between (p0, sigmap) and calibration structs
        for(int i=0; i<tree->GetEntries(); i++){
          tree->GetEntry(i);

          CalibratedToFScore calibration;

          calibration.beta_max_f1   = (float)_beta_max_f1;
          calibration.calibration_a = (float)_calibration_a;
          calibration.calibration_b = (float)_calibration_b;

          fScorerMap[std::make_pair(_p_min->at(0), _p_max->at(0))] = calibration;
        }

      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::ResetVariables()
      {
        fHitMap.clear();
        fTrackLengthExtra  = 0.0;
        fPropagateToEndCap = false;

        fTime = -1.0;
        fBeta = -1.0;
        fMass = -1.0;
      }

      //----------------------------------------------------------------------------
      bool ECALToFAlg::PrepareAlgo(const rec::Track* track, rec::TrackEnd track_end, float t0)
      {
        ResetVariables();

        fT0 = t0;
        fiTrackEnd = track_end;

        // Based on the track end associated, get the track position and fit parameters at that point
        if (fiTrackEnd == gar::rec::TrackEndBeg) {
          for (int i=0; i<5; ++i) fTrackPar[i] = track->TrackParBeg()[i];
          for (int i=0; i<3; ++i) fTrackEnd[i] = track->Vertex()[i];
          fTrackLength = track->LengthForward();
          fTrackMomentum = track->Momentum_beg();
        } else if (fiTrackEnd == gar::rec::TrackEndEnd) {
          for (int i=0; i<5; ++i) fTrackPar[i] = track->TrackParEnd()[i];
          for (int i=0; i<3; ++i) fTrackEnd[i] = track->End()[i];
          fTrackLength = track->LengthBackward();
          fTrackMomentum = track->Momentum_end();
        }

        if(fTrackPar[2] > 0) {
            fPhiMax = fTrackPar[3]+fMaxAngle*TMath::Pi();
        } else {
            fPhiMax = fTrackPar[3]-fMaxAngle*TMath::Pi();
        }

        return ComputeEntryPoint();

      }

      //----------------------------------------------------------------------------
      bool ECALToFAlg::ComputeEntryPoint()
      {

        auto helix_circle_intersections_to_wrap_entry = [this](float phi) {
            const float y0   = fTrackPar[0];
            const float z0   = fTrackPar[1];
            const float R    = 1/fTrackPar[2];
            const float phi0 = fTrackPar[3];
            const float r    = fECALInnerRadius;
            return this->helix_circle_intersections(phi, y0, z0, R, phi0, r);
        };

        ROOT::Math::WrappedFunction<std::function<float(float)>> helix_circle_intersections_wrapped_entry(helix_circle_intersections_to_wrap_entry);

        ROOT::Math::RootFinder rootFinderEntryBarrel(ROOT::Math::RootFinder::kBRENT);
        rootFinderEntryBarrel.SetFunction(helix_circle_intersections_wrapped_entry, fTrackPar[3], fPhiMax);

        bool  finderStatusBarrel = true;
        float EntryPointBarrel[3];
        float TrackLengthExtraBarrel = 9999.0;
        float root_entry_barrel = 0.0;

        try {
          rootFinderEntryBarrel.Solve();
          root_entry_barrel = rootFinderEntryBarrel.Root();
          
          helix(root_entry_barrel, fTrackPar, fTrackEnd, EntryPointBarrel, fT0);

          TrackLengthExtraBarrel = TMath::Abs((fTrackPar[3] - root_entry_barrel)/fTrackPar[2])*TMath::Sqrt(1+TMath::Power(TMath::Tan(fTrackPar[4]), 2));

        } catch (const std::exception &excpt) {
          finderStatusBarrel = false;
        }

        auto helix_x_intersections_to_wrap_entry = [this](float phi) {
            const float x0      = fTrackEnd[0];
            const float t0      = fT0;
            const float R       = 1/fTrackPar[2];
            const float lambda0 = fTrackPar[4];
            const float phi0    = fTrackPar[3];
            const float x       = ((x0 > 0.) ? 1. : ((x0 < 0.) ? -1. : 0.))*fECALStartX;
            return this->helix_x_intersections(phi, x0, t0, R, lambda0, phi0, x);
        };

        ROOT::Math::WrappedFunction<std::function<float(float)>> helix_x_intersections_wrapped_entry(helix_x_intersections_to_wrap_entry);

        ROOT::Math::RootFinder rootFinderEntryEndCap(ROOT::Math::RootFinder::kBRENT);
        rootFinderEntryEndCap.SetFunction(helix_x_intersections_wrapped_entry, fTrackPar[3], fPhiMax);

        bool  finderStatusEndCap = true;
        float EntryPointEndCap[3];
        float TrackLengthExtraEndCap = 9999.0;
        float root_entry_endcap = 0.0;

        try {
          rootFinderEntryEndCap.Solve();
          root_entry_endcap = rootFinderEntryEndCap.Root();
          
          helix(root_entry_endcap, fTrackPar, fTrackEnd, EntryPointEndCap, fT0);

          TrackLengthExtraEndCap = TMath::Abs((fTrackPar[3] - root_entry_endcap)/fTrackPar[2])*TMath::Sqrt(1+TMath::Power(TMath::Tan(fTrackPar[4]), 2));

        } catch (const std::exception &excpt) {
          finderStatusEndCap = false;
        }

        if ((finderStatusBarrel == false)&&(finderStatusEndCap == false)) {
          return false;
        }

        if ((TrackLengthExtraBarrel <= TrackLengthExtraEndCap)&&(finderStatusBarrel == true)) {
          // Using propagation to barrel
          fTrackLengthExtra = TrackLengthExtraBarrel;
          for (int i=0; i<3; ++i) fEntryPoint[i] = EntryPointBarrel[i];
          fEntryPhi = root_entry_barrel;
        } else if ((TrackLengthExtraEndCap < TrackLengthExtraBarrel)&&(finderStatusEndCap == true)) {
          // Using propagation to endcap
          fTrackLengthExtra = TrackLengthExtraEndCap;
          for (int i=0; i<3; ++i) fEntryPoint[i] = EntryPointEndCap[i];
          fEntryPhi = root_entry_endcap;
          fPropagateToEndCap = true;
        }

        fTrackLengthCorrected = fTrackLength+fTrackLengthExtra;

        return true;
      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::AddHits(std::vector<const rec::CaloHit*> hit_vec)
      {

        for (auto const& hit: hit_vec) {

          int DetID  = fFieldDecoder_ECAL->get(hit->CellID(), "system");
          int Layer  = hit->Layer();

          if (((DetID == 1)&&(Layer < 9))||((DetID == 2)&&(Layer < 7))) {
            if (fPropagateToEndCap) {
              PropagateHitEndCap(hit, Layer);
            } else {
              PropagateHitBarrel(hit, Layer);
            }
          }

        }

      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::PropagateHitEndCap(const rec::CaloHit* hit, int Layer)
      {

        TVector3 hitPosition(hit->Position());

        float xHit = hitPosition[0];
        float yHit = hitPosition[1];
        float zHit = hitPosition[2];
        float tHit = hit->Time().first;

        auto helix_x_intersections_to_wrap = [this, xHit](float phi) {
            const float x0      = fTrackEnd[0];
            const float t0      = fT0;
            const float R       = 1/fTrackPar[2];
            const float lambda0 = fTrackPar[4];
            const float phi0    = fTrackPar[3];
            const float x       = xHit;
            return this->helix_x_intersections(phi, x0, t0, R, lambda0, phi0, x);
        };

        ROOT::Math::WrappedFunction<std::function<float(float)>> helix_x_intersections_wrapped(helix_x_intersections_to_wrap);

        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
        rootFinder.SetFunction(helix_x_intersections_wrapped, fTrackPar[3], fPhiMax);

        float distance_hit_projected   = 0.0;
        float distance_projected_entry = 0.0;

        try {
          rootFinder.Solve();
          float root = rootFinder.Root();

          float projected[3];
          helix(root, fTrackPar, fTrackEnd, projected, fT0);

          distance_hit_projected   = std::hypot(xHit-projected[0], yHit-projected[1], zHit-projected[2]);
          distance_projected_entry = TMath::Abs((fEntryPhi - root)/fTrackPar[2])*TMath::Sqrt(1+TMath::Power(TMath::Tan(fTrackPar[4]), 2)); // compute arc length

          //float distance_hit_entry = std::hypot(distance_hit_projected, distance_projected_entry);

          if (fHitMap.find(Layer) == fHitMap.end()) {
              //fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_hit_entry, tHit);
              fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_projected_entry, tHit);
          } else {
              float old_distance_hit_projected = std::get<0>(fHitMap[Layer]);
              //if (distance_hit_projected <= old_distance_hit_projected) fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_hit_entry, tHit);
              if (distance_hit_projected <= old_distance_hit_projected) fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_projected_entry, tHit);
          }

        } catch (const std::exception &excpt) {
            distance_hit_projected   = -1.0;
            distance_projected_entry = -1.0;

        }

      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::PropagateHitBarrel(const rec::CaloHit* hit, int Layer)
      {

        TVector3 hitPosition(hit->Position());

        float xHit = hitPosition[0];
        float yHit = hitPosition[1];
        float zHit = hitPosition[2];
        float rHit = std::hypot(zHit-fTPCCent[2], yHit-fTPCCent[1]);
        float tHit = hit->Time().first;

        auto helix_circle_intersections_to_wrap = [this, rHit](float phi) {
            const float y0   = fTrackPar[0];
            const float z0   = fTrackPar[1];
            const float R    = 1/fTrackPar[2];
            const float phi0 = fTrackPar[3];
            const float r    = rHit;
            return this->helix_circle_intersections(phi, y0, z0, R, phi0, r);
        };

        ROOT::Math::WrappedFunction<std::function<float(float)>> helix_circle_intersections_wrapped(helix_circle_intersections_to_wrap);

        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
        rootFinder.SetFunction(helix_circle_intersections_wrapped, fTrackPar[3], fPhiMax);

        float distance_hit_projected   = 0.0;
        float distance_projected_entry = 0.0;

        try {
          rootFinder.Solve();
          float root = rootFinder.Root();

          float projected[3];
          helix(root, fTrackPar, fTrackEnd, projected, fT0);

          distance_hit_projected   = std::hypot(xHit-projected[0], yHit-projected[1], zHit-projected[2]);
          //distance_projected_entry = std::hypot(fEntryPoint[0]-projected[0], fEntryPoint[1]-projected[1], fEntryPoint[2]-projected[2]);
          distance_projected_entry = TMath::Abs((fEntryPhi - root)/fTrackPar[2])*TMath::Sqrt(1+TMath::Power(TMath::Tan(fTrackPar[4]), 2)); // compute arc length

          //float distance_hit_entry = std::hypot(distance_hit_projected, distance_projected_entry);

          if (fHitMap.find(Layer) == fHitMap.end()) {
              //fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_hit_entry, tHit);
              fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_projected_entry, tHit);
          } else {
              float old_distance_hit_projected = std::get<0>(fHitMap[Layer]);
              //if (distance_hit_projected <= old_distance_hit_projected) fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_hit_entry, tHit);
              if (distance_hit_projected <= old_distance_hit_projected) fHitMap[Layer] = std::make_tuple(distance_hit_projected, distance_projected_entry, tHit);
          }

        } catch (const std::exception &excpt) {
            distance_hit_projected   = -1.0;
            distance_projected_entry = -1.0;

        }

      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::ComputeArrivalTime()
      {

        if (fHitMap.size() == 0) {
          return; // stop if your HitMap is empty (no tile hits)
        }

        std::vector<float> distance_vec;
        std::vector<float> time_vec;
        std::vector<float> time_corrected_vec;
        for (auto& [key, value]: fHitMap) {

          distance_vec.push_back(std::get<1>(value)*10);
          time_vec.push_back(std::get<2>(value));
          time_corrected_vec.push_back(std::get<2>(value)-(std::get<1>(value)*10)/300);

        }

        if (fTimeMethod.compare("Earliest") == 0) {
          std::vector<float>::iterator time_min = std::min_element(time_vec.begin(), time_vec.end());
          float time_min_distance = distance_vec[std::distance(time_vec.begin(), time_min)];

          fTime = *time_min-time_min_distance/300;

        } else if (fTimeMethod.compare("Average") == 0) {
          fTime = std::accumulate(time_corrected_vec.begin(), time_corrected_vec.end(), 0.0) / time_corrected_vec.size();

        } else if (fTimeMethod.compare("Fit") == 0) {
          float distance_min = *std::min_element(distance_vec.begin(), distance_vec.end());
          float distance_max = *std::max_element(distance_vec.begin(), distance_vec.end());

          if (distance_vec.size() >= 2) {
            TGraph* graph = new TGraph(distance_vec.size(), &distance_vec[0], &time_vec[0]);
            TF1* linearFit = new TF1("linearFit", "pol1", distance_min, distance_max);
            graph->Fit(linearFit, "Q");

            fTime = linearFit->GetParameter(0);

          }
        }

        fBeta = beta(fTrackLengthCorrected, fTime);
        fMass = mass(fTrackMomentum, fTrackLengthCorrected, fTime);

        return;

      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::GetTime()
      {
        return fTime;
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::GetBeta()
      {
        return fBeta;
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::GetMass()
      {
        return (fBeta > 1.0) ? 0.0 : fMass;
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::helix_circle_intersections(float phi, float y0, float z0, float R, float phi0, float r)
      {
        return TMath::Power(y0-R*(TMath::Cos(phi)-TMath::Cos(phi0))-fTPCCent[1], 2)+TMath::Power(z0+R*(TMath::Sin(phi)-TMath::Sin(phi0))-fTPCCent[2], 2)-TMath::Power(r, 2);
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::helix_x_intersections(float phi, float x0, float t0, float R, float lambda0, float phi0, float x)
      {
        return x0 + ((x0 > 0.) ? 1. : ((x0 < 0.) ? -1. : 0.))*fDriftVelocity*t0 + R*TMath::Tan(lambda0)*(phi-phi0) - fTPCCent[0] - x;
      }

      //----------------------------------------------------------------------------
      void ECALToFAlg::helix(float phi, float *trackPar, float *trackEnd, float *projected, float t0)
      {
          
        projected[0] = trackEnd[0] + ((trackEnd[0] > 0.) ? 1. : ((trackEnd[0] < 0.) ? -1. : 0.))*fDriftVelocity*t0 + (1/trackPar[2])*TMath::Tan(trackPar[4])*(phi-trackPar[3]);
        projected[1] = trackPar[0] - (1/trackPar[2])*(TMath::Cos(phi) - TMath::Cos(trackPar[3]));
        projected[2] = trackPar[1] + (1/trackPar[2])*(TMath::Sin(phi) - TMath::Sin(trackPar[3]));
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::beta(float length, float time)
      {
        return length/(time*30); // length in cm and time in ns, so c = 30 cm/ns
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::mass(float momentum, float length, float time)
      {
        float b = beta(length, time);
        return momentum*TMath::Sqrt(1-TMath::Power(b, 2))/b;
      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::GetToFProtonScore()
      {

      if (fHitMap.size() == 0) {
          return -1.0; // stop if your HitMap is empty (no tile hits)
      }

      if (fBeta >= 1.0) {
        // if beta >= 1 then assume it's not a proton
        return 0.0;
      }

      if (fTrackMomentum >= 3.00) {
        // for momenta higher than 3.0 GeV we can't say nothing
        return 0.0;
      } else if (fTrackMomentum < 0.50) {
        // for momenta lower than 500 MeV we apply a simple cut at 0.6
        if (fBeta <= 0.6) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      float ProtonScore = 0.0;
      for (auto& [key, clf]: fScorerMap) {

          float p_min = std::stof(key.first);
          float p_max = std::stof(key.second);

          if ((fTrackMomentum >= p_min)&&(fTrackMomentum < p_max)) {
            // Apply corresponding probability calibration
            CalibratedToFScore calibration = fScorerMap[std::make_pair(key.first, key.second)];
            ProtonScore = Sigmoid(-(fBeta-calibration.beta_max_f1)*10, calibration.calibration_a, calibration.calibration_b);
            break;
          }
        }

        return ProtonScore;

      }

      //----------------------------------------------------------------------------
      float ECALToFAlg::Sigmoid(float x, float a, float b)
      {
        return 1.0/(1.0+TMath::Exp(a*x+b));
      }

    } // namespace alg
  } // namespace rec
} // namespace gar
