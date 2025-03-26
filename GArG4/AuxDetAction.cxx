//
//  AuxDetAction.cxx
//
//  Created by Brian Rebel on 11/22/16.
//

#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <regex>

#include "TGeoBBox.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"

// G4 includes
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4EmSaturation.hh"
#include "Geant4/G4Event.hh"
#include "Geant4/G4MaterialCutsCouple.hh"
#include "Geant4/G4NavigationHistory.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PrimaryParticle.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4VProcess.hh"
#include "Geant4/G4VUserPrimaryParticleInformation.hh"
#include "Geant4/G4Version.hh"

#include "GArG4/AuxDetAction.h"
#include "GArG4/ParticleListAction.h"

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
namespace gar {

  namespace garg4 {

    //-------------------------------------------------------------
    // Constructor.
    AuxDetAction::AuxDetAction(CLHEP::HepRandomEngine*, //  engine,
                               fhicl::ParameterSet const& pset)
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

      this->reconfigure(pset);
    }

    //-------------------------------------------------------------
    // Destructor.
    AuxDetAction::~AuxDetAction()
    {
      // Delete anything that we created with "new'.
    }

    //-------------------------------------------------------------
    void AuxDetAction::reconfigure(fhicl::ParameterSet const& pset)
    {
      fLArEnergyCut = pset.get<double>("LArEnergyCut", 0.);
      fLArVolumeName = pset.get<std::vector<std::string>>("LArVolumeName");
      fLArMaterial = pset.get<std::string>("LArMaterial", "LAr");

      fECALMaterial = pset.get<std::string>("ECALMaterial", "Scintillator");
      fECALVolumeName = pset.get<std::vector<std::string>>("ECALVolumeName");

      fTrackerScMaterial = pset.get<std::string>("TrackerScMaterial", "Scintillator");
      fTrackerScVolumeName = pset.get<std::vector<std::string>>("TrackerScVolumeName");

      fMuIDMaterial = pset.get<std::string>("MuIDMaterial", "Scintillator");
      fMuIDVolumeName = pset.get<std::vector<std::string>>("MuIDVolumeName");

      fApplyBirksLaw = pset.get<bool>("ApplyBirksLaw", true);

      // std::cout << "AuxDetAction: Name of the Volumes to track for the LArTPC" << std::endl;
      // for(unsigned int i = 0; i < fLArVolumeName.size(); i++) std::cout << fLArVolumeName.at(i) << " ";
      // std::cout << std::endl;

      fMinervaSegAlg = (gar::geo::seg::MinervaSegmentationAlg*)fGeo->MinervaSegmentationAlg();

      return;
    }

    //-------------------------------------------------------------
    void AuxDetAction::BeginOfEventAction(const G4Event*)
    {
      // Clear any previous information.
      m_ECALDeposits.clear();
      fECALDeposits.clear();

      m_TrackerScDeposits.clear();
      fTrackerScDeposits.clear();

      m_MuIDDeposits.clear();
      fMuIDDeposits.clear();

      fLArDeposits.clear();
    }

    //-------------------------------------------------------------
    void AuxDetAction::PreTrackingAction(const G4Track* /* track */) {}

    //-------------------------------------------------------------
    void AuxDetAction::PostTrackingAction(const G4Track* /*track*/) {}

    //-------------------------------------------------------------
    // With every step, check and see if we are in an auxdet and save the energy
    void AuxDetAction::SteppingAction(const G4Step* step)
    {
      MF_LOG_DEBUG("AuxDetAction") << "AuxDetAction::SteppingAction";

      if (fGeo->HasLArTPCDetector()) this->LArSteppingAction(step);

      if (fGeo->HasECALDetector()) this->ECALSteppingAction(step);

      if (fGeo->HasTrackerScDetector()) this->TrackerScSteppingAction(step);

      if (fGeo->HasMuonDetector()) this->MuIDSteppingAction(step);
    } // end of AuxDetAction::SteppingAction

    //------------------------------------------------------------------------------
    void AuxDetAction::EndOfEventAction(const G4Event*)
    {
      //sort per time the hits
      if (fGeo->HasLArTPCDetector()) { std::sort(fLArDeposits.begin(), fLArDeposits.end()); }

      if (fGeo->HasECALDetector()) {
        this->AddHits(m_ECALDeposits, fECALDeposits);
        std::sort(fECALDeposits.begin(), fECALDeposits.end());
      }

      if (fGeo->HasTrackerScDetector() && nullptr != fMinervaSegAlg) {
        fMinervaSegAlg->AddHitsMinerva(m_TrackerScDeposits, fTrackerScDeposits);
        std::sort(fTrackerScDeposits.begin(), fTrackerScDeposits.end());
      }

      if (fGeo->HasMuonDetector()) {
        this->AddHits(m_MuIDDeposits, fMuIDDeposits);
        std::sort(fMuIDDeposits.begin(), fMuIDDeposits.end());
      }
    }

    //------------------------------------------------------------------------------
    void AuxDetAction::LArSteppingAction(const G4Step* step)
    {
      MF_LOG_DEBUG("AuxDetAction") << "AuxDetAction::LArSteppingAction";

      // Get the pointer to the track
      G4Track* track = step->GetTrack();

      const CLHEP::Hep3Vector& start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector& stop = track->GetPosition();

      // If it's a null step, don't use it.
      if (start == stop) return;
      if (step->GetTotalEnergyDeposit() == 0) return;

      // check that we are in the correct material to record a hit
      const std::string VolumeName = this->GetVolumeName(track);
      if (std::find_if(fLArVolumeName.begin(), fLArVolumeName.end(), [&VolumeName](std::string s) {
            return VolumeName.find(s) != std::string::npos;
          }) == fLArVolumeName.end())
        return;

      // check the material
      auto pos = 0.5 * (start + stop);
      TGeoNode* node = fGeo->FindNode(
        pos.x() / CLHEP::cm, pos.y() / CLHEP::cm, pos.z() / CLHEP::cm); //Node in cm...

      if (nullptr == node) {
        MF_LOG_DEBUG("AuxDetAction::LArSteppingAction")
          << "Node not found in " << pos.x() << " mm " << pos.y() << " mm " << pos.z() << " mm";
        return;
      }

      std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
      if (!std::regex_match(volmaterial, std::regex(fLArMaterial))) return;

      // only worry about energy depositions larger than the minimum required
      if (step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fLArEnergyCut) {

        MF_LOG_DEBUG("AuxDetAction::LArSteppingAction")
          << "In volume " << VolumeName << " Material is " << volmaterial << " step size is "
          << step->GetStepLength() / CLHEP::cm << " cm and deposited "
          << step->GetTotalEnergyDeposit() << " MeV of energy with a minimum of " << fLArEnergyCut
          << " required.";

        // the step mid point is used for the position of the deposit
        auto midPoint =
          0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());

        // get the track id for this step
        auto trackID = ParticleListAction::GetCurrentTrackID();
        float time = step->GetPreStepPoint()->GetGlobalTime();

        fLArDeposits.emplace_back(trackID,
                                  time, //get the time of the first subhit
                                  step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV,
                                  midPoint.x() / CLHEP::cm,
                                  midPoint.y() / CLHEP::cm,
                                  midPoint.z() / CLHEP::cm,
                                  step->GetStepLength() / CLHEP::cm,
                                  (trackID > 0));

      } // end if enough energy to worry about this step
    }

    //------------------------------------------------------------------------------
    void AuxDetAction::ECALSteppingAction(const G4Step* step)
    {
      MF_LOG_DEBUG("AuxDetAction") << "AuxDetAction::ECALSteppingAction";

      // Get the pointer to the track
      G4Track* track = step->GetTrack();

      const CLHEP::Hep3Vector& start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector& stop = track->GetPosition();

      // If it's a null step, don't use it.
      if (start == stop) return;
      if (step->GetTotalEnergyDeposit() == 0) return;

      // check that we are in the correct material to record a hit
      const std::string VolumeName = this->GetVolumeName(track);

      if (std::find_if(
            fECALVolumeName.begin(), fECALVolumeName.end(), [&VolumeName](std::string s) {
              return VolumeName.find(s) != std::string::npos;
            }) == fECALVolumeName.end())
        return;

      // check the material
      auto pos = 0.5 * (start + stop);
      TGeoNode* node = fGeo->FindNode(
        pos.x() / CLHEP::cm, pos.y() / CLHEP::cm, pos.z() / CLHEP::cm); //Node in cm...

      if (nullptr == node) {
        MF_LOG_DEBUG("AuxDetAction::ECALSteppingAction")
          << "Node not found in " << pos.x() << " mm " << pos.y() << " mm " << pos.z() << " mm";
        return;
      }

      std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
      if (!std::regex_match(volmaterial, std::regex(fECALMaterial))) return;

      // only worry about energy depositions larger than the minimum required
      MF_LOG_DEBUG("AuxDetAction::ECALSteppingAction")
        << "In volume " << VolumeName << " Material is " << volmaterial << " step size is "
        << step->GetStepLength() / CLHEP::cm << " cm and deposited "
        << step->GetTotalEnergyDeposit() << " MeV of energy";

      // the step mid point is used for the position of the deposit
      G4ThreeVector G4Global =
        0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());

      //get layer and slice number from the volume name (easier for segmentation later)
      unsigned int layer = this->GetLayerNumber(VolumeName);   //get layer number
      unsigned int slice = this->GetSliceNumber(VolumeName);   // get slice number
      unsigned int det_id = this->GetDetNumber(VolumeName);    // 1 == Barrel, 2 = Endcap
      unsigned int stave = this->GetStaveNumber(VolumeName);   //get the stave number
      unsigned int module = this->GetModuleNumber(VolumeName); //get the module number

      MF_LOG_DEBUG("AuxDetAction::ECALSteppingAction")
        << "det_id " << det_id << " stave " << stave << " module " << module << " layer " << layer
        << " slice " << slice;

      //Transform from global coordinates to local coordinates in mm
      G4ThreeVector G4Local = this->globalToLocal(step, G4Global);
      //Transform in cm
      std::array<double, 3> G4Localcm = {
        G4Local.x() / CLHEP::cm, G4Local.y() / CLHEP::cm, G4Local.z() / CLHEP::cm};
      //Get cellID
      raw::CellID_t cellID = fGeo->GetCellID(
        node, det_id, stave, module, layer, slice, G4Localcm); //encoding the cellID on 64 bits

      //Don't correct for the strip or tile position yet
      double G4Pos[3] = {0., 0., 0.}; // in cm
      G4Pos[0] = G4Global.x() / CLHEP::cm;
      G4Pos[1] = G4Global.y() / CLHEP::cm;
      G4Pos[2] = G4Global.z() / CLHEP::cm;

      // get the track id for this step
      auto trackID = ParticleListAction::GetCurrentTrackID();
      float time = step->GetPreStepPoint()->GetGlobalTime();

      float edep = this->GetStepEnergy(step, true) * CLHEP::MeV / CLHEP::GeV;

      MF_LOG_DEBUG("AuxDetAction::ECALSteppingAction")
        << "Energy deposited " << step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV
        << " GeV in cellID " << cellID << " after Birks " << edep << " GeV";

      //Create a map of cellID to SimHit
      gar::sdp::CaloDeposit hit(trackID, time, edep, G4Pos, cellID);
      if (m_ECALDeposits.find(cellID) != m_ECALDeposits.end())
        m_ECALDeposits[cellID].push_back(hit);
      else {
        std::vector<gar::sdp::CaloDeposit> vechit;
        vechit.push_back(hit);
        m_ECALDeposits.emplace(cellID, vechit);
      }
    }

    //------------------------------------------------------------------------------
    void AuxDetAction::TrackerScSteppingAction(const G4Step* step)
    {
      MF_LOG_DEBUG("AuxDetAction") << "AuxDetAction::TrackerScSteppingAction";

      // Get the pointer to the track
      G4Track* track = step->GetTrack();

      const CLHEP::Hep3Vector& start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector& stop = track->GetPosition();

      // If it's a null step, don't use it.
      if (start == stop) return;
      if (step->GetTotalEnergyDeposit() == 0) return;

      // check that we are in the correct material to record a hit
      const std::string VolumeName = this->GetVolumeName(track);
      if (std::find_if(
            fTrackerScVolumeName.begin(), fTrackerScVolumeName.end(), [&VolumeName](std::string s) {
              return VolumeName.find(s) != std::string::npos;
            }) == fTrackerScVolumeName.end())
        return;

      // check the material
      auto pos = 0.5 * (start + stop);
      TGeoNode* node = fGeo->FindNode(
        pos.x() / CLHEP::cm, pos.y() / CLHEP::cm, pos.z() / CLHEP::cm); //Node in cm...

      if (nullptr == node) {
        MF_LOG_DEBUG("AuxDetAction::TrackerScSteppingAction")
          << "Node not found in " << pos.x() / CLHEP::cm << " cm " << pos.y() / CLHEP::cm << " cm "
          << pos.z() / CLHEP::cm << " cm in volume " << VolumeName;
        return;
      }

      std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
      if (!std::regex_match(volmaterial, std::regex(fTrackerScMaterial))) return;

      // only worry about energy depositions larger than the minimum required
      MF_LOG_DEBUG("AuxDetAction::TrackerScSteppingAction")
        << "In volume " << VolumeName << " Material is " << volmaterial << " step size is "
        << step->GetStepLength() / CLHEP::cm << " cm and deposited "
        << step->GetTotalEnergyDeposit() << " MeV of energy";

      // the step mid point is used for the position of the deposit
      G4ThreeVector G4Global =
        0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());

      //get layer and slice number from the volume name (easier for segmentation later)
      unsigned int layer = this->GetLayerNumber(VolumeName); //get layer number
      unsigned int slice = this->GetSliceNumber(VolumeName); // get slice number
      unsigned int det_id = 3;

      MF_LOG_DEBUG("AuxDetAction::TrackerScSteppingAction")
        << "det_id " << det_id << " layer " << layer << " slice " << slice;

      //Transform from global coordinates to local coordinates in mm
      G4ThreeVector G4Local = this->globalToLocal(step, G4Global);
      //Transform in cm
      std::array<double, 3> G4Localcm = {
        G4Local.x() / CLHEP::cm, G4Local.y() / CLHEP::cm, G4Local.z() / CLHEP::cm};
      //Get cellID
      raw::CellID_t cellID = fGeo->GetCellID(
        node, det_id, 0, 0, layer, slice, G4Localcm); //encoding the cellID on 64 bits

      //Don't correct for the strip or tile position yet
      double G4Pos[3] = {0., 0., 0.}; // in cm
      G4Pos[0] = G4Global.x() / CLHEP::cm;
      G4Pos[1] = G4Global.y() / CLHEP::cm;
      G4Pos[2] = G4Global.z() / CLHEP::cm;

      // get the track id for this step
      auto trackID = ParticleListAction::GetCurrentTrackID();
      float time = step->GetPreStepPoint()->GetGlobalTime();
      float edep = this->GetStepEnergy(step, true) * CLHEP::MeV / CLHEP::GeV;

      MF_LOG_DEBUG("AuxDetAction::TrackerScSteppingAction")
        << "Energy deposited " << step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV
        << " GeV in cellID " << cellID << " after Birks " << edep << " GeV";

      //Create a map of cellID to SimHit
      gar::sdp::CaloDeposit hit(trackID, time, edep, G4Pos, cellID);
      if (m_TrackerScDeposits.find(cellID) != m_TrackerScDeposits.end())
        m_TrackerScDeposits[cellID].push_back(hit);
      else {
        std::vector<gar::sdp::CaloDeposit> vechit;
        vechit.push_back(hit);
        m_TrackerScDeposits.emplace(cellID, vechit);
      }
    }

    //------------------------------------------------------------------------------
    void AuxDetAction::MuIDSteppingAction(const G4Step* step)
    {
      MF_LOG_DEBUG("AuxDetAction") << "AuxDetAction::MuIDSteppingAction";

      // Get the pointer to the track
      G4Track* track = step->GetTrack();

      const CLHEP::Hep3Vector& start = step->GetPreStepPoint()->GetPosition();
      const CLHEP::Hep3Vector& stop = track->GetPosition();

      // If it's a null step, don't use it.
      if (start == stop) return;
      if (step->GetTotalEnergyDeposit() == 0) return;

      // check that we are in the correct material to record a hit
      const std::string VolumeName = this->GetVolumeName(track);

      if (std::find_if(
            fMuIDVolumeName.begin(), fMuIDVolumeName.end(), [&VolumeName](std::string s) {
              return VolumeName.find(s) != std::string::npos;
            }) == fMuIDVolumeName.end())
        return;

      // check the material
      auto pos = 0.5 * (start + stop);
      TGeoNode* node = fGeo->FindNode(
        pos.x() / CLHEP::cm, pos.y() / CLHEP::cm, pos.z() / CLHEP::cm); //Node in cm...

      if (nullptr == node) {
        MF_LOG_DEBUG("AuxDetAction::MuIDSteppingAction")
          << "Node not found in " << pos.x() << " mm " << pos.y() << " mm " << pos.z() << " mm";
        return;
      }

      std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
      if (!std::regex_match(volmaterial, std::regex(fMuIDMaterial))) return;

      // only worry about energy depositions larger than the minimum required
      MF_LOG_DEBUG("AuxDetAction::MuIDSteppingAction")
        << "In volume " << VolumeName << " Material is " << volmaterial << " step size is "
        << step->GetStepLength() / CLHEP::cm << " cm and deposited "
        << step->GetTotalEnergyDeposit() << " MeV of energy";

      // the step mid point is used for the position of the deposit
      G4ThreeVector G4Global =
        0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());

      //get layer and slice number from the volume name (easier for segmentation later)
      unsigned int layer = this->GetLayerNumber(VolumeName);   //get layer number
      unsigned int slice = this->GetSliceNumber(VolumeName);   // get slice number
      unsigned int det_id = this->GetDetNumber(VolumeName);    // 1 == Barrel, 2 = Endcap
      unsigned int stave = this->GetStaveNumber(VolumeName);   //get the stave number
      unsigned int module = this->GetModuleNumber(VolumeName); //get the module number

      MF_LOG_DEBUG("AuxDetAction::MuIDSteppingAction")
        << "det_id " << det_id << " stave " << stave << " module " << module << " layer " << layer
        << " slice " << slice;

      //Transform from global coordinates to local coordinates in mm
      G4ThreeVector G4Local = this->globalToLocal(step, G4Global);
      //Transform in cm
      std::array<double, 3> G4Localcm = {
        G4Local.x() / CLHEP::cm, G4Local.y() / CLHEP::cm, G4Local.z() / CLHEP::cm};
      //Get cellID
      raw::CellID_t cellID = fGeo->GetCellID(
        node, det_id, stave, module, layer, slice, G4Localcm); //encoding the cellID on 64 bits

      //Don't correct for the strip or tile position yet
      double G4Pos[3] = {0., 0., 0.}; // in cm
      G4Pos[0] = G4Global.x() / CLHEP::cm;
      G4Pos[1] = G4Global.y() / CLHEP::cm;
      G4Pos[2] = G4Global.z() / CLHEP::cm;

      // get the track id for this step
      auto trackID = ParticleListAction::GetCurrentTrackID();
      float time = step->GetPreStepPoint()->GetGlobalTime();

      float edep = this->GetStepEnergy(step, true) * CLHEP::MeV / CLHEP::GeV;

      MF_LOG_DEBUG("AuxDetAction::MuIDSteppingAction")
        << "Energy deposited " << step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV
        << " GeV in cellID " << cellID << " after Birks " << edep << " GeV";

      //Create a map of cellID to SimHit
      gar::sdp::CaloDeposit hit(trackID, time, edep, G4Pos, cellID);
      if (m_MuIDDeposits.find(cellID) != m_MuIDDeposits.end())
        m_MuIDDeposits[cellID].push_back(hit);
      else {
        std::vector<gar::sdp::CaloDeposit> vechit;
        vechit.push_back(hit);
        m_MuIDDeposits.emplace(cellID, vechit);
      }
    }

    //------------------------------------------------------------------------------
    std::string AuxDetAction::GetVolumeName(const G4Track* track)
    {
      std::string VolName = track->GetVolume()->GetName();
      VolName.erase(VolName.length() - 3, 3);
      return VolName;
    }

    //------------------------------------------------------------------------------
    unsigned int AuxDetAction::GetLayerNumber(std::string volname)
    {
      return std::atoi((volname.substr(volname.find("layer_") + 6, 2)).c_str());
    }

    //------------------------------------------------------------------------------
    unsigned int AuxDetAction::GetSliceNumber(std::string volname)
    {
      return std::atoi((volname.substr(volname.find("slice") + 5, 1)).c_str());
    }

    //------------------------------------------------------------------------------
    unsigned int AuxDetAction::GetDetNumber(std::string volname)
    {
      unsigned int det_id = 0;

      if (volname.find("ECal") != std::string::npos) {
        if (volname.find("Barrel") != std::string::npos) det_id = 1;
        if (volname.find("Endcap") != std::string::npos) det_id = 2;
      }

      if (volname.find("Yoke") != std::string::npos) {
        if (volname.find("Barrel") != std::string::npos) det_id = 4;
        if (volname.find("Endcap") != std::string::npos) det_id = 4;
      }

      return det_id;
    }

    //------------------------------------------------------------------------------
    unsigned int AuxDetAction::GetStaveNumber(std::string volname)
    {
      return std::atoi((volname.substr(volname.find("_stave") + 6, 2)).c_str());
    }

    //------------------------------------------------------------------------------
    unsigned int AuxDetAction::GetModuleNumber(std::string volname)
    {
      return std::atoi((volname.substr(volname.find("_module") + 7, 2)).c_str());
    }

    //------------------------------------------------------------------------------
    G4ThreeVector AuxDetAction::globalToLocal(const G4Step* step, const G4ThreeVector& glob)
    {
      return step->GetPreStepPoint()
        ->GetTouchable()
        ->GetHistory()
        ->GetTopTransform()
        .TransformPoint(glob);
    }

    //------------------------------------------------------------------------------
    G4ThreeVector AuxDetAction::localToGlobal(const G4Step* step, const G4ThreeVector& loc)
    {
      return step->GetPreStepPoint()
        ->GetTouchable()
        ->GetHistory()
        ->GetTopTransform()
        .Inverse()
        .TransformPoint(loc);
    }

    //------------------------------------------------------------------------------
    float AuxDetAction::GetStepEnergy(const G4Step* step, bool birks)
    {
      if (birks) { return this->birksAttenuation(step); }
      else {
        return step->GetTotalEnergyDeposit();
      }
    }

    //------------------------------------------------------------------------------
    float AuxDetAction::birksAttenuation(const G4Step* step)
    {
#if G4VERSION_NUMBER >= 1001
      static G4EmSaturation s_emSaturation(1);
#else
      static G4EmSaturation s_emSaturation();
      s_emSaturation.SetVerbose(1);
#endif

#if G4VERSION_NUMBER >= 1030
      static bool s_initialised = false;
      if (not s_initialised) {
        s_emSaturation.InitialiseG4Saturation();
        s_initialised = true;
      }
#endif

      double energyDeposition = step->GetTotalEnergyDeposit();
      double length = step->GetStepLength();
      double niel = step->GetNonIonizingEnergyDeposit();
      const G4Track* trk = step->GetTrack();
      const G4ParticleDefinition* particle = trk->GetDefinition();
      const G4MaterialCutsCouple* couple = trk->GetMaterialCutsCouple();
      double engyVis =
        s_emSaturation.VisibleEnergyDeposition(particle, couple, length, energyDeposition, niel);
      return engyVis;
    }

    //------------------------------------------------------------------------------
    void AuxDetAction::AddHits(
      const std::map<gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit>> m_hits,
      std::vector<gar::sdp::CaloDeposit>& fDeposits)
    {
      //Loop over the hits in the map and add them together
      for (auto const& it : m_hits) {

        raw::CellID_t cellID = it.first;
        std::vector<gar::sdp::CaloDeposit> vechit = it.second;
        std::sort(vechit.begin(), vechit.end()); //sort per time

        float esum = 0.;
        float time = vechit.at(0).Time();
        int trackID = vechit.at(0).TrackID();
        double pos[3] = {vechit.at(0).X(), vechit.at(0).Y(), vechit.at(0).Z()};

        for (auto const& hit : vechit) {
          esum += hit.Energy();
        }

        fDeposits.emplace_back(trackID, time, esum, pos, cellID);
      }
    }
  } // garg4

} // gar
