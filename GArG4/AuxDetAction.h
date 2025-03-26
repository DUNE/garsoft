//
//  AuxDetAction.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/22/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef AuxDetAction_h
#define AuxDetAction_h

// NuTools includes
#include "nug4/G4Base/UserAction.h"

// GArSoft includes
#include "Geometry/ChannelMapAlgs/MinervaSegmentationAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/GeometryGAr.h"

#include "DetectorInfo/DetectorProperties.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"

#include "CLHEP/Random/RandGauss.h"

// Forward declarations.
class G4Event;
class G4Track;
class G4Step;
class G4EnergyLossForExtrapolator;

namespace gar {
  namespace garg4 {

    /// list of energy deposits from Geant4
    class AuxDetAction : public g4b::UserAction {

    public:
      // Standard constructors and destructors;
      AuxDetAction(CLHEP::HepRandomEngine* engine, fhicl::ParameterSet const& pset);
      virtual ~AuxDetAction();

      void reconfigure(fhicl::ParameterSet const& pset);

      // UserActions method that we'll override, to obtain access to
      // Geant4's particle tracks and trajectories.
      void BeginOfEventAction(const G4Event*);
      void EndOfEventAction(const G4Event*);
      void PreTrackingAction(const G4Track*);
      void PostTrackingAction(const G4Track*);
      void SteppingAction(const G4Step*);

      //  Returns the LArDeposit set accumulated during the current event.
      std::vector<gar::sdp::LArDeposit> const& LArDeposits() const { return fLArDeposits; }
      //  Returns the CaloDeposit set accumulated during the current event.
      std::vector<gar::sdp::CaloDeposit> const& CaloDeposits() const { return fECALDeposits; }
      //  Returns the CaloDeposit set accumulated during the current event.
      std::vector<gar::sdp::CaloDeposit> const& TrackerScDeposits() const
      {
        return fTrackerScDeposits;
      }
      //  Returns the CaloDeposit set accumulated during the current event.
      std::vector<gar::sdp::CaloDeposit> const& MuIDDeposits() const { return fMuIDDeposits; }

    private:
      void LArSteppingAction(const G4Step*);
      void ECALSteppingAction(const G4Step*);
      void TrackerScSteppingAction(const G4Step*);
      void MuIDSteppingAction(const G4Step*);
      void AddHits(const std::map<gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit>> m_hits,
                   std::vector<gar::sdp::CaloDeposit>& fDeposits);

      std::string GetVolumeName(const G4Track* track);

      unsigned int GetDetNumber(std::string volname);
      unsigned int GetStaveNumber(std::string volname);
      unsigned int GetModuleNumber(std::string volname);
      unsigned int GetLayerNumber(std::string volname);
      unsigned int GetSliceNumber(std::string volname);

      G4ThreeVector globalToLocal(const G4Step* step, const G4ThreeVector& glob);
      G4ThreeVector localToGlobal(const G4Step* step, const G4ThreeVector& loc);
      float GetStepEnergy(const G4Step* step, bool birks);
      float birksAttenuation(const G4Step* step);

      double
        fLArEnergyCut; ///< The minimum energy in GeV for a particle to be included in the list.
      std::string fLArMaterial;                ///< Material for the LArTPC
      std::vector<std::string> fLArVolumeName; ///< volume we will record energy depositions in

      std::string fECALMaterial;                ///< Material for the ECAL
      std::vector<std::string> fECALVolumeName; ///< volume we will record energy depositions in

      std::string fTrackerScMaterial; ///< Material for the TrackerSc (GArLite)
      std::vector<std::string>
        fTrackerScVolumeName; ///< volume we will record energy depositions in

      std::string fMuIDMaterial;                ///< Material for the MuID
      std::vector<std::string> fMuIDVolumeName; ///< volume we will record energy depositions in
      bool fApplyBirksLaw;

      std::vector<gar::sdp::LArDeposit> fLArDeposits; ///< energy fDeposits for the LArTPC

      std::map<gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit>> m_ECALDeposits;
      std::vector<gar::sdp::CaloDeposit> fECALDeposits; ///< energy fDeposits for the ECAL

      std::map<gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit>> m_TrackerScDeposits;
      std::vector<gar::sdp::CaloDeposit> fTrackerScDeposits; ///< energy fDeposits for the TrackerSc

      std::map<gar::raw::CellID_t, std::vector<gar::sdp::CaloDeposit>> m_MuIDDeposits;
      std::vector<gar::sdp::CaloDeposit> fMuIDDeposits; ///< energy fDeposits for the MuID

      const gar::geo::GeometryCore* fGeo;          ///< geometry information
      const detinfo::DetectorProperties* fDetProp; ///< detector properties
      const gar::geo::seg::MinervaSegmentationAlg* fMinervaSegAlg;
    };

  } // garg4
} // gar

#endif /* AuxDetAction_h */
