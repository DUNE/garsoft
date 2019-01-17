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
#include "nutools/G4Base/UserAction.h"

// GArSoft includes
#include "Geometry/GeometryCore.h"

#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"
#include "DetectorInfo/DetectorProperties.h"

#include "GArG4EmSaturation.h"

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
      AuxDetAction(CLHEP::HepRandomEngine*        engine,
        fhicl::ParameterSet     const& pset);
        virtual ~AuxDetAction();

        void reconfigure(fhicl::ParameterSet const& pset);

        // UserActions method that we'll override, to obtain access to
        // Geant4's particle tracks and trajectories.
        void BeginOfEventAction(const G4Event*);
        void EndOfEventAction  (const G4Event*);
        void PreTrackingAction (const G4Track*);
        void PostTrackingAction(const G4Track*);
        void SteppingAction    (const G4Step* );

        void LArSteppingAction(const G4Step* );
        void ECALSteppingAction(const G4Step* );

        std::string GetVolumeName(const G4Track *track);

        unsigned int GetDetNumber(std::string volname);
        unsigned int GetStaveNumber(std::string volname);
        unsigned int GetModuleNumber(std::string volname);
        unsigned int GetLayerNumber(std::string volname);
        unsigned int GetSliceNumber(std::string volname);

        G4ThreeVector globalToLocal(const G4Step* step, const G4ThreeVector& glob);
        G4ThreeVector localToGlobal(const G4Step* step, const ROOT::Math::XYZVector& loc);
        G4ThreeVector localToGlobal(const G4Step* step, const G4ThreeVector& loc);

        //  Returns the CaloDeposit set accumulated during the current event.
        std::vector<gar::sdp::CaloDeposit> const& CaloDeposits() const { return fECALDeposits; }

        //  Returns the LArDeposit set accumulated during the current event.
        std::vector<gar::sdp::LArDeposit> const& LArDeposits() const { return fLArDeposits; }

      private:

        std::string fECALMaterial;                             ///< Material for the ECAL
        std::string fECALEncoding_str;                         ///< Encoding for the ECAL cells

        double                             fLArEnergyCut;     ///< The minimum energy in GeV for a particle to be included in the list.
        std::string fLArMaterial;                             ///< Material for the LArTPC
        std::vector<std::string>           fLArVolumeName;    ///< volume we will record energy depositions in

        std::vector<gar::sdp::CaloDeposit> fECALDeposits;          ///< energy fDeposits for the ECAL
        std::vector<gar::sdp::LArDeposit> fLArDeposits;          ///< energy fDeposits for the LArTPC

        const gar::geo::GeometryCore*      fGeo;               ///< geometry information
        TGeoManager*                       fGeoManager;

        GArG4EmSaturation           fGArG4EmSaturation;     ///< Determines the visible energy (after Birks suppression) Modified to include Birks-Chou

        const detinfo::DetectorProperties*       fDetProp;  ///< detector properties
      };

    } // garg4
  } // gar


  #endif /* AuxDetAction_h */
