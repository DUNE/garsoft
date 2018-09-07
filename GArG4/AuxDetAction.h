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
// #include "SimulationDataProducts/AuxDetSimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"

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

        //  Returns the AuxDetSimChannel set accumulated during the current event.
        // std::set<gar::sdp::AuxDetSimChannel> const& AuxDetSimChannels() const { return fAuxDetSimChannels; }
        std::vector<gar::sdp::CaloDeposit> const& CaloDeposits() const { return fDeposits; }

      private:

        // std::set<gar::sdp::AuxDetSimChannel> fAuxDetSimChannels; ///< The accumulated information for hits in the event.
        double                               fEnergyCut;         ///< The minimum energy in GeV for a particle to
        ///< be included in the list.
        std::string fMaterialMatchString;                        ///< Material for the AuxDet
        std::vector<std::string>             fVolumeName; ///< volume we will record energy depositions in
        //unused      CLHEP::HepRandomEngine*              fEngine;            ///< random number engine
        std::vector<gar::sdp::CaloDeposit> fDeposits;   ///< energy fDeposits
      };

    } // garg4
  } // gar


  #endif /* AuxDetAction_h */
