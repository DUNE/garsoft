//
//  EnergyDepositAction.h
//
//  Created by Brian Rebel on 10/12/16.
//
/// This class implements the G4Base::UserAction interface in order to
/// accumulate a list of hits in the gaseous argon modeled by Geant4.
//

#ifndef GAR_EnergyDepositAction_h
#define GAR_EnergyDepositAction_h

// NuTools includes
#include "nutools/G4Base/UserAction.h"

// GArSoft includes
#include "Geometry/Geometry.h"
#include "SimulationDataProducts/EnergyDeposit.h"

//ART includes
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CLHEP/Random/RandGauss.h"

// Forward declarations.
class G4Event;
class G4Track;
class G4Step;
class G4EnergyLossForExtrapolator;

namespace gar {
  namespace garg4 {

    ///list of energy deposits from Geant4
    class EnergyDepositAction : public g4b::UserAction {

    public:

      // Standard constructors and destructors;
      EnergyDepositAction(CLHEP::HepRandomEngine*        engine,
        fhicl::ParameterSet     const& pset);
        virtual ~EnergyDepositAction();

        void reconfigure(fhicl::ParameterSet const& pset);

        // UserActions method that we'll override, to obtain access to
        // Geant4's particle tracks and trajectories.
        void BeginOfEventAction(const G4Event*);
        void EndOfEventAction  (const G4Event*);
        void PreTrackingAction (const G4Track*);
        void PostTrackingAction(const G4Track*);
        void SteppingAction    (const G4Step* );

        std::string GetVolumeName(const G4Track *track);

        //  Returns the EnergyDeposits accumulated during the current event.
        std::vector<gar::sdp::EnergyDeposit> const& EnergyDeposits() const { return fDeposits; }

      private:

        void AddEnergyDeposition(const G4Step* step);

        double                               fEnergyCut;  ///< The minimum energy in GeV for a deposit to
        ///< be included in the list.
        std::vector<std::string>             fVolumeName; ///< volume we will record energy depositions in
        std::string                          fMaterialMatchString; ///< Energy deposition will be recorded for materials that match this
        //unused CLHEP::HepRandomEngine*              fEngine;     ///< random number engine
        std::vector<gar::sdp::EnergyDeposit> fDeposits;   ///< energy fDeposits
      };

    } // garg4
  } // gar

  #endif /* GAR_EnergyDepositAction_h */
