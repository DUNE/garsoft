//
//  EnergyDepositAction.cxx
//
//  Created by Brian Rebel on 10/12/16.
//

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"

// G4 includes
#include "Geant4/G4Event.hh"
#include "Geant4/G4Track.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4PrimaryParticle.hh"
#include "Geant4/G4DynamicParticle.hh"
#include "Geant4/G4VUserPrimaryParticleInformation.hh"
#include "Geant4/G4Step.hh"
#include "Geant4/G4StepPoint.hh"
#include "Geant4/G4VProcess.hh"

#include "Geometry/Geometry.h"
#include "GArG4/EnergyDepositAction.h"
#include "GArG4/ParticleListAction.h"

namespace gar {

  namespace garg4 {

    //-------------------------------------------------------------
    // Constructor.
    EnergyDepositAction::EnergyDepositAction(CLHEP::HepRandomEngine*    engine,
      fhicl::ParameterSet const& pset)
      //: fEngine(engine)
      {
        this->reconfigure(pset);
      }

      //-------------------------------------------------------------
      // Destructor.
      EnergyDepositAction::~EnergyDepositAction()
      {
        // Delete anything that we created with "new'.
      }

      //-------------------------------------------------------------
      void EnergyDepositAction::reconfigure(fhicl::ParameterSet const& pset)
      {
        fEnergyCut  = pset.get< double >("EnergyCut");
        fVolumeName = pset.get< std::vector<std::string> >("GArVolumeName");
        fMaterialMatchString = pset.get< std::string >("GArMaterialMatchString");

        std::cout << "EnergyDepositAction: Name of the Volumes to track for the GArTPC" << std::endl;
        for(unsigned int i = 0; i < fVolumeName.size(); i++) std::cout << fVolumeName.at(i) << " ";
        std::cout << std::endl;

        return;
      }

      //-------------------------------------------------------------
      void EnergyDepositAction::BeginOfEventAction(const G4Event*)
      {
        // Clear any previous information.
        fDeposits.clear();
      }

      //-------------------------------------------------------------
      void EnergyDepositAction::PreTrackingAction(const G4Track* /*track*/)
      {
      }

      //-------------------------------------------------------------
      void EnergyDepositAction::PostTrackingAction( const G4Track* /*track*/)
      {
      }

      //-------------------------------------------------------------
      // With every step, handle energy deposition in gaseous argon
      void EnergyDepositAction::SteppingAction(const G4Step* step)
      {
        LOG_DEBUG("EnergyDepositAction")
        << "EnergyDepositAction::SteppingAction";

        art::ServiceHandle<geo::Geometry> geo;
        TGeoManager *geomanager = geo->ROOTGeoManager();

        // Get the pointer to the track
        G4Track *track = step->GetTrack();

        const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
        const CLHEP::Hep3Vector &stop  = track->GetPosition();

        // If it's a null step, don't use it.
        if(start == stop) return;
        if(step->GetTotalEnergyDeposit() == 0) return;

        // check that we are in the correct material to record a hit
        std::string VolumeName   = this->GetVolumeName(track);

        if( std::find( fVolumeName.begin(), fVolumeName.end(), VolumeName ) == fVolumeName.end() )
        return;

        // check the material
        auto pos = 0.5 * (start + stop);
        TGeoNode *node = geomanager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm

        if(!node){
          LOG_DEBUG("EnergyDepositAction")
          << "Node not found in "
          << pos.x() << " mm "
          << pos.y() << " mm "
          << pos.z() << " mm";
          return;
        }

        std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
        if ( ! std::regex_match(volmaterial, std::regex(fMaterialMatchString)) ) return;

        // only worry about energy depositions larger than the minimum required
        if(step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fEnergyCut ){

          LOG_DEBUG("EnergyDepositAction")
          << "In volume "
          << VolumeName
          << " step size is "
          << step->GetStepLength() / CLHEP::cm
          << " cm and deposited "
          << step->GetTotalEnergyDeposit()
          << " MeV of energy with a minimum of "
          << fEnergyCut
          << " required.";

          // save this deposition
          this->AddEnergyDeposition(step);

        } // end if enough energy to worry about this step

      }// end of EnergyDepositAction::SteppingAction

      //--------------------------------------------------------------------------
      void EnergyDepositAction::AddEnergyDeposition(const G4Step* step)
      {
        // get the track id for this step
        auto trackID  = ParticleListAction::GetCurrentTrackID();

        // first check that we are not dealing with ionization, ie delta rays
        // if we have one of those, set the trackID to be the ID of the parent
        //if(trackID < 0)
        //  LOG_VERBATIM("EnergyDepositAction")
        //  << "TrackID "
        //  << trackID
        //  << " was created by process "
        //  << step->GetTrack()->GetCreatorProcess()->GetProcessName();
        //else
        //  LOG_VERBATIM("EnergyDepositAction")
        //  << "TrackID "
        //  << trackID
        //  << " is a primary particle ";

        // the step mid point is used for the position of the deposit
        auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
        step->GetPostStepPoint()->GetPosition() );
        float time = step->GetPreStepPoint()->GetGlobalTime();

        fDeposits.emplace_back(trackID,
          time,
          step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV,
          midPoint.x() / CLHEP::cm,
          midPoint.y() / CLHEP::cm,
          midPoint.z() / CLHEP::cm,
          step->GetStepLength() / CLHEP::cm,
          (trackID > 0));

          return;
        }

        //--------------------------------------------------------------------------
        void EnergyDepositAction::EndOfEventAction(const G4Event*)
        {
          // sort the EnergyDeposit lists in each EnergyDeposits object
          std::sort(fDeposits.begin(), fDeposits.end());
        }

        //--------------------------------------------------------------------------
        std::string EnergyDepositAction::GetVolumeName(const G4Track *track)
        {
          std::string VolName = track->GetVolume()->GetName();
          VolName.erase(VolName.length()-3, 3);
          return VolName;
        }


      } // garg4

    } // gar
