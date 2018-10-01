//
//  AuxDetAction.cxx
//
//  Created by Brian Rebel on 11/22/16.
//

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoNode.h"
#include "TGeoBBox.h"

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

#include "GArG4/AuxDetAction.h"
#include "GArG4/ParticleListAction.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar {

  namespace garg4 {

    //-------------------------------------------------------------
    // Constructor.
    AuxDetAction::AuxDetAction(CLHEP::HepRandomEngine*    engine,
      fhicl::ParameterSet const& pset)
      //: fEngine(engine)
      {
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
        fEnergyCut = pset.get<double>("EnergyCut");
        fVolumeName = pset.get<std::vector<std::string>>("AuxDetVolumeName");
        fMaterialMatchString = pset.get<std::string>("AuxDetMaterial");
        return;
      }

      //-------------------------------------------------------------
      void AuxDetAction::BeginOfEventAction(const G4Event*)
      {
        // Clear any previous information.
        // fAuxDetSimChannels.clear();
        fDeposits.clear();
      }

      //-------------------------------------------------------------
      void AuxDetAction::PreTrackingAction(const G4Track* track)
      {
      }

      //-------------------------------------------------------------
      void AuxDetAction::PostTrackingAction( const G4Track* /*track*/)
      {
      }

      //-------------------------------------------------------------
      // With every step, check and see if we are in an auxdet and save the energy
      void AuxDetAction::SteppingAction(const G4Step* step)
      {
        LOG_DEBUG("AuxDetAction")
        << "AuxDetAction::SteppingAction";

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
        std::string volname = track->GetVolume()->GetName();
        std::string VolumeName;

        // check the material
        auto pos = 0.5 * (start + stop);
        TGeoNode *node = geomanager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm...

        if(!node){
          LOG_WARNING("AuxDetAction")
          << "Node not found in "
          << pos.x() << " mm "
          << pos.y() << " mm "
          << pos.z() << " mm";
          return;
        }

        std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();

        if ( ! std::regex_match(volmaterial, std::regex(fMaterialMatchString)) ) return;

        // only worry about energy depositions larger than the minimum required
        if(step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fEnergyCut){

          LOG_WARNING("AuxDetAction")
          << "In volume "
          << volname
          << " Material is "
          << volmaterial
          << " step size is "
          << step->GetStepLength() / CLHEP::cm
          << " and deposited "
          << step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV
          << " GeV of energy with a minimum of "
          << fEnergyCut
          << " required.";

          // the step mid point is used for the position of the deposit
          auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
          step->GetPostStepPoint()->GetPosition() );

          int CaloID = -1;
          if( (strncmp(volname.c_str(), "IBStrip_L2_vol_PV", 17) == 0) )
          CaloID = 1;//Inner Barrel
          if( (strncmp(volname.c_str(), "OBStrip_L2_vol_PV", 17) == 0) )
          CaloID = 2;//Outer Barrel
          if( (strncmp(volname.c_str(), "IECLayer_L2_vol_PV", 18) == 0) )
          CaloID = 3;//Endcap

          // get the track id for this step
          auto trackID  = ParticleListAction::GetCurrentTrackID();

          fDeposits.emplace_back(trackID,
            step->GetPreStepPoint()->GetGlobalTime(),//get the time of the first subhit
            step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV,
            midPoint.x() / CLHEP::cm,
            midPoint.y() / CLHEP::cm,
            midPoint.z() / CLHEP::cm,
            CaloID, 0, 0);

          } // end if enough energy to worry about this step

        }// end of AuxDetAction::SteppingAction

        //------------------------------------------------------------------------------
        void AuxDetAction::EndOfEventAction(const G4Event*)
        {
          //sort per time the hits
          std::sort(fDeposits.begin(), fDeposits.end());
        }

      } // garg4

    } // gar
