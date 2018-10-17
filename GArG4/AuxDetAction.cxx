//
//  AuxDetAction.cxx
//
//  Created by Brian Rebel on 11/22/16.
//

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>

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

#include "Geant4/G4MaterialCutsCouple.hh"

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
      {
        fGeo = gar::providerFrom<geo::Geometry>();
        fGeoManager = fGeo->ROOTGeoManager();
        fEmSaturation = new G4EmSaturation(0);

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
        fECALEnergyCut = pset.get<double>("ECALEnergyCut");
        fECALVolumeName = pset.get<std::vector<std::string>>("ECALVolumeName");
        fECALMaterial = pset.get<std::string>("ECALMaterial");

        fLArEnergyCut = pset.get<double>("LArEnergyCut");
        fLArVolumeName = pset.get<std::vector<std::string>>("LArVolumeName");
        fLArMaterial = pset.get<std::string>("LArMaterial");

        std::cout << "AuxDetAction: Name of the Volumes to track for the ECAL" << std::endl;
        for(unsigned int i = 0; i < fECALVolumeName.size(); i++) std::cout << fECALVolumeName.at(i) << " ";
        std::cout << std::endl;

        std::cout << "AuxDetAction: Name of the Volumes to track for the LArTPC" << std::endl;
        for(unsigned int i = 0; i < fLArVolumeName.size(); i++) std::cout << fLArVolumeName.at(i) << " ";
        std::cout << std::endl;

        return;
      }

      //-------------------------------------------------------------
      void AuxDetAction::BeginOfEventAction(const G4Event*)
      {
        // Clear any previous information.
        fECALDeposits.clear();
        fLArDeposits.clear();
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

        this->LArSteppingAction(step);

        this->ECALSteppingAction(step);

        }// end of AuxDetAction::SteppingAction

        //------------------------------------------------------------------------------
        void AuxDetAction::EndOfEventAction(const G4Event*)
        {
          //sort per time the hits
          std::sort(fECALDeposits.begin(), fECALDeposits.end());
          std::sort(fLArDeposits.begin(), fLArDeposits.end());
        }

        //------------------------------------------------------------------------------
        void AuxDetAction::LArSteppingAction(const G4Step* step)
        {
          LOG_DEBUG("AuxDetAction")
          << "AuxDetAction::LArSteppingAction";

          // Get the pointer to the track
          G4Track *track = step->GetTrack();

          const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
          const CLHEP::Hep3Vector &stop  = track->GetPosition();

          // If it's a null step, don't use it.
          if(start == stop) return;
          if(step->GetTotalEnergyDeposit() == 0) return;

          // check that we are in the correct material to record a hit
          std::string VolumeName   = this->GetVolumeName(track);

          if( std::find( fLArVolumeName.begin(), fLArVolumeName.end(), VolumeName ) == fLArVolumeName.end() )
          return;

          // check the material
          auto pos = 0.5 * (start + stop);
          TGeoNode *node = fGeoManager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm...

          if(!node){
            LOG_DEBUG("AuxDetAction::LArSteppingAction")
            << "Node not found in "
            << pos.x() << " mm "
            << pos.y() << " mm "
            << pos.z() << " mm";
            return;
          }

          std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
          if ( ! std::regex_match(volmaterial, std::regex(fLArMaterial)) ) return;

          // only worry about energy depositions larger than the minimum required
          if(step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fLArEnergyCut){

            LOG_DEBUG("AuxDetAction::LArSteppingAction")
            << "In volume "
            << VolumeName
            << " Material is "
            << volmaterial
            << " step size is "
            << step->GetStepLength() / CLHEP::cm
            << " cm and deposited "
            << step->GetTotalEnergyDeposit()
            << " MeV of energy with a minimum of "
            << fLArEnergyCut
            << " required.";

            // the step mid point is used for the position of the deposit
            auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
            step->GetPostStepPoint()->GetPosition() );

            // get the track id for this step
            auto trackID  = ParticleListAction::GetCurrentTrackID();
            float time = step->GetPreStepPoint()->GetGlobalTime();

            fLArDeposits.emplace_back(trackID,
              time,//get the time of the first subhit
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
          LOG_DEBUG("AuxDetAction")
          << "AuxDetAction::ECALSteppingAction";

          // Get the pointer to the track
          G4Track *track = step->GetTrack();

          const CLHEP::Hep3Vector &start = step->GetPreStepPoint()->GetPosition();
          const CLHEP::Hep3Vector &stop  = track->GetPosition();

          // If it's a null step, don't use it.
          if(start == stop) return;
          if(step->GetTotalEnergyDeposit() == 0) return;

          // check that we are in the correct material to record a hit
          std::string VolumeName   = this->GetVolumeName(track);

          if( std::find( fECALVolumeName.begin(), fECALVolumeName.end(), VolumeName ) == fECALVolumeName.end() )
          return;

          // check the material
          auto pos = 0.5 * (start + stop);
          TGeoNode *node = fGeoManager->FindNode(pos.x()/CLHEP::cm, pos.y()/CLHEP::cm, pos.z()/CLHEP::cm);//Node in cm...

          if(!node){
            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "Node not found in "
            << pos.x() << " mm "
            << pos.y() << " mm "
            << pos.z() << " mm";
            return;
          }

          std::string volmaterial = node->GetMedium()->GetMaterial()->GetName();
          if ( ! std::regex_match(volmaterial, std::regex(fECALMaterial)) ) return;

          // only worry about energy depositions larger than the minimum required
          if(step->GetTotalEnergyDeposit() * CLHEP::MeV / CLHEP::GeV > fECALEnergyCut){

            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "In volume "
            << VolumeName
            << " Material is "
            << volmaterial
            << " step size is "
            << step->GetStepLength() / CLHEP::cm
            << " cm and deposited "
            << step->GetTotalEnergyDeposit()
            << " MeV of energy with a minimum of "
            << fECALEnergyCut
            << " required.";

            // the step mid point is used for the position of the deposit
            auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
            step->GetPostStepPoint()->GetPosition() );

            unsigned int CaloID = 0;
            if( VolumeName == "IBStrip_L2_vol" )
            CaloID = 1;//Inner Barrel
            if( VolumeName == "OBStrip_L2_vol" )
            CaloID = 2;//Outer Barrel
            if( VolumeName == "IECLayer_L2_vol" )
            CaloID = 3;//Endcap

            // get the track id for this step
            auto trackID  = ParticleListAction::GetCurrentTrackID();
            float time = step->GetPreStepPoint()->GetGlobalTime();

            float edep = this->GetBirksAttenuatedEnergy(step);

            LOG_DEBUG("AuxDetAction::ECALSteppingAction")
            << "Energy deposited "
            << step->GetTotalEnergyDeposit()
            << " MeV after Birks "
            << edep
            << " MeV";

            fECALDeposits.emplace_back(trackID,
              time,//get the time of the first subhit
              edep * CLHEP::MeV / CLHEP::GeV,
              midPoint.x() / CLHEP::cm,
              midPoint.y() / CLHEP::cm,
              midPoint.z() / CLHEP::cm,
              CaloID, 0, 0);

            } // end if enough energy to worry about this step
        }

        //------------------------------------------------------------------------------
        std::string AuxDetAction::GetVolumeName(const G4Track *track)
        {
          std::string VolName = track->GetVolume()->GetName();
          VolName.erase(VolName.length()-3, 3);
          return VolName;
        }

        //------------------------------------------------------------------------------
        float AuxDetAction::GetBirksAttenuatedEnergy(const G4Step* step)
        {
          G4Track* track = step->GetTrack();
          float edep = fEmSaturation->VisibleEnergyDeposition(track->GetParticleDefinition(),
                                                 track->GetMaterialCutsCouple(),
                                                 step->GetStepLength(),
                                                 step->GetTotalEnergyDeposit(),
                                                 step->GetNonIonizingEnergyDeposit());

          return edep;
        }

      } // garg4

    } // gar
