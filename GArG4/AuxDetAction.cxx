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
        fContributions.clear();
        fEnergy.clear();
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
        for(unsigned int i = 0; i < fVolumeName.size(); i++)
        {
          if(volname.find(fVolumeName.at(i)) != std::string::npos){
            VolumeName = fVolumeName.at(i);
            break;
          }
          else
          return;
        }

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

        LOG_DEBUG("AuxDetAction")
        << "In volume "
        << VolumeName
        << " step size is "
        << step->GetStepLength()
        << " and deposited "
        << step->GetTotalEnergyDeposit()
        << " MeV of energy with a minimum of "
        << fEnergyCut
        << " required.";

        // only worry about energy depositions larger than the minimum required
        if(step->GetTotalEnergyDeposit() > fEnergyCut){

          // the step mid point is used for the position of the deposit
          auto midPoint = 0.5 * (step->GetPreStepPoint()->GetPosition() +
          step->GetPostStepPoint()->GetPosition() );

          int CaloID = -1;
          if( (strncmp(volname.c_str(), "IBTile_L2_vol_PV", 16) == 0) )
          CaloID = 1;//Inner Barrel
          if( (strncmp(volname.c_str(), "OBTile_L2_vol_PV", 16) == 0) )
          CaloID = 2;//Outer Barrel
          if( (strncmp(volname.c_str(), "IECTile_L2_vol_PV", 17) == 0) )
          CaloID = 3;//Endcap

          //Get the geometry information from the node
          TGeoBBox *tile = (TGeoBBox*)node->GetMotherVolume()->GetShape();
          const double origin[3] = {0, 0, 0.45};//in cm
          double master[3];
          //Get the position of the node in the world
          geomanager->LocalToMaster(origin, &master[0]);//in cm

          double x = master[0]*CLHEP::cm;
          double y = master[1]*CLHEP::cm;
          double z = master[2]*CLHEP::cm;

          LOG_DEBUG("AuxDetAction")
          << node->GetMotherVolume()->GetName() << " "
          << tile->GetDX()*CLHEP::cm*2 << " mm "
          << tile->GetDY()*CLHEP::cm*2 << " mm "
          << tile->GetDZ()*CLHEP::cm*2 << " mm, local "
          << origin[0]*CLHEP::cm << " "
          << origin[1]*CLHEP::cm << " "
          << origin[2]*CLHEP::cm << " master "
          << x << " mm "
          << y << " mm "
          << z << " mm";

          std::vector<gar::sdp::CaloDeposit> v;
          gar::sdp::CaloDeposit Contribhit(ParticleListAction::GetCurrentTrackID(), step->GetPreStepPoint()->GetGlobalTime(), step->GetTotalEnergyDeposit(), midPoint.x(), midPoint.y(), midPoint.z(), CaloID, v);
          fContributions[CaloID][x][y][z].push_back(Contribhit);
          fEnergy[CaloID][x][y][z] += step->GetTotalEnergyDeposit();

        } // end if enough energy to worry about this step

      }// end of AuxDetAction::SteppingAction

      //------------------------------------------------------------------------------
      void AuxDetAction::EndOfEventAction(const G4Event*)
      {
        //Check the map
        for(auto &id : fContributions)//std::map<int, std::map< double, std::map<double, std::map<double, std::vector<gar::sdp::CaloDeposit> > > > >
        {
          for(auto &x : id.second)//std::map< double, std::map<double, std::map<double, std::vector<gar::sdp::CaloDeposit> > > >
          {
            for(auto &y : x.second)//std::map<double, std::map<double, std::vector<gar::sdp::CaloDeposit> > >
            {
              for(auto &z : y.second)//std::map<double, std::vector<gar::sdp::CaloDeposit> >
              {
                //sort per time the subhits
                std::sort(z.second.begin(), z.second.end());

                //CaloID << x << y << z << vector of contributions
                fDeposits.emplace_back(-1,
                  z.second.at(0).Time(),//get the time of the first subhit
                  fEnergy[id.first][x.first][y.first][z.first],
                  x.first,
                  y.first,
                  z.first,
                  id.first,
                  z.second);
                }
              }
            }
          }

          // sort per time the hits
          std::sort(fDeposits.begin(), fDeposits.end());
        }


      } // garg4

    } // gar
