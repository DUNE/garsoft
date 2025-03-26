////////////////////////////////////////////////////////////////////////
/// \file  GArG4.cxx
/// \brief Use Geant4 to run the GArSoft detector simulation
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GARG4GArG4ANA_H
#define GARG4GArG4ANA_H

/// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "nug4/ParticleNavigation/ParticleList.h"

// GArSoft Includes
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"
#include "MCCheater/BackTracker.h"
#include "SimulationDataProducts/EnergyDeposit.h"

// ROOT includes
#include <TTree.h>

// C++ Includes
#include <cstring>
#include <iostream>
#include <sys/stat.h>

typedef struct {
  int run;
  int subrun;
  int event;
} EventInfo;

typedef struct {
  int trackID;
  int pdg;
  float x;
  float y;
  float z;
  float e;
  float dX;
  float t;
} EnergyDep;

typedef struct {
  int trackID;
  int pdg;
  float length;
  float energy;
} ParticleInfo;

namespace simb {
  class MCTruth;
}

namespace sim {
  class ParticleList;
}

///Geant4 interface
namespace gar {
  namespace garg4 {

    class GArG4Ana : public ::art::EDAnalyzer {
    public:
      /// Standard constructor and destructor for an FMWK module.
      explicit GArG4Ana(fhicl::ParameterSet const& pset);
      virtual ~GArG4Ana();

      void analyze(const ::art::Event& evt);
      void beginJob();
      void reconfigure(fhicl::ParameterSet const& pset);

    private:
      std::string fG4ModuleLabel;    ///< module label for the Geant
      std::string fTruthModuleLabel; ///< module label for the Geant

      TTree* fEDepTree;       ///< Tree to keep track of sanity check info
      TTree* fParticleTree;   ///< Tree to keep track of sanity check info
      EventInfo fEvt;         ///< Struct containing event identification
      EnergyDep fEDep;        ///< Struct containing energy deposition info
      ParticleInfo fPartInfo; ///< Struct containing particle info
    };

  } // namespace garg4

  namespace garg4 {

    //-----------------------------------------------------------------------
    // Constructor
    GArG4Ana::GArG4Ana(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
    {
      this->reconfigure(pset);
    }

    //-----------------------------------------------------------------------
    // Destructor
    GArG4Ana::~GArG4Ana() {}

    //-----------------------------------------------------------------------
    void GArG4Ana::beginJob()
    {
      ::art::ServiceHandle<::art::TFileService> tfs;

      fEDepTree = tfs->make<TTree>("EDepTree", "EDepTree");
      fParticleTree = tfs->make<TTree>("ParticleTree", "ParicleTree");

      std::string description("run/I:subrun/I:event/I");

      fEDepTree->Branch("info", &fEvt, description.c_str());
      fParticleTree->Branch("info", &fEvt, description.c_str());

      description = "trackID/I:pdg/I:x/F:y/F:z/F:e/F:dX/F:t/F";

      fEDepTree->Branch("edep", &fEDep, description.c_str());

      description = "trackID/I:pdg/I:length/F:energy/F";

      fParticleTree->Branch("part", &fPartInfo, description.c_str());
    }

    //-----------------------------------------------------------------------
    void GArG4Ana::reconfigure(fhicl::ParameterSet const& p)
    {
      fG4ModuleLabel = p.get<std::string>("GeantModuleLabel");
      fTruthModuleLabel = p.get<std::string>("TruthModuleLabel");

      return;
    }

    //-----------------------------------------------------------------------
    void GArG4Ana::analyze(const ::art::Event& evt)
    {
      // set the event id info
      fEvt.run = evt.run();
      fEvt.subrun = evt.subRun();
      fEvt.event = evt.id().event();

      // get the list of particles from this event
      art::ServiceHandle<cheat::BackTracker> bt;

      // get the energy depositions for this event
      auto edepsCol = evt.getValidHandle<std::vector<gar::sdp::EnergyDeposit>>(fG4ModuleLabel);

      if (edepsCol.failedToGet()) {
        MF_LOG_VERBATIM("GArG4Ana") << "failed to get energy depositions, bail";
        return;
      }

      sim::ParticleList* pList = bt->GetParticleList();

      sim::ParticleList::iterator itr = pList->begin();
      for (; itr != pList->end(); ++itr) {
        auto part = itr->second;
        if (!part) continue;

        fPartInfo.trackID = part->TrackId();
        fPartInfo.pdg = part->PdgCode();
        fPartInfo.energy = part->E();
        fPartInfo.length = part->Trajectory().TotalLength();

        MF_LOG_VERBATIM("GArG4Ana")
          << part->TrackId() << " " << part->PdgCode() << " " << part->E() << " "
          << part->NumberTrajectoryPoints() << " " << part->Trajectory().TotalLength();

        fParticleTree->Fill();
      }

      // loop over the energy deposition collections to fill the tree
      for (auto edep : *edepsCol) {

        // get the MCParticle for this track ID
        auto part = bt->TrackIDToParticle(edep.TrackID());

        fEDep.trackID = edep.TrackID();
        fEDep.pdg = part->PdgCode();
        fEDep.x = edep.X();
        fEDep.y = edep.Y();
        fEDep.z = edep.Z();
        fEDep.dX = edep.dX();
        fEDep.t = edep.Time();
        fEDep.e = edep.Energy();

        MF_LOG_DEBUG("GArG4Ana") << "pos: (" << fEDep.x << ", " << fEDep.y << ", " << fEDep.z
                                 << ") e: " << fEDep.e << " t: " << fEDep.t << " dX: " << fEDep.dX;

        // make the tree flat in terms of the energy depositions
        fEDepTree->Fill();

      } // end loop over collection of EnergyDeposits

      return;
    } // analyze

  } // namespace garg4

  namespace garg4 {

    DEFINE_ART_MODULE(GArG4Ana)

  } // namespace garg4

} // gar

#endif // GARG4GARG4H
