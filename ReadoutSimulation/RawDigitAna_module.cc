//
//  RawDigitAna_module.cc
//
//  Created by Brian Rebel on 3/2/17.
//

#ifndef GAR_ROSIM_RAWDIGITANA
#define GAR_ROSIM_RAWDIGITANA

/// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib_except/exception.h"

// GArSoft Includes
#include "MCCheater/BackTracker.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "Geometry/Geometry.h"
#include "RawDataProducts/RawDigit.h"
#include "CoreUtils/ServiceUtil.h"

// ROOT includes
#include <TTree.h>

// C++ Includes
#include <iostream>
#include <cstring>
#include <sys/stat.h>

typedef struct{
  int run;
  int subrun;
  int event;
} EventInfo;

typedef struct{
  int   trackID;
  int   pdg;
  float x;
  float y;
  float z;
  float e;
  float dX;
  float t;
} EnergyDep;

typedef struct{
  int   channel;
  float x;
  float y;
  float z;
} ChannelInfo;

namespace simb{
  class MCTruth;
}

namespace sim{
  class ParticleList;
}

namespace gar {
  namespace rosim {
    
    class RawDigitAna : public ::art::EDAnalyzer{
    public:
      
        /// Standard constructor and destructor for an FMWK module.
      explicit RawDigitAna(fhicl::ParameterSet const& pset);
      virtual ~RawDigitAna();
      
      void analyze (const ::art::Event& evt);
      void beginJob();
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      std::string  fReadoutModuleLabel; ///< module label for the Geant
      TTree*       fChannelTree;        ///< Tree to keep track of sanity check info
      EventInfo    fEvt;                ///< Struct containing event identification
      EnergyDep    fEDep;               ///< Struct containing energy deposition info
      ChannelInfo  fChannelInfo;        ///< Struct containing channel info
    };
    
  } // namespace rosim
  
  namespace rosim {
    
    //-----------------------------------------------------------------------
    // Constructor
    RawDigitAna::RawDigitAna(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
    {
      this->reconfigure(pset);
    }
    
    //-----------------------------------------------------------------------
    // Destructor
    RawDigitAna::~RawDigitAna()
    {
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::beginJob()
    {
      ::art::ServiceHandle<::art::TFileService> tfs;
      
      fChannelTree = tfs->make<TTree>("ChannelTree", "ChannelTree");
      
      std::string description("run/I:subrun/I:event/I");
      
      fChannelTree->Branch("info", &fEvt, description.c_str());
      
      description = "trackID/I:pdg/I:x/F:y/F:z/F:e/F:dX/F:t/F";
      
      fChannelTree->Branch("edep", &fEDep, description.c_str());
      
      description = "channel/I:x/F:y/F:z/F";
      
      fChannelTree->Branch("chan", &fChannelInfo, description.c_str());
      
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::reconfigure(fhicl::ParameterSet const& p)
    {
      fReadoutModuleLabel = p.get< std::string >("ReadoutModuleLabel");
      
      return;
    }
    
    //-----------------------------------------------------------------------
    void RawDigitAna::analyze(const ::art::Event& evt)
    {
      // set the event id info
      fEvt.run    = evt.run();
      fEvt.subrun = evt.subRun();
      fEvt.event  = evt.id().event();
      
      // get the list of particles from this event
      auto bt = gar::providerFrom<cheat::BackTracker>();
      
      // get the raw digits for this event
      auto digCol = evt.getValidHandle<std::vector<gar::raw::RawDigit> >(fReadoutModuleLabel);
      
      if( digCol.failedToGet() ){
        LOG_VERBATIM("RawDigitAna")
        << "failed to get raw digits, bail";
        return;
      }
 
      // get the find many for the energy depositions
      ::art::FindMany<gar::sdp::EnergyDeposit> fmEnergyDep(digCol, evt, fReadoutModuleLabel);
      
      if( !fmEnergyDep.isValid() ){
        LOG_WARNING("RawDigitAna")
        << "Unable to find valid association between RawDigits and "
        << "energy deposits, no analysis in this module is possible";
        
        return;
      }

      auto geo = gar::providerFrom<geo::Geometry>();

      float        xyz[3]     = {0.};
      float        xyzChan[3] = {0.};
      unsigned int chan       = 0;
      
      for(size_t d = 0; d < digCol->size(); ++d){
        std::vector<const gar::sdp::EnergyDeposit*> edepsCol;
        fmEnergyDep.get(d, edepsCol);
        if(edepsCol.size() < 1) continue;
        
	LOG_DEBUG("RawDigitAna")
        << "There are "
        << edepsCol.size()
        << " energy depositions for channel "
        << (*digCol)[d].Channel();
        
        geo->ChannelToPosition((*digCol)[d].Channel(), xyz);
        
        // fill the channel information
        fChannelInfo.channel = (int)(*digCol)[d].Channel();
        fChannelInfo.x       = xyz[0];
        fChannelInfo.y       = xyz[1];
        fChannelInfo.z       = xyz[2];
        
        // loop over the energy deposition collections to fill the tree
        for(auto edep : edepsCol){
        
          // get the MCParticle for this track ID
          auto part = bt->TrackIDToParticle(edep->TrackID());
          
          fEDep.trackID = edep->TrackID();
          fEDep.pdg     = part->PdgCode();
          fEDep.x       = edep->X();
          fEDep.y       = edep->Y();
          fEDep.z       = edep->Z();
          fEDep.dX      = edep->dX();
          fEDep.t       = edep->Time();
          fEDep.e       = edep->Energy();
          
	  LOG_DEBUG("RawDigitAna")
          << "pos: ("
          << fEDep.x
          << ", "
          << fEDep.y
          << ", "
          << fEDep.z
          << ") e: "
          << fEDep.e
          << " t: "
          << fEDep.t
          << " dX: "
          << fEDep.dX;
          
          // closure test for the channel map
          xyz[0] = fEDep.x;
          xyz[1] = fEDep.y;
          xyz[2] = fEDep.z;
          
          chan = geo->NearestChannel(xyz);
	  if (chan != geo->GapChannelNumber())
	    {
               geo->ChannelToPosition(chan, xyzChan);
          
//          if(std::abs(xyz[1] - fChannelInfo.y) > 0.33)
//            LOG_VERBATIM("RawDigitAna")
//            << "Channel "
//            << fChannelInfo.channel
//            << " / "
//            << chan
//            << " is off from the energy deposit in y: ("
//            << fChannelInfo.y
//            << ", "
//            << fChannelInfo.z
//            << ") vs ("
//            << xyz[1]
//            << ", "
//            << xyz[2]
//            << ") "
//            << std::abs(xyz[1] - fChannelInfo.y);
//
//          if(std::abs(xyz[2] - fChannelInfo.z) > 10.33)
//            LOG_VERBATIM("RawDigitAna")
//            << "Channel "
//            << fChannelInfo.channel
//            << " / "
//            << chan
//            << " is off from the energy deposit in z: ("
//            << fChannelInfo.y
//            << ", "
//            << fChannelInfo.z
//            << ") vs ("
//            << xyz[1]
//            << ", "
//            << xyz[2]
//            << ") "
//            << std::abs(xyz[2] - fChannelInfo.z);

	    }

          // make the tree flat in terms of the energy depositions
          fChannelTree->Fill();
          
        } // end loop over collection of EnergyDeposits
      } // end loop over RawDigit collection
      
      return;
    } // analyze
    
  } // namespace rosim
  
  namespace rosim {
    
    DEFINE_ART_MODULE(RawDigitAna)
    
  } // namespace rosim
  
} // gar

#endif // GAR_ROSIM_RAWDIGITANA

