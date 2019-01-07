////////////////////////////////////////////////////////////////////////
//
//
// \file: BackTracker_service.cc
//
// brebel@fnal.gov
//
////////////////////////////////////////////////////////////////////////

// GArSoft includes
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nutools/ParticleNavigation/EmEveIdCalculator.h"

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"

namespace gar{
  namespace cheat{
    
    //----------------------------------------------------------------------
    BackTracker::BackTracker(fhicl::ParameterSet     const& pset,
                             ::art::ActivityRegistry      & reg)
    : BackTrackerCore(pset)
    {
      reg.sPreProcessEvent.watch(this, &BackTracker::Rebuild);
    }
    
    //----------------------------------------------------------------------
    BackTracker::~BackTracker()
    {
    }

    //----------------------------------------------------------------------
    void BackTracker::beginJob()
    {

      return;
    }

    //----------------------------------------------------------------------
    void BackTracker::Rebuild(::art::Event const& evt, art::ScheduleContext)
    {
      RebuildNoSC(evt);
    }
    
    void BackTracker::RebuildNoSC(::art::Event const& evt)
    {
      // do nothing if this is data
      if(evt.isRealData()) return;
      
      // we have a new event, so clear the channel to EnergyDeposit collection
      for(auto vec : fChannelToEDepCol) vec.clear();
      fChannelToEDepCol.clear();
      
      // get the particles from the event
      auto pHandle = evt.getValidHandle<std::vector<simb::MCParticle> >(fG4ModuleLabel);
      
      // first check to see if we got called to early, that is the particles are
      // and sim channels are not made yet in a MC production job
      // if that is the case, we'll take care of it later
      if(pHandle.failedToGet()){
        LOG_WARNING("BackTracker")
        << "failed to get handle to simb::MCParticle from "
        << fG4ModuleLabel
        << ", return";
        return;
      }
      
      // Clear out anything remaining from previous calls to Rebuild
      fParticleList.clear();
      fMCTruthList .clear();
      
      ::art::FindOneP<simb::MCTruth> fo(pHandle, evt, fG4ModuleLabel);
      
      if( fo.isValid() ){
        for(size_t p = 0; p < pHandle->size(); ++p){
          
          simb::MCParticle *part = new simb::MCParticle(pHandle->at(p));
          fParticleList.Add(part);
          
          // get the simb::MCTruth associated to this sim::ParticleList
          try{
            ::art::Ptr<simb::MCTruth> mct = fo.at(p);
            if(fMCTruthList.size() < 1) fMCTruthList.push_back(mct);
            else{
              // check that we are not adding a simb::MCTruth twice to the collection
              // we know that all the particles for a given simb::MCTruth are put into the
              // collection of particles at the same time, so we can just check that the
              // current ::art::Ptr has a different id than the last one put
              if(!(mct == fMCTruthList.back())) fMCTruthList.push_back(mct);
            }
            
            // fill the track id to mctruth index map
            fTrackIDToMCTruthIndex[pHandle->at(p).TrackId()] = fMCTruthList.size() - 1;
          }
          catch(cet::exception &ex){
            LOG_WARNING("BackTracker")
            << "unable to find MCTruth from ParticleList "
            << "created in "
            << fG4ModuleLabel
            << " any attempt to get the MCTruth objects from "
            << "the backtracker will fail\n"
            << "message from caught exception:\n"
            << ex;
          }
        }// end loop over particles to get MCTruthList
      }// end if fo.isValid()
      
      fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);
      
      LOG_DEBUG("BackTracker")
      << "BackTracker has "
      << GetSetOfTrackIDs().size()
      << " tracks.  The particles are:\n"
      << fParticleList;
      
      // look for the RawDigit collection from the event and create a FindMany mapping
      // between the digits and energy deposits if it exists.  Really just map between the
      // Channel numbers and the energy deposits
      art::Handle<std::vector<gar::raw::RawDigit> > digCol;
      evt.getByLabel(fRawDataLabel, digCol);
      
      if( !digCol.isValid() ){
        //LOG_WARNING("BackTracker")
        //<< "Unable to find valid collection of RawDigits in "
        //<< fRawDataLabel
        //<< " no backtracking of hits will be possible";
        
        return;
      }
      
      ::art::FindMany<gar::sdp::EnergyDeposit> fmEnergyDep(digCol, evt, fRawDataLabel);
      
      if( !fmEnergyDep.isValid() ){
        LOG_WARNING("BackTracker")
        << "Unable to find valid association between RawDigits and "
        << "energy deposits, no backtracking of hits will be possible";
        
        return;
      }
      
      // figure out how large the channel to energy deposit collection should be
      auto geo = gar::providerFrom<geo::Geometry>();
      fChannelToEDepCol.resize(geo->NChannels());
      
      LOG_DEBUG("BackTrackerServ")
      << "There are "
      << geo->NChannels()
      << " "
      << fChannelToEDepCol.size()
      << " channels in the geometry";
      
      std::vector<const gar::sdp::EnergyDeposit*> eDeps;
      
      for(size_t d = 0; d < digCol->size(); ++d){
        fmEnergyDep.get(d, eDeps);
        if(eDeps.size() < 1) continue;
        
        fChannelToEDepCol[ (*digCol)[d].Channel() ].swap(eDeps);
      }
      
      //      for(auto mc : fMCTruthList)
      //        LOG_DEBUG("BackTracker") << *(mc);
      
      
      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      auto garProp = gar::providerFrom<detinfo::GArPropertiesService>();
      
      fInverseVelocity = 1. / detProp->DriftVelocity(detProp->Efield(),
                                                     detProp->Temperature());
      fLongDiffConst   = std::sqrt(2. * garProp->LongitudinalDiffusion());
      
      
      return;
    }
    
    DEFINE_ART_SERVICE(BackTracker)

  } // namespace cheat
} // namespace gar
