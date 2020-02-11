//
//  BackTrackerCore.cxx
//  Created by Brian Rebel on 3/13/17.
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nutools/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "SimulationDataProducts/sim.h"
#include "CoreUtils/ServiceUtil.h"
#include "MCCheater/BackTrackerCore.h"
#include "Geometry/Geometry.h"

#include "TMath.h"

namespace gar{
  namespace cheat{
    
    //--------------------------------------------------------------------------
    BackTrackerCore::BackTrackerCore(fhicl::ParameterSet const& pset)
      : fClocks(nullptr), fGeo(nullptr)
    {
      fClocks = gar::providerFrom<detinfo::DetectorClocksService>();
      fGeo     = gar::providerFrom<geo::Geometry>();
      this->reconfigure(pset);
      
      return;
    }

    //--------------------------------------------------------------------------
    BackTrackerCore::~BackTrackerCore()
    {
      
    }
    
    //----------------------------------------------------------------------
    void BackTrackerCore::reconfigure(const fhicl::ParameterSet& pset)
    {
      fG4ModuleLabel        = pset.get<std::string>("G4ModuleLabel",           "geant");
      fRawDataLabel         = pset.get<std::string>("RawDataLabel",            "daq");
      fMinHitEnergyFraction = pset.get<double     >("MinimumHitEnergyFraction", 0.1);
      fDisableRebuild       = pset.get<bool       >("DisableRebuild",           false);
    }
    
    //----------------------------------------------------------------------
    const simb::MCParticle* BackTrackerCore::TrackIDToParticle(int const& id) const
    {
      auto part_it = fParticleList.find(id);
      
      if(part_it == fParticleList.end()){
        LOG_WARNING("BackTracker")
        << "can't find particle with track id "
        << id
        << " in sim::ParticleList returning null pointer";
        
        return nullptr;
      }
      
      return part_it->second;
    }
    
    //----------------------------------------------------------------------
    const simb::MCParticle* BackTrackerCore::TrackIDToMotherParticle(int const& id) const
    {
      // get the mother id from the particle navigator
      // the EveId was adopted in the Rebuild method
      
      return this->TrackIDToParticle(fParticleList.EveId(std::abs(id)));
    }
    
    //----------------------------------------------------------------------
    ::art::Ptr<simb::MCTruth> const& BackTrackerCore::TrackIDToMCTruth(int const& id) const
    {
      // find the entry in the MCTruth collection for this track id
      size_t mct = fTrackIDToMCTruthIndex.find(std::abs(id))->second;
      
      if(/* mct < 0 || */ mct > fMCTruthList.size() )
        throw cet::exception("BackTracker")
        << "attempting to find MCTruth index for "
        << "out of range value: "
        << mct
        << "/" << fMCTruthList.size();
      
      return fMCTruthList[mct];
    }
    
    //----------------------------------------------------------------------
    ::art::Ptr<simb::MCTruth> const& BackTrackerCore::ParticleToMCTruth(const simb::MCParticle* p) const
    {
      return this->TrackIDToMCTruth(p->TrackId());
    }
    
    //----------------------------------------------------------------------
    std::vector<const simb::MCParticle*> BackTrackerCore::MCTruthToParticles(::art::Ptr<simb::MCTruth> const& mct) const
    {
      std::vector<const simb::MCParticle*> ret;
      
      // sim::ParticleList::value_type is a pair (track ID, particle pointer)
      for(auto const& TrackIDpair: fParticleList) {
        if( TrackIDToMCTruth(TrackIDpair.first) == mct )
          ret.push_back(TrackIDpair.second);
      }
      
      return ret;
    }
    
    //----------------------------------------------------------------------
    std::vector<HitIDE> BackTrackerCore::HitToTrackID(::art::Ptr<gar::rec::Hit> const& hit) const
    {
      std::vector<HitIDE> ides;
      
      this->ChannelToTrackID(ides,
                             hit->Channel(),
                             hit->StartTime(),
                             hit->EndTime());
      
      return ides;
    }

    //----------------------------------------------------------------------
    std::vector<HitIDE> BackTrackerCore::HitToTrackID(gar::rec::Hit const& hit) const
    {
      std::vector<HitIDE> ides;
      
      this->ChannelToTrackID(ides,
                             hit.Channel(),
                             hit.StartTime(),
                             hit.EndTime());
      
      return ides;
    }
    
    //----------------------------------------------------------------------
    std::vector<std::vector<::art::Ptr<gar::rec::Hit>>> const BackTrackerCore::TrackIDsToHits(std::vector<::art::Ptr<gar::rec::Hit>> const& allhits,
                                                                                              std::vector<int>                       const& tkIDs) const
    {
      // returns a subset of the hits in the allhits collection that are matched
      // to MC particles listed in tkIDs
      
      // temporary vector of TrackIDs and Ptrs to hits so only one
      // loop through the (possibly large) allhits collection is needed
      std::vector< std::pair<int, ::art::Ptr<gar::rec::Hit> > > hitList;
      std::vector<HitIDE> hids;
      for(auto hit : allhits){
        hids.clear();
        
        this->ChannelToTrackID(hids,
                               hit->Channel(),
                               hit->StartTime(),
                               hit->EndTime());
        
        for(auto const& hid : hids) {
          for(auto const& itkid : tkIDs) {
            if(hid.trackID == itkid) {
              if(hid.energyFrac > fMinHitEnergyFraction)
                hitList.push_back(std::make_pair(itkid, hit));
            }
          } // itkid
        } // itid
      } // itr
      
      // now build the truHits vector that will be returned to the caller
      std::vector< std::vector< ::art::Ptr<gar::rec::Hit> > > truHits;
      
      // temporary vector containing hits assigned to one MC particle
      std::vector<::art::Ptr<gar::rec::Hit>> tmpHits;
      for(auto const& itkid : tkIDs) {
        tmpHits.clear();
        
        for(auto itr : hitList) {
          if(itkid == itr.first) tmpHits.push_back(itr.second);
        }
        truHits.push_back(tmpHits);
      }
      
      return truHits;
    }
    
    //----------------------------------------------------------------------
    // plist is assumed to have adopted the appropriate EveIdCalculator prior to
    // having been passed to this method. It is likely that the EmEveIdCalculator is
    // the one you always want to use
    std::vector<HitIDE> BackTrackerCore::HitToEveID(::art::Ptr<gar::rec::Hit> const& hit) const
    {
      std::vector<HitIDE> ides = this->HitToTrackID(hit);
      
      // make a map of evd ID values and fraction of energy represented by
      // that eve id in this hit
      std::map<int, HitIDE> eveToE;
      
      for(auto const& hid : ides){
        eveToE[fParticleList.EveId( hid.trackID )] = hid;
      }
      
      // now fill the eveides vector from the map
      std::vector<HitIDE> eveides;
      eveides.reserve(eveToE.size());
      for(auto itr : eveToE){
        eveides.emplace_back(itr.first,
                             itr.second.energyFrac,
                             itr.second.numElectrons);
      }
      
      return eveides;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTrackerCore::GetSetOfEveIDs() const
    {
      std::set<int> eveIDs;
      
      sim::ParticleList::const_iterator plitr = fParticleList.begin();
      while(plitr != fParticleList.end() ){
        int eveID = fParticleList.EveId((*plitr).first);
        // look to see if this eveID is already in the set
        if( eveIDs.find(eveID) == eveIDs.end() ) eveIDs.insert(eveID);
        plitr++;
      }
      
      return eveIDs;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTrackerCore::GetSetOfTrackIDs() const
    {
      // fParticleList::value_type is a pair (track, particle pointer)
      std::set<int> trackIDs;
      for (auto const& pl : fParticleList)
        trackIDs.insert(pl.first);
      
      return trackIDs;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTrackerCore::GetSetOfEveIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits) const
    {
      std::set<int> eveIDs;
      
      //std::vector< ::art::Ptr<gar::rec::Hit> >::const_iterator itr = hits.begin();
      for(auto itr : hits){
        
        // get the eve ids corresponding to this hit
        const std::vector<HitIDE> ides = HitToEveID(itr);
        
        // loop over the ides and extract the track ids
        for(auto const& ide : ides) eveIDs.insert(ide.trackID);
        
      }
      
      return eveIDs;
    }
    
    //----------------------------------------------------------------------
    std::set<int> BackTrackerCore::GetSetOfTrackIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits) const
    {
      std::set<int>       trackIDs;
      std::vector<HitIDE> ides;
      
      for(auto itr : hits ){
        
        ides.clear();
        
        this->ChannelToTrackID(ides,
                               itr->Channel(),
                               itr->StartTime(),
                               itr->EndTime());
        
        // loop over the ides and extract the track ids
        for(auto const& hid : ides) {
          trackIDs.insert(hid.trackID);
        }
      }
      
      return trackIDs;
    }
    
    //----------------------------------------------------------------------
    double BackTrackerCore::HitCollectionPurity(std::set<int>                                   trackIDs,
                                                std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                                bool                                            weightByCharge) const
    {
      // get the list of EveIDs that correspond to the hits in this collection
      // if the EveID shows up in the input list of trackIDs, then it counts
      float total   = 0;
      float desired = 0.;
      float weight  = 1.;
      
      // don't have to check the view in the hits collection because
      // those are assumed to be from the object we are testing and will
      // the correct view by definition then.
      for(auto hit : hits){
        
        weight = (weightByCharge) ? hit->Signal() : 1.;
        
        std::vector<HitIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        total += weight;
        
        // don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(auto const& tid : hitTrackIDs){
          if(trackIDs.find(tid.trackID) != trackIDs.end()){
            desired += weight;
            break;
          }
        }
        
      }// end loop over hits
      
      double purity = 0.;
      
      if(total > 0) purity = desired/total;
      
      return purity;
    }
    
    //----------------------------------------------------------------------
    double BackTrackerCore::HitCollectionEfficiency(std::set<int>                                   trackIDs,
                                                    std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                                    std::vector< ::art::Ptr<gar::rec::Hit> > const& allhits,
                                                    bool                                            weightByCharge) const
    {
      // get the list of EveIDs that correspond to the hits in this collection
      // and the energy associated with the desired trackID
      float desired = 0.;
      float total   = 0.;
      float weight  = 0.;
      
      // don't have to check the view in the hits collection because
      // those are assumed to be from the object we are testing and will
      // the correct view by definition then.
      for(auto hit : hits){
        
        std::vector<HitIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        weight = (weightByCharge) ? hit->Signal() : 1.;
        
        // don't worry about hits where the energy fraction for the chosen
        // trackID is < 0.1
        // also don't double count if this hit has more than one of the
        // desired track IDs associated with it
        for(auto tid : hitTrackIDs){
          if(trackIDs.find(tid.trackID) != trackIDs.end() &&
             tid.energyFrac             >= fMinHitEnergyFraction){
            desired += weight;
            break;
          }
        }
      }// end loop over hits
      
      // now figure out how many hits in the whole collection are associated with this id
      for(auto hit : allhits){
        
        weight = (weightByCharge) ? hit->Signal() : 1.;
        
        std::vector<HitIDE> hitTrackIDs = this->HitToTrackID(hit);
        
        for(auto const& tid : hitTrackIDs){
          // don't worry about hits where the energy fraction for the chosen
          // trackID is < 0.1
          // also don't double count if this hit has more than one of the
          // desired track IDs associated with it
          if(trackIDs.find(tid.trackID) != trackIDs.end() &&
             tid.energyFrac             >= fMinHitEnergyFraction){
            total += weight;
            break;
          }
        }
        
      }// end loop over all hits
      
      double efficiency = 0.;
      if(total > 0.) efficiency = desired/total;
      
      return efficiency;
    }
    
    //--------------------------------------------------------------------------
    std::map<size_t, const sdp::EnergyDeposit*> BackTrackerCore::ChannelTDCToEnergyDeposit(raw::Channel_t channel) const
    {
      if(channel > fChannelToEDepCol.size())
        throw cet::exception("BackTrackerCore")
        << "Requested channel "
        << channel
        << " for backtracking TDCs is out of range, bail";
      
      LOG_DEBUG("BackTracker")
      << "There are "
      << fChannelToEDepCol[channel].size()
      << " energy depositions for channel "
      << channel;
      
      float chanpos[3]={0.,0.,0};
      fGeo->ChannelToPosition(channel, chanpos);

      std::map<size_t, const sdp::EnergyDeposit*> edeps;
      
      // loop over all the energy deposits for this channel
      size_t centralTDC = 0.;

      for(auto const* edep : fChannelToEDepCol[channel]){
        centralTDC = (size_t)fClocks->TPCG4Time2TDC(edep->Time() + TMath::Abs(edep->X()-chanpos[0]) * fInverseVelocity);
        
        LOG_DEBUG("BackTrackerCore")
        << "centralTDC: "
        << centralTDC
        << " time "
        << edep->Time()
        << " energy "
        << edep->Energy();
        
        edeps[centralTDC] = edep;
      }
      
      return edeps;
    }
    
    //----------------------------------------------------------------------
    void BackTrackerCore::ChannelToTrackID(std::vector<HitIDE>      & hitIDEs,
                                           raw::Channel_t      const& channel,
                                           double              const  start,
                                           double              const  end) const
    {
      
      // loop over the electrons in the channel and grab those that are in time
      // with the identified hit start and stop times
    
      // get the track ids corresponding to this hit
      unsigned int start_tdc = (unsigned int)start;
      unsigned int end_tdc   = (unsigned int)end;
      // clang though these were unnecessary as start_tdc and end_tdc are unsigned
      //if(start_tdc < 0) start_tdc = 0;
      //if(end_tdc   < 0) end_tdc   = 0;
      
      hitIDEs.clear();
      
      double totalE = 0.;
      
      if(fChannelToEDepCol.size() < channel){
        LOG_WARNING("BackTracker")
        << "Attempting to find energy deposits for channel "
        << channel
        << " while there are only "
        << fChannelToEDepCol.size()
        << " channels available, return with empty result.";
        
        return;
      }
      
      auto chanEDeps = fChannelToEDepCol[channel];
      
      if(chanEDeps.size() < 1){
        LOG_WARNING("BackTracker")
        << "No sdp::EnergyDeposits for selected channel: "
        << channel
        << ". There is no way to backtrack the given hit, return"
        << " an empty vector.";
        return;
      }
      
      // first get the total energy represented by all track ids for
      // this channel and range of tdc values
      unsigned short tdc         = 0;
      int            prevTrackID = chanEDeps.front()->TrackID();
      float          energy      = 0.;
      
      for(auto const& edep : chanEDeps){
        
        if(edep->TrackID() == gar::sdp::NoParticleId) continue;
        
        tdc = fClocks->TPCG4Time2TDC( edep->Time() );
        
        // only worry about TDC values in the correct range, stop looking if
        // we go past the upper limit.  The IDEs are sorted by channel, then
        // TDC value, then track ID when fChannelToIDEs is filled.
        if     (tdc <  start_tdc) continue;
        else if(tdc == start_tdc){
          prevTrackID = edep->TrackID();
          totalE      = 0.;
          energy      = 0.;
        }
        if(tdc > end_tdc  ) break;
        
        if(prevTrackID != edep->TrackID()){
          
          // we are on a new TrackID, so store the information for the
          // previous one.  For now set the total energy value to 0
          // and the fractional energy to be just the total electrons
          // for this TrackID
          hitIDEs.emplace_back(prevTrackID, 0., energy);
          
          prevTrackID = edep->TrackID();
          energy      = 0.;
        }
        
        totalE += edep->Energy();
        energy += edep->Energy();
        
      }
      
      // loop over the hitIDEs to set the fractional energy for each TrackID
      // and the total energy as well
      if(totalE < 1.e-5) totalE = 1.;
      
      for(size_t i = 0; i < hitIDEs.size(); ++i)
        hitIDEs[i].energyFrac = hitIDEs[i].numElectrons / totalE;
      
      return;
    }
    
    //----------------------------------------------------------------------
    std::vector<float> BackTrackerCore::HitToXYZ(::art::Ptr<gar::rec::Hit> const& hit) const
    {
      std::vector<float> xyz;
      
      // There is a data member which is a vector of all energy deposits for each channel
      // in the event.  Let's look to see if the vector is non-empty for this one
      // If it is empty, it could be a noise hit
      if(fChannelToEDepCol.size() < hit->Channel()){
        LOG_WARNING("BackTracker")
        << "Attempting to back track a hit from a channel without any corresponding EnergyDeposits, "
        << "return empty vector: "
        << fChannelToEDepCol.size()
        << "/"
        << hit->Channel();
        
        return xyz;
      }
      
      if(fChannelToEDepCol[hit->Channel()].size() < 1){
        LOG_WARNING("BackTracker")
        << "Attempting to back track a hit without any corresponding EnergyDeposits, "
        << "return empty vector";
        return xyz;
      }
      
      // We know time extent of the hit, so that should allow us to narrow down
      // the EnergyDeposits
      // Since multiple times correspond to a single hit, we need to
      // average over all contributing positions
      float avX    = 0.;
      float avY    = 0.;
      float avZ    = 0.;
      float wgtSum = 0.;
      
      unsigned int start = (unsigned int)hit->StartTime();
      unsigned int tdc   = 0;
      unsigned int end   = (unsigned int)hit->EndTime();
      
      for(auto const& edep : fChannelToEDepCol[hit->Channel()]){
        
        // check the TDC value of the deposit to make sure it is in the
        // correct range
        // TODO: strictly speaking, the TDC retrieved this way may not be correct
        // because we are getting a single TDC value for the EnergyDeposit, but
        // the simulation breaks the deposits into clusters of electrons to
        // account for diffusion.  This method is approximately correct
        tdc = fClocks->TPCG4Time2TDC( edep->Time() );
        if(tdc < start || tdc > end) continue;
        
        avX += edep->X() * edep->Energy();
        avY += edep->Y() * edep->Energy();
        avZ += edep->Z() * edep->Energy();
        
        wgtSum += edep->Energy();
      }
      
      if(wgtSum == 0.){
        LOG_WARNING("BackTracker")
        << "no energy recorded for any EnergyDeposit belonging to this hit,"
        << " return an empty vector";
        return xyz;
      }
      
      xyz.push_back(avX / wgtSum);
      xyz.push_back(avY / wgtSum);
      xyz.push_back(avZ / wgtSum);
      
      return xyz;
    }

  } // cheat
} // gar
