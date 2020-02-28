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
      fGeo    = gar::providerFrom<geo::Geometry>();

      fG4ModuleLabel        = pset.get<std::string>("G4ModuleLabel",           "geant");
      fRawDataLabel         = pset.get<std::string>("RawDataLabel",            "daq");
      fMinHitEnergyFraction = pset.get<double     >("MinimumHitEnergyFraction", 0.1);
      fDisableRebuild       = pset.get<bool       >("DisableRebuild",           false);
      
      return;
    }

    //--------------------------------------------------------------------------
    BackTrackerCore::~BackTrackerCore() {}
    



    
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
    art::Ptr<simb::MCTruth> const&
	BackTrackerCore::ParticleToMCTruth(const simb::MCParticle* p) const {

		int id = p->TrackId();
		// find the entry in the MCTruth collection for this track id
      	size_t mct = fTrackIDToMCTruthIndex.find(std::abs(id))->second;
      
		if( mct > fMCTruthList.size() ) {
			throw cet::exception("BackTracker")
				<< "attempting to find MCTruth index for out-of-range value: "
				<< mct << "/" << fMCTruthList.size();
		}
		return fMCTruthList[mct];
	}



    //----------------------------------------------------------------------
    std::vector<HitIDE>
	BackTrackerCore::HitToHitIDEs(art::Ptr<gar::rec::Hit> const& hit) const {

		std::vector<HitIDE> ides;
		ides = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());
		return ides;
	}



    //----------------------------------------------------------------------
    std::vector<HitIDE>
	BackTrackerCore::HitToHitIDEs(gar::rec::Hit const& hit) const {

		std::vector<HitIDE> ides;
		ides = this->ChannelToHitIDEs(hit.Channel(),hit.StartTime(),hit.EndTime());
		return ides;
    }



	//----------------------------------------------------------------------
	std::vector<HitIDE>
	BackTrackerCore::ChannelToHitIDEs(raw::Channel_t const& channel,
                                    double const  start, double const stop) const {


		// loop over the electrons in the channel and grab those that are in time
		// with the identified hit start and stop times

		// get the track ids corresponding to this hit
		unsigned int start_tdc = (unsigned int)start;
		unsigned int stop_tdc   = (unsigned int)stop;

		std::vector<HitIDE> hitIDEs;

		double totalE = 0.;
      
		if(fChannelToEDepCol.size() < channel){
			LOG_WARNING("BackTracker")
				<< "Attempting to find energy deposits for channel " << channel
				<< " while there are only " << fChannelToEDepCol.size()
				<< " channels available, returning an empty vector.";
			return hitIDEs;
		}

		auto chanEDeps = fChannelToEDepCol[channel];

		if(chanEDeps.size() < 1){
			LOG_WARNING("BackTracker")
				<< "No sdp::EnergyDeposits for selected channel: " << channel
				<< ". There is no way to backtrack the given hit, returning"
				<< " an empty vector.";
			return hitIDEs;
		}

		std::map<int,size_t> tidmap;

		// first get the total energy represented by all track ids for
		// this channel and range of tdc values
		unsigned short tdc = 0;

		float chanpos[3]={0.,0.,0};
		fGeo->ChannelToPosition(channel, chanpos);

		for (auto const& edep : chanEDeps) {

			if (edep->TrackID() == gar::sdp::NoParticleId) continue;
        
			tdc = fClocks->TPCG4Time2TDC( edep->Time() +TMath::Abs(edep->X()-chanpos[0]) *fInverseVelocity);
			if ( tdc >= start_tdc && tdc <= stop_tdc) {
				float energy = edep->Energy();
				totalE += energy;
				int tid = edep->TrackID();
				auto tidmapiter = tidmap.find(tid);
				if (tidmapiter == tidmap.end()) {
					hitIDEs.emplace_back(tid,0,energy);
					tidmap[tid] = hitIDEs.size()-1;
				} else {
					hitIDEs.at(tidmapiter->second).numElectrons += energy;
				}
			}
		}

		// loop over the hitIDEs to set the fractional energy for each TrackID
		// and the total energy as well
		if (totalE < 1.e-5) totalE = 1.;

		for (size_t i = 0; i < hitIDEs.size(); ++i) {
			hitIDEs[i].energyFrac = hitIDEs[i].numElectrons / totalE;
		}
		return hitIDEs;
	}

    //----------------------------------------------------------------------
    double BackTrackerCore::HitCollectionPurity(simb::MCParticle* p,
		std::vector<art::Ptr<gar::rec::Hit> > const& hits,
		bool weightByCharge) const {

		int trackIDin = p->TrackId();


		float weight;
		float desired = 0.0;
		float total   = 0.0;
		for (auto hit : hits) {
			std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

			weight = (weightByCharge) ? hit->Signal() : 1.0;
			total += weight;
        
			for (auto iHitIDE : hitIDEs) {
				if (iHitIDE.trackID    == trackIDin &&
					iHitIDE.energyFrac >= fMinHitEnergyFraction) {
					desired += weight;
					break;
				}
			}
		} // end loop over hits
      
		double purity = 0.;
		if(total > 0) purity = desired/total;
		return purity;
	}

    //----------------------------------------------------------------------
    double BackTrackerCore::HitCollectionEfficiency(simb::MCParticle* p,
		std::vector<art::Ptr<gar::rec::Hit> > const& hits,
        std::vector<art::Ptr<gar::rec::Hit> > const& allhits,
        bool weightByCharge) const {

		int trackIDin = p->TrackId();

		float weight;
		float desired = 0.0;
		for (auto hit : hits) {
			std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

			weight = (weightByCharge) ? hit->Signal() : 1.0;

			// Don't double count if this hit has more than one of the
			// desired track IDs associated with it
			for (auto iHitIDE : hitIDEs) {
				if (iHitIDE.trackID    == trackIDin &&
					iHitIDE.energyFrac >= fMinHitEnergyFraction) {
					desired += weight;
					break;
				}
			}
		} // end loop over hits

		// Next figure out how many hits in the whole collection are associated with this id
		float total   = 0.0;
		for (auto hit : allhits) {
			std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

			weight = (weightByCharge) ? hit->Signal() : 1.0;
        
			for (auto iHitIDE : hitIDEs) {
				if (iHitIDE.trackID    == trackIDin &&
					iHitIDE.energyFrac >= fMinHitEnergyFraction) {
					total += weight;
					break;
				}
			}
        } // end loop over all hits

		double efficiency = 0.;
		if (total > 0.) efficiency = desired/total;
		return efficiency;
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
        
        hids = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());
        
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
      std::vector<HitIDE> ides = this->HitToHitIDEs(hit);
      
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
    std::set<int> BackTrackerCore::GetSetOfTrackIDs() const
    {
      // fParticleList::value_type is a pair (track, particle pointer)
      std::set<int> trackIDs;
      for (auto const& pl : fParticleList)
        trackIDs.insert(pl.first);
      
      return trackIDs;
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
