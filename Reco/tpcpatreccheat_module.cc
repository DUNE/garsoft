////////////////////////////////////////////////////////////////////////
// Class:       tpcpatreccheat
// Plugin Type: producer (art v3_00_00)
// File:        tpcpatreccheat_module.cc
//
// copied from tpcpatrec2 and modified to just use the backtracker to cheat.  No vector hits are required
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>
#include <set>

#include "TMath.h"
#include "TVector3.h"

#include "Geant4/G4ThreeVector.hh"
#include "nutools/MagneticField/MagneticField.h"

// GArSoft Includes
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "MCCheater/BackTracker.h"
#include "Geometry/Geometry.h"

namespace gar {
  namespace rec {

    class tpcpatreccheat : public art::EDProducer {
    public:
      explicit tpcpatreccheat(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpcpatreccheat(tpcpatreccheat const&) = delete;
      tpcpatreccheat(tpcpatreccheat&&) = delete;
      tpcpatreccheat& operator=(tpcpatreccheat const&) = delete;
      tpcpatreccheat& operator=(tpcpatreccheat&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      // Declare member data here.

      std::string fHitLabel;     ///< label of module making hits
      std::string fTPCClusterLabel; /// label of module for input TPC clusters
      int    fPrintLevel;           ///< debug printout:  0: none, 1: selected, 2: all
      size_t fMinNumTPCClusters;           ///< minimum number of hits for a patrec track

      float  fSortTransWeight;      ///< for use in the hit sorting algorithm -- transverse distance weight factor
      float  fSortDistBack;         ///< for use in the hit sorting algorithm -- how far to go back before raising the distance figure of merit
      unsigned int fInitialTPNTPCClusters; ///< number of hits to use for initial trackpar estimate, if present

      // rough estimate of track parameters
      int makepatrectrack(std::vector<gar::rec::TPCCluster> &hits, gar::rec::TrackPar &trackpar);

    };


    tpcpatreccheat::tpcpatreccheat(fhicl::ParameterSet const& p) : EDProducer{p}  
      {
        fHitLabel           = p.get<std::string>("HitLabel","hit");
        fTPCClusterLabel    = p.get<std::string>("TPCClusterLabel","tpccluster");
        fPrintLevel         = p.get<int>("PrintLevel",0);
        fMinNumTPCClusters  = p.get<size_t>("MinNumTPCClusters",20);
        fInitialTPNTPCClusters    = p.get<unsigned int>("InitialTPNTPCClusters",100);
        fSortTransWeight    = p.get<float>("SortTransWeight",0.1);
        fSortDistBack       = p.get<float>("SortDistBack",2.0);

        art::InputTag hitTag(fHitLabel);
        art::InputTag TPCClusterTag(fTPCClusterLabel);

        consumes< std::vector<gar::rec::Hit> >(hitTag);
        consumes< std::vector<gar::rec::TPCCluster> >(TPCClusterTag);
        consumes< art::Assns<gar::rec::TPCCluster, gar::rec::Hit> >(TPCClusterTag);

        produces< std::vector<gar::rec::Track> >();
        produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();
        produces< art::Assns<gar::rec::VecHit, gar::rec::Track> >();
      }

    void tpcpatreccheat::produce(art::Event& e)
    {

      if( e.isRealData() ){
        throw cet::exception("tpcpatreccheat")
        << "Attempting to cheat the track pattern recognition for real data, "
        << "that will never work";
        return;
      }

      auto TPCClusterHandle = e.getValidHandle< std::vector<gar::rec::TPCCluster> >(fTPCClusterLabel);
      auto const& TPCClusters = *TPCClusterHandle;

      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::VecHit,gar::rec::Track> > vhTrkAssns(new ::art::Assns<gar::rec::VecHit,gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());
      const art::FindManyP<gar::rec::Hit> HitsFromTPCClusters(TPCClusterHandle,e,fTPCClusterLabel);

      auto bt = gar::providerFrom<cheat::BackTracker>();
      auto geo = gar::providerFrom<geo::Geometry>();
      float xtpccent = geo->TPCXCent();

      // Do everything separately for the two sides -- do not make a single track out of pieces on either side of the cathode

      for (int iside=-1; iside <= 1; iside += 2)
	{
	  std::map<int, std::set<size_t> > trackIDclusmap;
	  for (size_t iclus = 0; iclus<TPCClusters.size(); ++iclus)
	    {
	      float totcharge = 0;
	      std::map<int, float> chargemap;  // trackID to charge map -- charge is in electrons
	      auto hits = HitsFromTPCClusters.at(iclus);
	      for (size_t ihit=0; ihit< hits.size(); ++ihit)
		{
		  unsigned int channel = hits.at(ihit)->Channel();
		  float chanpos[3];
		  geo->ChannelToPosition(channel,chanpos);
		  int hitside = 1;
		  if (chanpos[0] < xtpccent)
		    {
		      hitside = -1;
		    }
		  if (hitside != iside) continue;
		  auto hitides = bt->HitToTrackID(hits.at(ihit));
		  //std::cout << "number of hitides: " << hitides.size() << std::endl;
		  for (size_t i=0; i<hitides.size(); ++i)
		    {
		      float charge = hitides.at(i).numElectrons;
		      totcharge += charge;
		      int trackid = TMath::Abs(hitides.at(i).trackID);
		      auto cmf = chargemap.find(trackid);
		      if (cmf == chargemap.end())
			{
			  chargemap[trackid]=charge;
			}
		      else
			{
			  cmf->second += charge;
			}
		    }
		}
	      if (totcharge > 0)
		{
		  float maxcharge = 0;
		  int ibesttrackid = 0;
		  for (auto const& cm : chargemap)
		    {
		      if (cm.second > maxcharge)
			{
			  maxcharge = cm.second;
			  ibesttrackid = cm.first;
			}
		    }
		  trackIDclusmap[ibesttrackid].emplace(iclus);
		}
	    } // end loop over TPC Clusters for grouping them by trackID

	  // loop over the map associating trackIDs to TPC clusters and for each entry
	  // make a local list of TPCClusters for each track and find initial track parameters

	  for (auto const& tic : trackIDclusmap)
	    {
  	      std::vector<gar::rec::TPCCluster> TPCClusters_thistrack;
	      std::vector<art::Ptr<gar::rec::TPCCluster> > TPCClusterptrs;
	      for (auto const& iclus : tic.second)
		{
		  TPCClusters_thistrack.push_back(TPCClusters.at(iclus));
		  TPCClusterptrs.push_back(TPCClusterPtrMaker(iclus));
		}


	      if (TPCClusters_thistrack.size() >= fMinNumTPCClusters)
		{
		  gar::rec::TrackPar trackpar;
		  if ( makepatrectrack(TPCClusters_thistrack,trackpar) == 0 )
		    {
		      trkCol->push_back(trackpar.CreateTrack());
		      auto const trackpointer = trackPtrMaker(trkCol->size()-1);
		      for (size_t iTPCCluster=0; iTPCCluster<TPCClusters_thistrack.size(); ++iTPCCluster)
			{
			  TPCClusterTrkAssns->addSingle(TPCClusterptrs.at(iTPCCluster),trackpointer);
			}
		    }
		}
	    }
	}
      e.put(std::move(trkCol));
      e.put(std::move(vhTrkAssns));
      e.put(std::move(TPCClusterTrkAssns));
    }


    // rough estimate of track parameters -- both ends

    int tpcpatreccheat::makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar)
    {
      // track parameters:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end


      float lengthforwards = 0;
      std::vector<int> hlf;
      float lengthbackwards = 0;
      std::vector<int> hlb;

      gar::rec::sort_TPCClusters_along_track(trackTPCClusters,hlf,hlb,fPrintLevel,lengthforwards,lengthbackwards,fSortTransWeight,fSortDistBack);

      std::vector<float> tparbeg(6,0);
      float xother = 0;
      if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlf, tparbeg[2], tparbeg[4], 
					       tparbeg[3], tparbeg[5], tparbeg[0], tparbeg[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0) 
	{
	  return 1;
	}

      std::vector<float> tparend(6,0);
      if ( gar::rec::initial_trackpar_estimate(trackTPCClusters, hlb, tparend[2], tparend[4], 
					       tparend[3], tparend[5], tparend[0], tparend[1], xother, fInitialTPNTPCClusters, fPrintLevel) != 0)
	{
	  return 1;
	}

      // no chisquare or covariance in patrec tracks
      float covmatbeg[25];
      float covmatend[25];
      for (size_t i=0; i<25; ++i) // no covmat in patrec tracks
	{
	  covmatend[i] = 0;
	  covmatbeg[i] = 0;
	}

      trackpar.setNTPCClusters(trackTPCClusters.size());
      trackpar.setTime(0);
      trackpar.setChisqForwards(0);
      trackpar.setChisqBackwards(0);
      trackpar.setLengthForwards(lengthforwards);
      trackpar.setLengthBackwards(lengthbackwards);
      trackpar.setCovMatBeg(covmatbeg);
      trackpar.setCovMatEnd(covmatend);
      trackpar.setTrackParametersBegin(tparbeg.data());
      trackpar.setXBeg(tparbeg[5]);
      trackpar.setTrackParametersEnd(tparend.data());
      trackpar.setXEnd(tparend[5]);

      return 0;
    }


    DEFINE_ART_MODULE(tpcpatreccheat)

    
  } // namespace rec
} // namespace gar
