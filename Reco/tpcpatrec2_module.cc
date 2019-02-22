////////////////////////////////////////////////////////////////////////
// Class:       tpcpatrec2
// Plugin Type: producer (art v3_00_00)
// File:        tpcpatrec2_module.cc
//
// Generated at Tue Feb  5 08:57:00 2019 by Thomas Junk using cetskelgen
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

namespace gar {
  namespace rec {

    class tpcpatrec2 : public art::EDProducer {
    public:
      explicit tpcpatrec2(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpcpatrec2(tpcpatrec2 const&) = delete;
      tpcpatrec2(tpcpatrec2&&) = delete;
      tpcpatrec2& operator=(tpcpatrec2 const&) = delete;
      tpcpatrec2& operator=(tpcpatrec2&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      // Declare member data here.

      std::string fVecHitLabel;     ///< label of module to get the vector hits and associations with hits
      int    fPrintLevel;           ///< debug printout:  0: none, 1: selected, 2: all
      float  fVecHitMatchCos;       ///< matching condition for pairs of vector hits cos angle between directions
      float  fVecHitMatchPos;       ///< matching condition for pairs of vector hits -- 3D distance (cm)
      float  fVecHitMatchPEX;       ///< matching condition for pairs of vector hits -- miss distance (cm)
      float  fVecHitMatchEta;       ///< matching condition for pairs of vector hits -- eta match (cm)
      float  fVecHitMatchLambda;    ///< matching condition for pairs of vector hits -- dLambda (radians)
      unsigned int fInitialTPNTPCClusters; ///< number of hits to use for initial trackpar estimate, if present
      size_t fMinNumTPCClusters;           ///< minimum number of hits for a patrec track

      // criteria for associating vector hits together to form clusters
      bool vhclusmatch(const std::vector<gar::rec::VecHit> &vechits, std::vector<size_t> &cluster, size_t vh);

      // rough estimate of track parameters
      int makepatrectrack(std::vector<gar::rec::TPCCluster> &hits, gar::rec::TrackPar &trackpar);

    };


    tpcpatrec2::tpcpatrec2(fhicl::ParameterSet const& p) : EDProducer{p}  
      {
        fVecHitLabel       = p.get<std::string>("VecHitLabel","vechit");
        fPrintLevel        = p.get<int>("PrintLevel",0);
	fVecHitMatchCos    = p.get<float>("VecHitMatchCos",0.9);
	fVecHitMatchPos    = p.get<float>("VecHitMatchPos",20.0);
	fVecHitMatchPEX    = p.get<float>("VecHitMatchPEX",5.0);
	fVecHitMatchEta    = p.get<float>("VecHitMatchEta",1.0);
	fVecHitMatchLambda = p.get<float>("VecHitMatchLambda",0.1);
        fInitialTPNTPCClusters    = p.get<unsigned int>("InitialTPNTPCClusters",100);
        fMinNumTPCClusters        = p.get<size_t>("MinNumTPCClusters",20);

        art::InputTag vechitTag(fVecHitLabel);
        consumes< std::vector<gar::rec::VecHit> >(vechitTag);
        consumes< art::Assns<gar::rec::TPCCluster, gar::rec::VecHit> >(vechitTag);
        produces< std::vector<gar::rec::Track> >();
        produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();
        produces< art::Assns<gar::rec::VecHit, gar::rec::Track> >();
      }

    void tpcpatrec2::produce(art::Event& e)
    {
      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::VecHit,gar::rec::Track> > vhTrkAssns(new ::art::Assns<gar::rec::VecHit,gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      auto vechitHandle = e.getValidHandle< std::vector<gar::rec::VecHit> >(fVecHitLabel);
      auto const& vechits = *vechitHandle;

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const vhPtrMaker = art::PtrMaker<gar::rec::VecHit>(e, vechitHandle.id());

      const art::FindManyP<gar::rec::TPCCluster> TPCClustersFromVecHits(vechitHandle,e,fVecHitLabel);

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);


      // stitch together vector hits into tracks
      // question -- do we need to iterate this, first looking for small-angle matches, and then
      // loosening up?  Also need to address tracks that get stitched across the primary vertex -- may need to split
      // these after other tracks have been found.

      std::vector< std::vector< size_t > > vhclusters;  // change this to use just indices of vh's

      for (size_t ivh = 0; ivh< vechits.size(); ++ ivh)
	{
	  if (fPrintLevel>1)
	    {
	      std::cout << " vhprint " << vechits[ivh].Position()[0] << " " <<  vechits[ivh].Position()[1] << " " <<  vechits[ivh].Position()[2] << " " <<
		vechits[ivh].Direction()[0] << " " <<  vechits[ivh].Direction()[1] << " " <<  vechits[ivh].Direction()[2] << " "  << std::endl;
	    }

	  std::vector<size_t> clusmatchlist;
	  for (size_t iclus=0; iclus<vhclusters.size(); ++iclus)
	    {
	      if (vhclusmatch(vechits,vhclusters[iclus],ivh))
		{
		  clusmatchlist.push_back(iclus);
		}
	    }
	  if (clusmatchlist.size() == 0)
	    {
	      std::vector<size_t> newclus;
	      newclus.push_back(ivh);
	      vhclusters.push_back(newclus);
	    }
	  else if (clusmatchlist.size() == 1)
	    {
	      vhclusters[clusmatchlist[0]].push_back(ivh);
	    }
	  else   // multiple matches -- merge clusters togetehr
	    {
	      for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		{
		  for (size_t ivh2=0; ivh2<vhclusters[clusmatchlist[icm]].size(); ++ivh2)
		    {
		      vhclusters[clusmatchlist[0]].push_back(vhclusters[clusmatchlist[icm]][ivh2]);
		    }
		}
	      // remove the merged vh clusters, using the new indexes after removing earlier ones
	      for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		{
		  vhclusters.erase(vhclusters.begin() + (clusmatchlist[icm]-icm+1));
		}
	    }
	}

      // make a local list of TPCClusters for each track and find initial track parameters
      for (size_t iclus=0; iclus < vhclusters.size(); ++iclus)
	{
	  std::vector<gar::rec::TPCCluster> TPCClusters;
	  std::vector<art::Ptr<gar::rec::TPCCluster> > TPCClusterptrs;
	  for (size_t ivh=0; ivh<vhclusters[iclus].size(); ++ivh)
	    {
	      for (size_t iTPCCluster=0; iTPCCluster<TPCClustersFromVecHits.at(vhclusters[iclus][ivh]).size(); ++iTPCCluster)
		{
		  TPCClusterptrs.push_back(TPCClustersFromVecHits.at(vhclusters[iclus][ivh]).at(iTPCCluster));
		  TPCClusters.push_back(*TPCClusterptrs.back());
		}
	    }

	  if (TPCClusters.size() >= fMinNumTPCClusters)
	    {
	      gar::rec::TrackPar trackpar;
               if ( makepatrectrack(TPCClusters,trackpar) == 0 )
	         {
	            trkCol->push_back(trackpar.CreateTrack());
		    auto const trackpointer = trackPtrMaker(trkCol->size()-1);
	            for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
		      {
			TPCClusterTrkAssns->addSingle(TPCClusterptrs.at(iTPCCluster),trackpointer);
		      }
		    for (size_t ivh=0; ivh<vhclusters[iclus].size(); ++ivh)
		      {
	                 auto const vhpointer = vhPtrMaker(ivh);
	                 vhTrkAssns->addSingle(vhpointer,trackpointer);
		      }
	         }
	    }
	}

      e.put(std::move(trkCol));
      e.put(std::move(vhTrkAssns));
      e.put(std::move(TPCClusterTrkAssns));

    }


    bool tpcpatrec2::vhclusmatch(const std::vector<gar::rec::VecHit> &vechits, std::vector<size_t> &cluster, size_t vhindex)
    {
      gar::rec::VecHit vh = vechits[vhindex];
      TVector3 vhdir(vh.Direction());
      TVector3 vhpos(vh.Position());

      for (size_t ivh=0; ivh<cluster.size(); ++ivh)
	{
	  gar::rec::VecHit vhtest = vechits[cluster[ivh]];
	  TVector3 vhtestdir(vhtest.Direction());
	  TVector3 vhtestpos(vhtest.Position());

	  //std::cout << "Testing vh " << ivh << " in a cluster of size: " << cluster.size() << std::endl;

	  // require the two VH's directions to point along each other -- use dot product

	  if (TMath::Abs((vhdir).Dot(vhtestdir)) < fVecHitMatchCos)
	    {
	      // std::cout << " Dot failure: " << TMath::Abs((vhdir).Dot(vhtestdir)) << std::endl;
	      continue;
	    }

	  // require the positions to be within fVecHitMatchPos of each other

	  if ((vhpos-vhtestpos).Mag() > fVecHitMatchPos)
	    {
	      //std::cout << " Pos failure: " << (vhpos-vhtestpos).Mag() << std::endl;
	      continue;
	    }

	  // require the extrapolation of one VH's line to another VH's center to match up.  Do for
	  // both VH's.

	  if ( ((vhpos-vhtestpos).Cross(vhdir)).Mag() > fVecHitMatchPEX )
	    {
	      //std::cout << "PEX failure: " << ((vhpos-vhtestpos).Cross(vhdir)).Mag() << std::endl;
	      continue;
	    }
	  if ( ((vhpos-vhtestpos).Cross(vhtestdir)).Mag() > fVecHitMatchPEX )
	    {
	      //std::cout << "PEX failure: " << ((vhpos-vhtestpos).Cross(vhtestdir)).Mag() << std::endl;
	      continue;
	    }

	  //--------------------------
	  // compute a 2D eta

	  // normalized direction vector for the VH under test, just the components
	  // perpendicular to X

	  TVector3 vhdp(vhdir);
	  vhdp.SetX(0);
	  float norm = vhdp.Mag();
	  if (norm > 0) vhdp *= (1.0/norm);

	  // same for the VH in the cluster under test

	  TVector3 vhcp(vhtestdir);
	  vhcp.SetX(0);
	  norm = vhcp.Mag();
	  if (norm > 0) vhcp *= (1.0/norm);

	  float relsign = 1.0;
	  if (vhdp.Dot(vhcp) < 0) relsign = -1;

	  TVector3 dcent = vhpos-vhtestpos;
	  dcent.SetX(0);

	  TVector3 avgdir1 = 0.5*(vhdp + relsign*vhcp);
	  float amag = avgdir1.Mag();
	  if (amag != 0) avgdir1 *= 1.0/amag;
	  float eta = (dcent.Cross(avgdir1)).Mag();

	  if ( eta > fVecHitMatchEta )
	    {
	      //std::cout << "Eta failure: " << eta1 << " " << eta2 << std::endl;
	      continue;
	    }

	  //----------------
	  // lambda requirement

	  float vhpd = TMath::Sqrt( TMath::Sq(vhdir.Y()) + TMath::Sq(vhdir.Z()) );
	  float vhxd = TMath::Abs( vhdir.X() );
	  float vhlambda = TMath::Pi()/2.0;
	  if (vhpd >0) vhlambda = TMath::ATan(vhxd/vhpd);

	  float cvhpd = TMath::Sqrt( TMath::Sq(vhtestdir.Y()) + TMath::Sq(vhtestdir.Z()) );
	  float cvhxd = TMath::Abs( vhtestdir.X() );
	  float cvhlambda = TMath::Pi()/2.0;
	  if (cvhpd >0) cvhlambda = TMath::ATan(cvhxd/cvhpd);

	  if ( TMath::Abs(vhlambda - cvhlambda) > fVecHitMatchLambda )
	    {
	      //std::cout << "dlambda  failure: " << vhlambda << " " << cvhlambda << std::endl;
	      continue;
	    }

	  if ( vhdir.Dot(vhtestdir) * vhdir.X() * vhtestdir.X() < 0 &&
	       TMath::Abs(vhdir.X()) > 0.01 && TMath::Abs(vhtestdir.X()) > 0.01)
	    {
	      //std::cout << "lambda sign failure" << std::endl;
	      continue;
	    }

	  //std::cout << " vh cluster match " << std::endl;
	  return true;
	}
      return false;
    }

    // rough estimate of track parameters -- both ends

    int tpcpatrec2::makepatrectrack(std::vector<gar::rec::TPCCluster> &trackTPCClusters, gar::rec::TrackPar &trackpar)
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

      gar::rec::sort_TPCClusters_along_track(trackTPCClusters,hlf,hlb,fPrintLevel,lengthforwards,lengthbackwards);

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

    DEFINE_ART_MODULE(tpcpatrec2)

  } // namespace rec
} // namespace gar
