////////////////////////////////////////////////////////////////////////
// Class:       tracker1
// Plugin Type: producer (art v2_11_02)
// File:        tracker1_module.cc
//
// Generated at Thu Jun 14 15:47:04 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Sort hits in X and add nearby hits.
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
#include "cetlib_except/exception.h"
#include "art/Persistency/Common/PtrMaker.h"

#include <memory>
#include <set>

#include "TVectorF.h"
#include "TMatrix.h"
#include "TMath.h"
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TVector3.h"

#include "Geant4/G4ThreeVector.hh"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/MagneticField/MagneticField.h"

// GArSoft Includes
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Geometry/GeometryGAr.h"

namespace gar {
  namespace rec {

    class tracker1 : public art::EDProducer {
    public:
      explicit tracker1(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tracker1(tracker1 const &) = delete;
      tracker1(tracker1 &&) = delete;
      tracker1 & operator = (tracker1 const &) = delete;
      tracker1 & operator = (tracker1 &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      float  fHitRCut;              ///< only take hits within rcut of the center of the detector
      size_t fPatRecAlg;            ///< 1: x-sorted patrec.  2: vector-hit patrec
      size_t fPatRecLookBack1;      ///< n hits to look backwards to make a linear extrapolation
      size_t fPatRecLookBack2;      ///< extrapolate from lookback1 to lookback2 and see how close the new hit is to the line
      float  fHitResolYZ;           ///< resolution in cm of a hit in YZ (pad size)
      float  fHitResolX;            ///< resolution in cm of a hit in X (drift direction)
      float  fSigmaRoad;            ///< how many sigma away from a track a hit can be and still add it during patrec

      float  fMaxVecHitLen;         ///< maximum vector hit length in patrec alg 2, in cm
      float  fVecHitRoad;           ///< max dist from a vector hit to a hit to assign it. for patrec alg 2.  in cm.
      float  fVecHitMatchCos;       ///< matching condition for pairs of vector hits cos angle between directions
      float  fVecHitMatchPos;       ///< matching condition for pairs of vector hits -- 3D distance (cm)
      float  fVecHitMatchPEX;       ///< matching condition for pairs of vector hits -- miss distance (cm)
      float  fVecHitMatchEta;       ///< matching condition for pairs of vector hits -- eta match (cm)
      float  fVecHitMatchLambda;    ///< matching condition for pairs of vector hits -- dLambda (radians)
      unsigned int    fVecHitMinHits;        ///< minimum number of hits on a vector hit for it to be considered

      float  fKalCurvStepUncSq;     ///< constant uncertainty term on each step of the Kalman fit -- squared, for curvature
      float  fKalPhiStepUncSq;      ///< constant uncertainty term on each step of the Kalman fit -- squared, for phi
      float  fKalLambdaStepUncSq;   ///< constant uncertainty term on each step of the Kalman fit -- squared, for lambda

      //float fXGapToEndTrack;      ///< how big a gap must be before we end a track and start a new one (unused for now)
      unsigned int fMinNumHits;     ///< minimum number of hits to define a track
      unsigned int fInitialTPNHits; ///< number of hits to use for initial trackpar estimate, if present
      std::string fHitLabel;        ///< label of module creating hits
      int fPrintLevel;              ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      int fTrackPass;               ///< which pass of the tracking to save as the tracks in the event
      int fDumpTracks;              ///< 0: do not print out tracks, 1: print out tracks

      std::string fSortOrder;      ///< switch to tell what way to sort hits before presenting them to the fitter

      float fHitResolYZinFit;      ///< Hit resolution parameter to use in fit
      float fRoadYZinFit;          ///< cut in cm for dropping hits from tracks in fit

      std::string fFirstPassFitType; ///< helix or Kalman -- which fitter to call for first-pass tracks
      std::string fSecondPassFitType; ///< helix or Kalman -- which fitter to call for second-pass tracks

      int initial_trackpar_estimate(art::ValidHandle<std::vector<Hit> > &hitHandle,
				    std::vector<std::vector<int> >      &hitlist,
				    std::vector<int>                    &hsi,
				    int itrack,
				    bool isForwards,
				    float &curvature_init,
				    float &lambda_init,
				    float &phi_init,
				    float &xpos,
				    float &ypos,
				    float &zpos,
				    float &x_other_end);

      int KalmanFit( art::ValidHandle<std::vector<Hit> > &hitHandle,
		     std::vector<std::vector<int> > &hitlist,
		     std::vector<int> &hsi,
		     int itrack,
		     bool isForwards,
		     std::vector<float> &trackparatend,
		     float &chisquared,
		     float &length,
		     float *covmat,    // 5x5 covariance matrix
		     std::set<int> &unused_hits);

      int KalmanFitBothWays(art::ValidHandle<std::vector<Hit> > &hitHandle,
			    std::vector<std::vector<int> > &hitlist,
			    std::vector<int> &hsi,
			    int itrack,
			    std::set<int> &unused_hits,
			    TrackPar &trackpar
			    );

      int FitHelix(art::ValidHandle<std::vector<Hit> > &hitHandle,
		   std::vector<std::vector<int> > &hitlist,
		   std::vector<int> &hsi,
		   int itrack,
		   bool isForwards,
		   std::set<int> &unused_hits,
		   TrackPar &trackpar
		   );

      size_t ifob(size_t ihit, size_t nhits, bool isForwards);

      float capprox(float x1,float y1,
		    float x2,float y2,
		    float x3,float y3,
		    float &xc, float &yc);  ///< initial guess of curvature calculator -- from ALICE.  Also returns circle center

      float capprox2(float y0, float z0, float y1, float z1, float y2, float z2);  //  -- returns abs value of curvature

      typedef struct{
	TVector3 pos;
	TVector3 dir;
	std::vector<size_t> hitindex;
      } vechit_t;

      bool vh_hitmatch(TVector3 &hpvec, int ihit, vechit_t &vechit, const std::vector<gar::rec::Hit> &hits, std::vector<int> &hsi);
      void fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir);
      void fitline(std::vector<double> &x, std::vector<double> &y, double &lambda, double &intercept);
      bool vhclusmatch(std::vector<vechit_t> &cluster, vechit_t &vh);

    };

    // constructor

    tracker1::tracker1(fhicl::ParameterSet const & p)
    {

      fHitRCut             = p.get<float>("HitRCut",240);
      fPatRecAlg           = p.get<size_t>("PatRecAlg",2);
      fPatRecLookBack1     = p.get<size_t>("PatRecLookBack1",5);
      fPatRecLookBack2     = p.get<size_t>("PatRecLookBack2",10);
      if (fPatRecLookBack1 == fPatRecLookBack2)
	{
	  throw cet::exception("tracker1_module: PatRecLookBack1 and PatRecLookBack2 are the same");
	}
      fHitResolYZ        = p.get<float>("HitResolYZ",1.0); // TODO -- think about what this value is
      fHitResolX         = p.get<float>("HitResolX",0.5);  // this is probably much better
      fSigmaRoad         = p.get<float>("SigmaRoad",5.0);
      fMinNumHits        = p.get<unsigned int>("MinNumHits",20);
      fHitLabel          = p.get<std::string>("HitLabel","hit");
      fPrintLevel        = p.get<int>("PrintLevel",0);
      fTrackPass         = p.get<int>("TrackPass",2);
      fDumpTracks        = p.get<int>("DumpTracks",2);
      fHitResolYZinFit   = p.get<float>("HitResolYZinFit",4.0);
      fRoadYZinFit       = p.get<float>("RoadYZinFit",1.0);
      fFirstPassFitType  = p.get<std::string>("FirstPassFitType","helix");
      fSecondPassFitType = p.get<std::string>("SecondPassFitType","Kalman");
      fMaxVecHitLen      = p.get<float>("MaxVecHitLen",10.0);
      fVecHitRoad        = p.get<float>("VecHitRoad",5.0);
      fVecHitMatchCos    = p.get<float>("VecHitMatchCos",0.9);
      fVecHitMatchPos    = p.get<float>("VecHitMatchPos",20.0);
      fVecHitMatchPEX    = p.get<float>("VecHitMatchPEX",5.0);
      fVecHitMatchEta    = p.get<float>("VecHitMatchEta",1.0);
      fVecHitMatchLambda = p.get<float>("VecHitMatchLambda",0.1);
      fVecHitMinHits     = p.get<unsigned int>("VecHitMinHits",3);

      fSortOrder         = p.get<std::string>("SortOrder","AlongLength");
      fInitialTPNHits    = p.get<int>("InitialTPNHits",100);

      fKalCurvStepUncSq  = p.get<float>("KalCurvStepUncSq",1.0E-9);
      fKalPhiStepUncSq   = p.get<float>("KalPhiStepUncSq",1.0E-9);
      fKalLambdaStepUncSq = p.get<float>("KalLambdaStepUncSq",1.0E-9);

      art::InputTag itag(fHitLabel);
      consumes< std::vector<gar::rec::Hit> >(itag);
      produces< std::vector<gar::rec::Track> >();
      produces< art::Assns<gar::rec::Hit, gar::rec::Track> >();
    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::Hit,gar::rec::Track> > hitTrkAssns(new ::art::Assns<gar::rec::Hit,gar::rec::Track>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const hitPtrMaker = art::PtrMaker<gar::rec::Hit>(e, hitHandle.id());

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::ServiceHandle<geo::GeometryGAr> geo;
      double xtpccent = geo->TPCXCent();
      double ytpccent = geo->TPCYCent();
      double ztpccent = geo->TPCZCent();
      TVector3 tpccent(xtpccent,ytpccent,ztpccent);
      TVector3 xhat(1,0,0);

      // make an array of hit indices sorted by hit X position
      std::vector<float> hitx;
      for (size_t i=0; i<hits.size(); ++i)
	{
	  hitx.push_back(hits[i].Position()[0]);
	}
      std::vector<int> hsi(hitx.size());
      TMath::Sort((int) hitx.size(),hitx.data(),hsi.data());

      float roadsq = fSigmaRoad*fSigmaRoad;

      // record which hits we have assigned to which tracks.  Make a forwards and a backwards hit list
      // as the sorting is different

      std::vector< std::vector<int> > hitlist;

      float resolSq = fHitResolYZ*fHitResolYZ;

      // initial patrec algorithm -- sort hits in x and look for clusters in y and z

      if (fPatRecAlg == 1)
	{

	  // do this twice, once going forwards through the hits, once backwards, and then collect hits
	  // in groups and split them when either the forwards or the backwards list says to split them.

	  std::vector< std::vector<int> > hitlistf;
	  std::vector<int> trackf(hits.size());
	  for (size_t ihit=0; ihit<hits.size(); ++ihit)
	    {
	      const float *hpos = hits[hsi[ihit]].Position();
	      TVector3 hpvec(hpos);
	      if ( ((hpos - tpccent).Cross(xhat)).Mag() > fHitRCut ) continue;  // skip hits if they are too far away from center as the
	      // last few pad rows may have distorted hits

	      float bestsignifs = -1;
	      int ibest = -1;
	      for (size_t itcand = 0; itcand < hitlistf.size(); ++itcand)
		{
		  float signifs = 1E9;
		  size_t hlsiz = hitlistf[itcand].size();
		  if (hlsiz > fPatRecLookBack1 && hlsiz > fPatRecLookBack2)
		    {
		      TVector3 pt1( hits[hsi[hitlistf[itcand][hlsiz-fPatRecLookBack1]]].Position() );
		      TVector3 pt2( hits[hsi[hitlistf[itcand][hlsiz-fPatRecLookBack2]]].Position() );
		      TVector3 uv = pt1-pt2;
		      uv *= 1.0/uv.Mag();
		      signifs = ((hpvec-pt1).Cross(uv)).Mag2()/resolSq;
		    }
		  else // not enough hits, just look how close we are to the last one
		    {
		      const float *cpos = hits[hsi[hitlistf[itcand].back()]].Position();
		      signifs = (TMath::Sq( (hpos[1]-cpos[1]) ) +
				 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
		    }
		  if (bestsignifs < 0 || signifs < bestsignifs)
		    {
		      bestsignifs = signifs;
		      ibest = itcand;
		    }
		}
	      if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
		{
		  ibest = hitlistf.size();
		  std::vector<int> vtmp;
		  hitlistf.push_back(vtmp);
		}
	      hitlistf[ibest].push_back(ihit);
	      trackf[ihit] = ibest;
	    }

	  std::vector< std::vector<int> > hitlistb;
	  std::vector<int> trackb(hits.size());
	  for (int ihit=hits.size()-1; ihit >= 0; --ihit)
	    {
	      const float *hpos = hits[hsi[ihit]].Position();
	      TVector3 hpvec(hpos);
	      if ( ((hpos - tpccent).Cross(xhat)).Mag() > fHitRCut ) continue;  // skip hits if they are too far away from center as the
	      // last few pad rows may have distorted hits

	      float bestsignifs = -1;
	      int ibest = -1;
	      for (size_t itcand = 0; itcand < hitlistb.size(); ++itcand)
		{
		  float signifs = 1E9;
		  size_t hlsiz = hitlistb[itcand].size();
		  if (hlsiz > fPatRecLookBack1 && hlsiz > fPatRecLookBack2)
		    {
		      TVector3 pt1( hits[hsi[hitlistb[itcand][hlsiz-fPatRecLookBack1]]].Position() );
		      TVector3 pt2( hits[hsi[hitlistb[itcand][hlsiz-fPatRecLookBack2]]].Position() );
		      TVector3 uv = pt1-pt2;
		      uv *= 1.0/uv.Mag();
		      signifs = ((hpvec-pt1).Cross(uv)).Mag2()/resolSq;
		    }
		  else // not enough hits, just look how close we are to the last one
		    {
		      const float *cpos = hits[hsi[hitlistb[itcand].back()]].Position();
		      signifs = (TMath::Sq( (hpos[1]-cpos[1]) ) +
				 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
		    }

		  if (bestsignifs < 0 || signifs < bestsignifs)
		    {
		      bestsignifs = signifs;
		      ibest = itcand;
		    }
		}
	      if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
		{
		  ibest = hitlistb.size();
		  std::vector<int> vtmp;
		  hitlistb.push_back(vtmp);
		}
	      hitlistb[ibest].push_back(ihit);
	      trackb[ihit] = ibest;
	    }

	  // make a list of tracks that is the set of disjoint subsets

	  for (size_t itrack=0; itrack<hitlistf.size(); ++itrack)
	    {
	      int itrackl = 0;
	      for (size_t ihit=0; ihit<hitlistf[itrack].size(); ++ihit)
		{
		  int ihif = hitlistf[itrack][ihit];
		  if (ihit == 0 || (trackb[ihif] != itrackl))
		    {
		      std::vector<int> vtmp;
		      hitlist.push_back(vtmp);
		      itrackl = trackb[ihif];
		    }
		  hitlist.back().push_back(ihif);
		}
	    }

	}

      // second try at a patrec algorithm -- find "vector hits"
      // start with hits sorted in x.  At least that limits how far through the list we must search
      // make a vector of vector hits, and then piece them together to form track candidates.

      else if (fPatRecAlg == 2)
	{
	  std::vector<vechit_t> vechits;
	  std::vector<vechit_t> vhtmp;

	  for (size_t ihit=0; ihit<hits.size(); ++ihit)
	    {
	      const float *hpos = hits[hsi[ihit]].Position();
	      TVector3 hpvec(hpos);
	      if ( ((hpos - tpccent).Cross(xhat)).Mag() > fHitRCut ) continue;  // skip hits if they are too far away from center as the
	      // last few pad rows may have distorted hits

	      bool matched=false;
	      for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
		{
		  //std::cout << "testing hit: " << ihit << " with vector hit  " << ivh << std::endl;
		  if (vh_hitmatch(hpvec, ihit, vhtmp[ivh], hits, hsi ))  // updates vechit with this hit if matched
		    {
		      matched = true;
		      break;
		    }
		}
	      if (!matched)   // make a new vechit if we haven't found one yet
		{
		  vechit_t vh;
		  vh.pos.SetXYZ(hpos[0],hpos[1],hpos[2]);
		  vh.dir.SetXYZ(0,0,0);      // new vechit with just one hit; don't know the direction yet
		  vh.hitindex.push_back(ihit);
		  vhtmp.push_back(vh);
		  //std::cout << "Created a new vector hit with one hit: " << hpos[0] << " " << hpos[1] << " " << hpos[2] << std::endl;
		}
	    }

	  // trim the list of vechits down to only those with more than two hits

	  for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
	    {
	      if (vhtmp[ivh].hitindex.size() >= (size_t)fVecHitMinHits)
		{
		  vechits.push_back(vhtmp[ivh]);
		}
	    }

	  // stitch together vector hits into tracks
	  // question -- do we need to iterate this, first looking for small-angle matches, and then
	  // loosening up?  Also need to address tracks that get stitched across the primary vertex -- may need to split
	  // these after other tracks have been found.

	  std::vector< std::vector< vechit_t > > vhclusters;

	  for (size_t ivh = 0; ivh< vechits.size(); ++ ivh)
	    {
	      //std::cout << " vhprint " << vechits[ivh].pos.X() << " " <<  vechits[ivh].pos.Y() << " " <<  vechits[ivh].pos.Z() << " " <<
	      //vechits[ivh].dir.X() << " " <<  vechits[ivh].dir.Y() << " " <<  vechits[ivh].dir.Z() << " " << vechits[ivh].hitindex.size() << std::endl;

	      std::vector<size_t> clusmatchlist;
	      for (size_t iclus=0; iclus<vhclusters.size(); ++iclus)
		{
		  if (vhclusmatch(vhclusters[iclus],vechits[ivh]))
		    {
		      clusmatchlist.push_back(iclus);
		    }
		}
	      if (clusmatchlist.size() == 0)
		{
		  std::vector<vechit_t> newclus;
		  newclus.push_back(vechits[ivh]);
		  vhclusters.push_back(newclus);
		}
	      else if (clusmatchlist.size() == 1)
		{
		  vhclusters[clusmatchlist[0]].push_back(vechits[ivh]);
		}
	      else   // multiple matches -- merge clusters togetehr
		{
		  for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		    {
		      for (size_t ivh=0; ivh<vhclusters[clusmatchlist[icm]].size(); ++ivh)
			{
			  vhclusters[clusmatchlist[0]].push_back(vhclusters[clusmatchlist[icm]][ivh]);
			}
		    }
		  // remove the merged vh clusters, using the new indexes after removing earlier ones
		  for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		    {
		      vhclusters.erase(vhclusters.begin() + (clusmatchlist[icm]-icm+1));
		    }
		}
	    }


	  // populate the hit list with hits from the vector hits in clusters.

	  for (size_t iclus=0; iclus < vhclusters.size(); ++iclus)
	    {
	      std::vector<int> vtmp;
	      hitlist.push_back(vtmp);
	      for (size_t ivh=0; ivh < vhclusters[iclus].size(); ++ ivh)
		{
		  for (size_t ihit=0; ihit < vhclusters[iclus][ivh].hitindex.size(); ++ihit)
		    {
		      hitlist.back().push_back(vhclusters[iclus][ivh].hitindex[ihit]);
		      //std::cout << "Added a hit " << hitlist.back().size() << " to track: " << iclus << std::endl;
		    }
		}
	      // re-sort the hitlist in x

	      std::sort(hitlist[iclus].begin(), hitlist[iclus].end(),
			[&hsi](int a, int b ) { return (hsi[a] > hsi[b]);});
	    }
	}

      else
	{
	  throw cet::exception("tracker1_module.cc: ununderstood PatRecAlg: ") << fPatRecAlg;
	}

      size_t ntracks = hitlist.size();

      // do a first pass of fitting the tracks

      std::vector<TrackPar> firstpass_tracks;
      std::vector<int> firstpass_tid;
      std::vector<TrackPar> secondpass_tracks;
      std::vector<int> secondpass_tid;

      if (fDumpTracks > 0)
	{
	  int ntracktmp = 0;
	  for (size_t itrack=0;itrack<ntracks;++itrack)
	    {
	      if (hitlist[itrack].size() >= fMinNumHits) ntracktmp++;
	    }
	  std::cout << "Trkdump: " << ntracktmp << std::endl;
	}

      std::vector<int> whichtrack(hits.size(),-1);

      for (size_t itrack=0; itrack<ntracks; ++itrack)
	{
	  size_t nhits = hitlist[itrack].size();
	  if ( nhits >= fMinNumHits)
	    {
	      for (size_t ihit=0; ihit<nhits; ++ihit)
		{
		  whichtrack[hitlist[itrack][ihit]] = itrack;  // fill this here so we only get the ones passing the nhits cut
		}

	      if (fPrintLevel)
		{
		  std::cout << "Starting a new Pass1 track: " << itrack << " Number of hits: " << nhits << std::endl;
		}

	      TrackPar trackparams;
	      std::set<int> unused_hits;

	      int retcode=0;
	      if (fFirstPassFitType == "helix")
		{
	          retcode = FitHelix(hitHandle,hitlist,hsi,itrack,true,unused_hits,trackparams);
		}
	      else if (fFirstPassFitType == "Kalman")
		{
	          retcode = KalmanFitBothWays(hitHandle,hitlist,hsi,itrack,unused_hits,trackparams);
		}
	      else
		{
		  throw cet::exception("Tracker1") << "Invalid first-pass fit type: " << fFirstPassFitType;
		}
	      if (retcode != 0) continue;

	      firstpass_tracks.push_back(trackparams);
	      firstpass_tid.push_back(itrack);

	      if (fDumpTracks > 0)
		{
		  std::cout << "Trkdump: " << itrack << std::endl;
		  std::cout << "Trkdump: " << trackparams.getTrackParametersBegin()[5] << std::endl;
		  for (int i=0; i<5;++i) std::cout << "Trkdump: " << trackparams.getTrackParametersBegin()[i] << std::endl;
		  std::cout << "Trkdump: " << trackparams.getTrackParametersEnd()[5] << std::endl;
		  for (int i=0; i<5;++i) std::cout << "Trkdump: " << trackparams.getTrackParametersEnd()[i] << std::endl;
		  std::cout << "Trkdump: " << nhits << std::endl;
		  for (size_t ihit=0;ihit<nhits;++ihit)
		    {
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[0] << std::endl;
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[1] << std::endl;
		      std::cout << "Trkdump: " << hits[hsi[hitlist[itrack][ihit]]].Position()[2] << std::endl;
		    }
		}
	    }
	}

      // Rearrange the hit lists -- ask ourselves which track each hit is best assigned to.  Make a new hit list, hitlist2
      // start only with hits that were assigned to tracks in pass1 (i.e. don't try to add new hits that made short, faraway tracks)
      // todo -- change this last requirement to an absolute distance cut.  Debug the track parameter extrapolation first.

      std::vector< std::vector<int> > hitlist2(firstpass_tracks.size());
      if (firstpass_tracks.size() > 0 && fTrackPass > 1)
	{
	  for (size_t ihit=0; ihit< hits.size(); ++ihit)
	    {
	      if (whichtrack[ihit] < 0) continue;
	      const float *hpos = hits[hsi[ihit]].Position();
	      float mindist = 0;
	      size_t ibest = 0;
	      for (size_t itrack=0; itrack<firstpass_tracks.size(); ++itrack)
		{
		  float dist = firstpass_tracks[itrack].DistXYZ(hpos);
		  if (itrack == 0 || dist < mindist)
		    {
		      mindist = dist;
		      ibest = itrack;
		    }
		}
	      hitlist2[ibest].push_back(ihit);
	    }

	  size_t ntracks2 = hitlist2.size();
	  for (size_t itrack=0; itrack<ntracks2; ++itrack)
	    {
	      size_t nhits = hitlist2[itrack].size();
	      if ( nhits >= fMinNumHits)
		{
		  if (fPrintLevel)
		    {
		      std::cout << "Starting a new Pass2 track: " << itrack << " Number of hits: " << nhits << std::endl;
		    }

		  TrackPar trackparams;
		  std::set<int> unused_hits;

		  int retcode=0;
		  if (fFirstPassFitType == "helix")
		    {
		      retcode = FitHelix(hitHandle,hitlist2,hsi,itrack,true,unused_hits,trackparams);
		    }
		  else if (fFirstPassFitType == "Kalman")
		    {
		      retcode = KalmanFitBothWays(hitHandle,hitlist2,hsi,itrack,unused_hits,trackparams);
		    }
		  else
		    {
		      throw cet::exception("Tracker1") << "Invalid first-pass fit type: " << fFirstPassFitType;
		    }
		  if (retcode != 0) continue;

		  secondpass_tracks.push_back(trackparams);
		  secondpass_tid.push_back(itrack);
		}
	    }
	}
      //  Remove stray hits.  Dig through unassociated hits and try to make extra tracks out of them.
      //  May need to wait until vertex finding is done so we know where to concentrate the effort

      // currently -- put second-pass tracks and associations with hits in the event

      if (fTrackPass == 1)
	{
	  for (size_t itrack=0; itrack<firstpass_tracks.size(); ++itrack)
	    {
	      trkCol->push_back(firstpass_tracks[itrack].CreateTrack());
	      auto const trackpointer = trackPtrMaker(itrack);
	      for (size_t ihit=0; ihit<hitlist[firstpass_tid[itrack]].size(); ++ ihit)
		{
		  auto const hitpointer = hitPtrMaker(hsi[hitlist[firstpass_tid[itrack]][ihit]]);
		  hitTrkAssns->addSingle(hitpointer,trackpointer);
		}
	    }
	}
      else if (fTrackPass == 2)
	{
	  for (size_t itrack=0; itrack<secondpass_tracks.size(); ++itrack)
	    {
	      trkCol->push_back(secondpass_tracks[itrack].CreateTrack());
	      auto const trackpointer = trackPtrMaker(itrack);
	      for (size_t ihit=0; ihit<hitlist2[secondpass_tid[itrack]].size(); ++ ihit)
		{
		  auto const hitpointer = hitPtrMaker(hsi[hitlist2[secondpass_tid[itrack]][ihit]]);
		  hitTrkAssns->addSingle(hitpointer,trackpointer);
		}
	    }
	}
      else
	{
	  throw cet::exception("Tracker1") << "Invalid track pass number requested: " << fTrackPass;
	}

      e.put(std::move(trkCol));
      e.put(std::move(hitTrkAssns));
    }

    //_____________________________________________________________________________
    float tracker1::capprox(float x1,float y1,
			    float x2,float y2,
			    float x3,float y3,
	                    float &xc, float &yc)
    {
      //-----------------------------------------------------------------
      // Initial approximation of the track curvature -- copied from ALICE
      // here x is y and y is z for us
      //-----------------------------------------------------------------
      x3 -=x1;
      x2 -=x1;
      y3 -=y1;
      y2 -=y1;
      //
      float det = x3*y2-x2*y3;
      if (TMath::Abs(det)<1e-10){
	return 100;
      }
      //
      float u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
      float x0 = x3*0.5-y3*u;
      float y0 = y3*0.5+x3*u;
      float c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
      xc = x0 + x1;
      yc = y0 + y1;
      if (det<0) c2*=-1;
      return c2;
    }

    // redo this using y and z as coordinates, and also return a positive value.  The caller needs to decide the
    // screw orientation of the helix -- not enough to go on with just y and z coordinates

    float tracker1::capprox2(float y0, float z0, float y1, float z1, float y2, float z2)
    {
      float A = y1*(z2 - z0) - z1*(y2 - y0) + (y2*z0 - z2*y0);

      float B = -( (y1*y1 + z1*z1)*(z2 - z0)
		   -z1*( (y2*y2 + z2*z2) - (y0*y0 + z0*z0) )
		   + ( (y2*y2 + z2*z2)*z0 - z2*(y0*y0 + z0*z0) ) );

      float C = (y1*y1 + z1*z1)*(y2 - y0)
	-y1*( (y2*y2 + z2*z2) - (y0*y0 + z0*z0))
	+ ( (y2*y2 + z2*z2)*y0 - y2*(y0*y0 + z0*z0) );

      float D = - ( (y1*y1 + z1*z1)*(y2*z0 - z2*y0)
		    - y1*( (y2*y2 + z2*z2)*z0 - z2*(y0*y0 + z0*z0) )
		    + z1*( (y2*y2 + z2*z2)*y0 - y2*(y0*y0 + z0*z0) ));

      if (TMath::Abs(A) < 1E-10)
	{ return 0;}

      float yc = -B/(2.0*A);
      float zc = -C/(2.0*A);
      float rs = yc*yc + zc*zc -D/A;

      if (fPrintLevel > 1)
	{
	  std::cout << " In capprox2, A, B, C, D, rs: " << A << " " << B << " " << C << " " << D << " " << rs << std::endl;
	}
      if (rs <= 0)
	{ throw cet::exception("tracker1capprox2: negative input to sqrt"); }
      float curv = 1.0/TMath::Sqrt(rs);
      return curv;
    }


    int tracker1::KalmanFitBothWays(art::ValidHandle<std::vector<Hit> > &hitHandle,
				    std::vector<std::vector<int> > &hitlist,
				    std::vector<int> &hsi,
				    int itrack,
				    std::set<int> &unused_hits,
				    TrackPar &trackpar
				    )

    {
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end

      // the "forward" fit is just in increasing x.  Track parameters are at the end of the fit


      // re-sort hits before giving them to the fitter.  Sorting in X may scramble the hit order if the track is in the Y,Z plane
      // start calling them tracks here.

      std::vector<std::vector<int> > hlf;
      std::vector<std::vector<int> > hlb;

      if (fSortOrder == "AlongLength")
	{

	  hlf = hitlist;  // only sort the track we're working on but pass the whole hitlist to the Kalman filter for historical reasons
	  hlb = hitlist;

          auto const& hits = *hitHandle;

	  float cmin[3];  // min x, y, and z coordinates over all hits
	  float cmax[3];  // max x, y, and z coordinates over all hits
	  size_t ihex[6];  // index of hit which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)

	  for (size_t ihit=0; ihit<hitlist[itrack].size(); ++ihit)
	    {
	      for (int i=0; i<3; ++i)
		{
		  float c = hits[hsi[hitlist[itrack][ihit]]].Position()[i];
		  if (ihit==0)
		    {
		      cmin[i] = c;
		      cmax[i] = c;
		      ihex[i] = 0;
		      ihex[i+3] = 0;
		    }
		  else
		    {
		      if (c<cmin[i])
			{
			  cmin[i] = c;
			  ihex[i] = ihit;
			}
		      if (c>cmax[i])
			{
			  cmax[i] = c;
			  ihex[i+3] = ihit;
			}
		    }
		}
	    }
	  // now we have six hits that have the min and max x, y, and z values.  Find out which of these six
	  // hits has the biggest sum of distances to all the other hits (the most extreme)
	  float sumdmax = 0;
	  size_t imax = 0;
	  for (size_t i=0; i<6; ++i)
	    {
	      float sumd = 0;
	      TVector3 poshc(hits[hsi[hitlist[itrack][ihex[i]]]].Position());
	      for (size_t ihit=0; ihit<hitlist[itrack].size(); ++ihit)
		{
		  TVector3 hp(hits[hsi[hitlist[itrack][ihit]]].Position());
		  sumd += (poshc - hp).Mag();
		}
	      if (sumd > sumdmax)
		{
		  sumdmax = sumd;
		  imax = i;
		}
	    }

	  //  Use this hit as a starting point -- find the closest hit to the last
	  //  and add it to the newly sorted list hls.  Change -- sort hits in order of how
	  //  far they are from the first hit.  Prevents oscillations in position on sort order.
	  //  This can be optimized to just sort an arry of distances using TMath::Sort.

	  std::vector<int> hls;
	  hls.push_back(hitlist[itrack][ihex[imax]]);
	  TVector3 lpos(hits[hsi[hls[0]]].Position());
	  for (size_t inh=1;inh<hitlist[itrack].size();++inh)
	    {
	      float dmin=0;
	      float jmin=-1;
	      for (size_t jh=0;jh<hitlist[itrack].size();++jh)
		{
		  bool found = false;
		  for (size_t kh=0;kh<hls.size();++kh)
		    {
		      if (hls[kh] == hitlist[itrack][jh])
			{
			  found = true;
			  break;
			}
		    }
		  if (found) continue;   // skip if we've already assigned this hit on this track
		  TVector3 hpos(hits[hsi[hitlist[itrack][jh]]].Position());
		  float d=(hpos-lpos).Mag();
		  if (jmin == -1)
		    {
		      jmin = jh;
		      dmin = d;
		    }
		  else
		    {
		      if (d<dmin)
			{
			  jmin = jh;
			  dmin = d;
			}
		    }
		}
	      //  std::cout << "dmin: " << dmin << std::endl;
	      hls.push_back(hitlist[itrack][jmin]);
	    }
	  // replace our hit list with our newly sorted hit list.

	  if (fPrintLevel>2)
	    {
	      std::cout << "Itrack: " << itrack << std::endl;
	      for (size_t ihit=0; ihit<hitlist[itrack].size(); ++ihit)
		{
		  printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
			 hitlist[itrack][ihit],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[0],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[1],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[2],
			 hls[ihit],
			 hits[hsi[hls[ihit]]].Position()[0],
			 hits[hsi[hls[ihit]]].Position()[1],
			 hits[hsi[hls[ihit]]].Position()[2]);
		}
	    }
	  hlf[itrack] = hls;

	  // now go backwards -- start at the end hit and use that as a starting point

	  hls.clear();
	  hls.push_back(hlf[itrack].back());
	  TVector3 lpos2(hits[hsi[hls[0]]].Position());
	  for (size_t inh=1;inh<hitlist[itrack].size();++inh)
	    {
	      float dmin=0;
	      float jmin=-1;
	      for (size_t jh=0;jh<hitlist[itrack].size();++jh)
		{
		  bool found = false;
		  for (size_t kh=0;kh<hls.size();++kh)
		    {
		      if (hls[kh] == hitlist[itrack][jh])
			{
			  found = true;
			  break;
			}
		    }
		  if (found) continue;   // skip if we've already assigned this hit on this track
		  TVector3 hpos(hits[hsi[hitlist[itrack][jh]]].Position());
		  float d=(hpos-lpos2).Mag();
		  if (jmin == -1)
		    {
		      jmin = jh;
		      dmin = d;
		    }
		  else
		    {
		      if (d<dmin)
			{
			  jmin = jh;
			  dmin = d;
			}
		    }
		}
	      //  std::cout << "dmin: " << dmin << std::endl;
	      hls.push_back(hitlist[itrack][jmin]);
	    }
	  // replace our hit list with our newly sorted hit list.

	  if (fPrintLevel>2)
	    {
	      std::cout << "Itrack: " << itrack << std::endl;
	      for (size_t ihit=0; ihit<hitlist[itrack].size(); ++ihit)
		{
		  printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
			 hitlist[itrack][ihit],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[0],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[1],
			 hits[hsi[hitlist[itrack][ihit]]].Position()[2],
			 hls[ihit],
			 hits[hsi[hls[ihit]]].Position()[0],
			 hits[hsi[hls[ihit]]].Position()[1],
			 hits[hsi[hls[ihit]]].Position()[2]);
		}
	    }
	  hlb[itrack] = hls;

	}
      else if (fSortOrder == "X")
	{
	  hlf = hitlist;
	  std::sort(hlf[itrack].begin(), hlf[itrack].end(),
		    [&hsi](int a, int b ) { return (hsi[a] > hsi[b]);});
	  hlb = hitlist;
	  std::sort(hlb[itrack].begin(), hlb[itrack].end(),
		    [&hsi](int a, int b ) { return (hsi[a] < hsi[b]);});
	}

      std::vector<float> tparend(6);
      float covmatend[25];
      float chisqforwards = 0;
      float lengthforwards = 0;
      int retcode = KalmanFit(hitHandle,hlf,hsi,itrack,true,tparend,chisqforwards,lengthforwards,covmatend,unused_hits);
      if (retcode != 0) return 1;

      // the "backwards" fit is in decreasing x.  Track paramters are at the end of the fit, the other end of the track

      std::vector<float> tparbeg(6);
      float covmatbeg[25];
      float chisqbackwards = 0;
      float lengthbackwards = 0;

      retcode = KalmanFit(hitHandle,hlb,hsi,itrack,false,tparbeg,chisqbackwards,lengthbackwards,covmatbeg,unused_hits);
      if (retcode != 0) return 1;

      size_t nhits=0;
      if (hitlist[itrack].size()>unused_hits.size())
	{ nhits = hitlist[itrack].size()-unused_hits.size(); }
      trackpar.setNHits(nhits);
      trackpar.setTime(0);
      trackpar.setChisqForwards(chisqforwards);
      trackpar.setChisqBackwards(chisqbackwards);
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

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    // KalmanFit does a forwards or backwards Kalman fit using the sorted hit list
    // variables:  x is the independent variable
    // 0: y
    // 1: z
    // 2: curvature
    // 3: phi
    // 4: lambda

    int tracker1::KalmanFit( art::ValidHandle<std::vector<Hit> > &hitHandle,
			     std::vector<std::vector<int> >      &hitlist,
			     std::vector<int>                    &hsi,
			     int itrack,
			     bool isForwards,
			     std::vector<float> &trackparatend,
			     float &chisquared,
			     float &length,
			     float *covmat,                     // 5x5 covariance matrix
			     std::set<int> &unused_hits)
    {

      // set some default values in case we return early

      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();
      chisquared = 0;
      length = 0;
      for (size_t i=0; i<5; ++i) trackparatend[i] = 0;
      for (size_t i=0; i<25; ++i) covmat[i] = 0;

      float roadsq = fRoadYZinFit*fRoadYZinFit;

      // estimate curvature, lambda, phi, xpos from the initial track parameters
      float curvature_init=0.1;
      float phi_init = 0;
      float lambda_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;
      float x_other_end = 0;
      if ( initial_trackpar_estimate(hitHandle,
				     hitlist,
				     hsi,
				     itrack,
				     isForwards,
				     curvature_init,
				     lambda_init,
				     phi_init,
				     xpos_init,
				     ypos_init,
				     zpos_init,
				     x_other_end) != 0)
	{
	  //std::cout << "kalman fit failed on initial trackpar estimate" << std::endl;
	  return 1;
	}

      // Kalman fitter variables

      float xpos = xpos_init;

      TMatrixF P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      P[0][0] = TMath::Sq(1); // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1); // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5); // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // lambda uncertainty

      TMatrixF PPred(5,5);

      // per-step additions to the covariance matrix
      TMatrixF Q(5,5);
      Q.Zero();
      Q[2][2] = fKalCurvStepUncSq;     // allow for some curvature uncertainty between points
      Q[3][3] = fKalPhiStepUncSq;      // phi
      Q[4][4] = fKalLambdaStepUncSq;   // lambda

      // uncertainties on the measured points  (big for now)
      TMatrixF R(2,2);
      R.Zero();
      R[0][0] = TMath::Sq(fHitResolYZinFit);  // in cm^2
      R[1][1] = TMath::Sq(fHitResolYZinFit);  // in cm^2

      // add the hits and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
      TMatrixF F(5,5);
      TMatrixF FT(5,5);
      TVectorF parvec(5);
      parvec[0] = ypos_init;
      parvec[1] = zpos_init;
      parvec[2] = curvature_init;
      parvec[3] = phi_init;
      parvec[4] = lambda_init;
      TVectorF predstep(5);

      TMatrixF H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
      TMatrixF HT(5,2);

      TVectorF z(2);
      TVectorF ytilde(2);
      TVectorF hx(2);
      TMatrixF S(2,2);
      TMatrixF K(5,2);

      TMatrixF I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;

      for (size_t ihit=1; ihit<nhits; ++ihit)
	{

	  size_t ihf = ifob(ihit,nhits,isForwards);
	  float xh = hits[hsi[hitlist[itrack][ihf]]].Position()[0];
	  float yh = hits[hsi[hitlist[itrack][ihf]]].Position()[1];
	  float zh = hits[hsi[hitlist[itrack][ihf]]].Position()[2];

	  if (fPrintLevel > 0)
	    {
	      std::cout << std::endl;
	      std::cout << "Adding a new hit: " << xh << " " << yh << " " << zh << std::endl;
	    }

	  // for readability

	  float curvature = parvec[2];
	  float phi = parvec[3];
	  float lambda = parvec[4];

	  // update prediction to the plane containing x.  Maybe we need to find the closest point on the helix to the hit we are adding,
	  // and not necessarily force it to be at this x

	  F.Zero();

	  // y = yold + slope*dx*Sin(phi).   F[0][i] = dy/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

	  float slope = TMath::Tan(lambda);
	  if (slope != 0)
	    {
	      slope = 1.0/slope;
	    }
	  else
	    {
	      slope = 1E9;
	    }

	  // relocate dx to be the location along the helix of the closest point.  Linearize for now near xpos.
	  // old calc

	  float dx = xh - xpos;

	  float dxdenom = slope*slope/(fHitResolYZ*fHitResolYZ) + 1.0/(fHitResolX*fHitResolX);
	  float dxnum = (slope/(fHitResolYZ*fHitResolYZ))*( (yh - parvec[0])*TMath::Sin(phi) + (zh - parvec[1])*TMath::Cos(phi) )
	    + (xh - xpos)/(fHitResolX*fHitResolX);
	  dx = dxnum/dxdenom;
	  if (dx == 0) dx = 1E-3;
	  //std::cout << "dxdenom, dxnum: " << dxdenom << " " << dxnum << std::endl;
	  //std::cout << "Track pos: " << xpos << " " << parvec[0] << " " << parvec[1] << " " << " Hit pos: " << xh << " " << yh << " " << zh << std::endl;
	  //std::cout << "dx old and new: " << xh - xpos << " " << dx << std::endl;


	  //TODO check this -- are these the derivatives?

	  // y = yold + dx*slope*TMath::Sin(phi)
	  // slope = cot(lambda), so dslope/dlambda = -csc^2(lambda) = -1 - slope^2
	  F[0][0] = 1.;
	  F[0][3] = dx*slope*TMath::Cos(phi);
	  F[0][4] = dx*TMath::Sin(phi)*(-1.0-slope*slope);

	  // z = zold + slope*dx*Cos(phi)
	  F[1][1] = 1.;
	  F[1][3] = -dx*slope*TMath::Sin(phi);
	  F[1][4] = dx*TMath::Cos(phi)*(-1.0-slope*slope);

	  // curvature = old curvature -- doesn't change but put in an uncertainty
	  F[2][2] = 1.;

	  // phi = old phi + curvature*slope*dx
	  // need to take the derivative of a product here
	  F[3][2] = dx*slope;
	  F[3][3] = 1.;
	  F[3][4] = dx*curvature*(-1.0-slope*slope);

	  // lambda -- same -- but put in an uncertainty in case it changes
	  F[4][4] = 1.;

	  // predicted step

	  if (fPrintLevel > 1)
	    {
	      std::cout << "F Matrix: " << std::endl;
	      F.Print();
	      std::cout << "P Matrix: " << std::endl;
	      P.Print();
	    }
	  if (fPrintLevel > 0)
	    {
	      std::cout << "x: " << xpos << " dx: " << dx <<  std::endl;
	      std::cout << " Parvec:   y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
	    }

	  predstep = parvec;
	  predstep[0] += slope*dx*TMath::Sin(phi);  // update y
	  predstep[1] += slope*dx*TMath::Cos(phi);  // update z
	  predstep[3] += slope*dx*curvature;        // update phi

	  if (fPrintLevel > 1)
	    {
	      std::cout << " Predstep: y " << predstep[0] << " z " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
	    }
	  // equations from the extended Kalman filter
	  FT.Transpose(F);
	  PPred = F*P*FT + Q;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "PPred Matrix: " << std::endl;
	      PPred.Print();
	    }

	  ytilde[0] = yh - predstep[0];
	  ytilde[1] = zh - predstep[1];
	  float ydistsq = ytilde.Norm2Sqr();
	  if (ydistsq > roadsq)
	    {
	      unused_hits.insert(ihf);
	      continue;
	    }
	  chisquared += ytilde.Norm2Sqr()/TMath::Sq(fHitResolYZ);
	  if (fPrintLevel > 0)
	    {
	      std::cout << "ytilde (residuals): " << std::endl;
	      ytilde.Print();
	    }
	  if (fPrintLevel > 1)
	    {
	      std::cout << "H Matrix: " << std::endl;
	      H.Print();
	    }

	  HT.Transpose(H);
	  S = H*PPred*HT + R;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "S Matrix: " << std::endl;
	      S.Print();
	    }

	  S.Invert();
	  if (fPrintLevel > 1)
	    {
	      std::cout << "Inverted S Matrix: " << std::endl;
	      S.Print();
	    }

	  K = PPred*HT*S;
	  if (fPrintLevel > 1)
	    {
	      std::cout << "K Matrix: " << std::endl;
	      K.Print();
	    }

	  float yprev = parvec[0];
	  float zprev = parvec[1];
	  parvec = predstep + K*ytilde;
	  P = (I-K*H)*PPred;
	  xpos = xpos + dx;
	  //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;

	  length += TMath::Sqrt( dx*dx + TMath::Sq(parvec[0]-yprev) + TMath::Sq(parvec[1]-zprev) );
	}

      for (size_t i=0; i<5; ++i)
	{
	  trackparatend[i] = parvec[i];
	}
      trackparatend[5] = xpos;  // tack this on so we can specify where the track endpoint is
      if (fPrintLevel > 1)
	{
	  std::cout << "Track params at end (y, z, curv, phi, lambda) " << trackparatend[0] << " " << trackparatend[1] << " " <<
	    trackparatend[2] << " " << trackparatend[3] <<" " << trackparatend[4] << std::endl;
	  S.Print();
	}

      // just for visualization of the initial track parameter guesses.  Comment out when fitting tracks

      //trackparatend[0] = ypos_init;
      //trackparatend[1] = zpos_init;
      //trackparatend[2] = curvature_init;
      //trackparatend[3] = phi_init;
      //trackparatend[4] = lambda_init;
      //trackparatend[5] = xpos_init;


      size_t icov=0;
      for (size_t i=0; i<5; ++i)
	{
	  for (size_t j=0; j<5; ++j)
	    {
	      covmat[icov] = P[i][j];
	    }
	}

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------


    size_t tracker1::ifob(size_t ihit, size_t nhits, bool isForwards)
    {

      return ihit;  // disable ifob -- hit list now encodes track direction

      //if (ihit >= nhits)
      //	{
      //  throw cet::exception("Tracker1") << "Invalid hit index in ifob: " << ihit << " nhits: " << nhits;
      //	}
      //if (isForwards)
      //{
      //  return ihit;
      //}
      //else
      //{
      //  return (nhits - ihit - 1);
      //	}
    }


    //--------------------------------------------------------------------------------------------------------------

    int tracker1::initial_trackpar_estimate(art::ValidHandle<std::vector<Hit> > &hitHandle,
					    std::vector<std::vector<int> >      &hitlist,
					    std::vector<int>                    &hsi,
					    int itrack,
					    bool isForwards,
					    float &curvature_init,
					    float &lambda_init,
					    float &phi_init,
					    float &xpos,
					    float &ypos,
					    float &zpos,
					    float &x_other_end)
    {
      // form a rough guess of track parameters

      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();

      size_t farhit_index = TMath::Min(nhits-1, (size_t) fInitialTPNHits);
      size_t inthit_index = farhit_index/2;

      size_t firsthit = ifob(0,nhits,isForwards);
      //size_t inthit = ifob(fMinNumHits/2,nhits,isForwards);
      //size_t farhit = ifob(fMinNumHits-1,nhits,isForwards);
      size_t inthit = ifob(inthit_index,nhits,isForwards);
      size_t farhit = ifob(farhit_index,nhits,isForwards);
      size_t lasthit = ifob(nhits-1,nhits,isForwards);

      float trackbeg[3] = {hits[hsi[hitlist[itrack][firsthit]]].Position()[0],
			   hits[hsi[hitlist[itrack][firsthit]]].Position()[1],
			   hits[hsi[hitlist[itrack][firsthit]]].Position()[2]};

      float tp1[3] = {hits[hsi[hitlist[itrack][inthit]]].Position()[0],
		      hits[hsi[hitlist[itrack][inthit]]].Position()[1],
		      hits[hsi[hitlist[itrack][inthit]]].Position()[2]};

      float tp2[3] = {hits[hsi[hitlist[itrack][farhit]]].Position()[0],
		      hits[hsi[hitlist[itrack][farhit]]].Position()[1],
		      hits[hsi[hitlist[itrack][farhit]]].Position()[2]};

      if (fPrintLevel>1)
	{
	  std::cout << "Hit Dump in initial_trackpar_estimate: " << std::endl;
	  for (size_t i=0;i<nhits;++i)
	    {
	      size_t ihf = ifob(i,nhits,isForwards);
	      std::cout << i << " : " <<
		hits[hsi[hitlist[itrack][ihf]]].Position()[0] << " " <<
		hits[hsi[hitlist[itrack][ihf]]].Position()[1] << " " <<
		hits[hsi[hitlist[itrack][ihf]]].Position()[2] << std::endl;
	    }
	}
      if (fPrintLevel>0)
	{
	  std::cout << "isForwards: " << isForwards << std::endl;
	  std::cout << "first hit: " << firsthit << ", inter hit: " << inthit << " " << " far hit: " << farhit << std::endl;
	  std::cout << "in the hit list: " << hsi[hitlist[itrack][firsthit]] << " " << hsi[hitlist[itrack][inthit]] << " " << hsi[hitlist[itrack][farhit]] << std::endl;
	  std::cout << "First hit x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
	  std::cout << "Inter hit x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
	  std::cout << "Far   hit x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
	}

      xpos = trackbeg[0];
      ypos = trackbeg[1];
      zpos = trackbeg[2];
      x_other_end = hits[hsi[hitlist[itrack][lasthit]]].Position()[0];


      float ycc=0;
      float zcc=0;
      curvature_init = capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2],ycc,zcc);
      //std::cout << " inputs to trackpar circ fit (y,z): " << trackbeg[1] << " " << trackbeg[2] << " : "
      //	    << tp1[1] << " " << tp1[2] << " : " << tp2[1] << " " << tp2[2] << std::endl;
      //std::cout << "curvature output: " << curvature_init << std::endl;

      phi_init = TMath::ATan2( trackbeg[2] - zcc, ycc - trackbeg[1] );
      float phi2 = phi_init;
      if (curvature_init<0) phi_init += TMath::Pi();
      float radius_init = 10000;
      if (curvature_init != 0) radius_init = 1.0/curvature_init;

      float dx1 = tp2[0] - xpos;
      if (dx1 != 0)
	{
	  float dphi2 = TMath::ATan2(tp2[2]-zcc,ycc-tp2[1])-phi2;
	  if (dphi2 > TMath::Pi()) dphi2 -= 2.0*TMath::Pi();
	  if (dphi2 < -TMath::Pi()) dphi2 += 2.0*TMath::Pi();
	  lambda_init = TMath::ATan(1.0/((radius_init/dx1)*dphi2));
	}
      else
	{
	  //std::cout << "initial track par estimate failure" << std::endl;
	  lambda_init = 0;
	  return 1;
	} // got fMinNumHits all at exactly the same value of x (they were sorted).  Reject track.

      if (fPrintLevel>0)
	{
	  std::cout << "phi calc: dz, dy " << tp2[2]-trackbeg[2] << " " <<  tp2[1]-trackbeg[1] << std::endl;
	  std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << phi_init << " " << lambda_init << std::endl;
	}
      return 0;
    }

    //--------------------------------------------------------------
    // the isForwards switch is only used to select which end to use to estimate the initial track parameters. Since this is a single helix
    // fit, the track parameters are the same, except for an evaluation of x, y, and z at the beginning and the end.
    // to do -- drop some hits if they have bad chisquared

    int tracker1::FitHelix(art::ValidHandle<std::vector<Hit> > &hitHandle,
			   std::vector<std::vector<int> > &hitlist,
			   std::vector<int> &hsi,
			   int itrack,
			   bool isForwards,
			   std::set<int> &unused_hits,
			   TrackPar &trackpar
			   )
    {
      auto const& hits = *hitHandle;
      size_t nhits = hitlist[itrack].size();
      if (nhits < fMinNumHits) return 1;

      // estimate curvature, lambda, phi, xpos from the initial track parameters
      float curvature_init=0.1;
      float phi_init = 0;
      float lambda_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;
      float x_other_end=0;
      if ( initial_trackpar_estimate(hitHandle,
				     hitlist,
				     hsi,
				     itrack,
				     isForwards,
				     curvature_init,
				     lambda_init,
				     phi_init,
				     xpos_init,
				     ypos_init,
				     zpos_init,
				     x_other_end) != 0)
	{
	  return 1;
	}

      float tpi[5] = {ypos_init, zpos_init, curvature_init, phi_init, lambda_init};
      float covmat[25] = {0};

      // syntax from $ROOTSYS/tutorials/fit/fitCircle.C

      auto chi2Function = [&](const Double_t *par) {
	//minimisation function computing the sum of squares of residuals
	// looping at the graph points

	float tpl[5] = { (float) par[0], (float) par[1], (float) par[2], (float) par[3], (float) par[4] };

	// only need this to compute chisquared, so set the track parameters the same at the beginning and end,
	// and no covmat is needed (defined outside to be zero)

	TrackPar tpar(0,0,nhits,xpos_init,tpl,covmat,0,xpos_init,tpl,covmat,0,0);

	float c2sum = 0;
	for (size_t ihit=0; ihit<nhits; ++ihit)
	  {
	    TVector3 hitpos(hits[hsi[hitlist[itrack][ihit]]].Position()[0],
			    hits[hsi[hitlist[itrack][ihit]]].Position()[1],
			    hits[hsi[hitlist[itrack][ihit]]].Position()[2]);
	    TVector3 helixpos = tpar.getPosAtX(hitpos.X(),isForwards);
	    //std::cout << hitpos.X() << " " << hitpos.Y() << " " << hitpos.Z() << " " << helixpos.X() << " " << helixpos.Y() << " " << helixpos.Z() << std::endl;
	    c2sum += (hitpos-helixpos).Mag2();
	  }
	double dc2sum = c2sum;
	return dc2sum;
      };

      // wrap chi2 funciton in a function object for the fit
      // 5 is the number of fit parameters (size of array par)
      ROOT::Math::Functor fcn(chi2Function,5);
      ROOT::Fit::Fitter  fitter;

      double pStart[5];
      for (size_t i=0; i<5; ++i) pStart[i] = tpi[i];
      fitter.SetFCN(fcn, pStart);
      fitter.Config().ParSettings(0).SetName("y0");
      fitter.Config().ParSettings(1).SetName("z0");
      fitter.Config().ParSettings(2).SetName("curvature");
      fitter.Config().ParSettings(3).SetName("phi0");
      fitter.Config().ParSettings(4).SetName("lambda");

      // do the fit
      bool ok = fitter.FitFCN();
      if (!ok) {
	LOG_WARNING("gar::rec::tracker1") << "Helix Fit failed";
	return(1);
      }

      const ROOT::Fit::FitResult & result = fitter.Result();
      //result.Print(std::cout);
      float fity0 = result.Value(0);
      float fitz0 = result.Value(1);
      float fitcurvature = result.Value(2);
      float fitphi0 = result.Value(3);
      float fitlambda = result.Value(4);
      float tpfit[5] = {fity0,fitz0,fitcurvature,fitphi0,fitlambda};
      float chisqmin = result.MinFcnValue();

      float stmp = TMath::Tan(fitlambda);
      //float stmp = TMath::Tan(fTrackParametersBegin[4]);
      if (stmp != 0)
	{
	  stmp = 1.0/stmp;
	}
      else
	{
	  stmp = 1E9;
	}

      float tracklength = TMath::Abs(xpos_init - x_other_end) * TMath::Sqrt( 1.0 + stmp*stmp );
      float covmatfit[25];
      for (size_t i=0; i<5; ++i)
	{
	  for (size_t j=0; j<5; ++j)
	    {
	      covmatfit[5*i + j] = result.CovMatrix(i,j);
	    }
	}

      trackpar.setNHits(nhits);
      trackpar.setTime(0);
      trackpar.setChisqForwards(chisqmin);
      trackpar.setChisqBackwards(chisqmin);
      trackpar.setLengthForwards(tracklength);
      trackpar.setLengthBackwards(tracklength);
      trackpar.setCovMatBeg(covmatfit);  // todo -- put in covariance matrices at both ends properly.
      trackpar.setCovMatEnd(covmatfit);
      if (isForwards)
	{
	  trackpar.setTrackParametersBegin(tpfit);
	  trackpar.setXBeg(xpos_init);
	  trackpar.setXEnd(x_other_end);
	  TVector3 xyzend = trackpar.getPosAtX(x_other_end,true);
	  float yend = xyzend[1];
	  float zend = xyzend[2];
          float tp_other_end[5] = {yend,zend,fitcurvature,fitphi0,fitlambda};
	  trackpar.setTrackParametersEnd(tp_other_end);
	}
      else
	{
	  trackpar.setTrackParametersEnd(tpfit);
	  trackpar.setXEnd(xpos_init);
	  trackpar.setXBeg(x_other_end);
	  TVector3 xyzend = trackpar.getPosAtX(x_other_end,true);
	  float yend = xyzend[1];
	  float zend = xyzend[2];
          float tp_other_end[5] = {yend,zend,fitcurvature,fitphi0,fitlambda};
	  trackpar.setTrackParametersBegin(tp_other_end);
	}

      return 0;
    }

    // see if a hit is consistent with a vector hit and add it if it is.
    // fit lines in y vs x and z vs x

    bool tracker1::vh_hitmatch(TVector3 &hpvec, int ihit, tracker1::vechit_t &vechit, const std::vector<gar::rec::Hit> &hits, std::vector<int> &hsi)
    {
      bool retval = false;

      float dist = 1E6;
      for (size_t iht=0; iht<vechit.hitindex.size(); ++iht)
	{
	  TVector3 ht(hits[hsi[vechit.hitindex[iht]]].Position());
	  float d = (hpvec - ht).Mag();
	  if (d>fMaxVecHitLen) return retval;
	  dist = TMath::Min(dist,d);
	}

      if (vechit.hitindex.size() > 1)
	{
          dist = ((hpvec - vechit.pos).Cross(vechit.dir)).Mag();
	  //std::cout << " Distance cross comparison: " << dist << std::endl;
	}

      if (dist < fVecHitRoad)  // add hit to vector hit if we have a match
	{
	  //std::cout << "matched a hit to a vh" << std::endl;
	  std::vector<TVector3> hplist;
	  for (size_t i=0; i< vechit.hitindex.size(); ++i)
	    {
	      hplist.emplace_back(hits[hsi[vechit.hitindex[i]]].Position());
	    }
	  hplist.push_back(hpvec);
	  fitlinesdir(hplist,vechit.pos,vechit.dir);
	  vechit.hitindex.push_back(ihit);
	  //std::cout << "vechit now has " << hplist.size() << " hits" << std::endl;
	  //std::cout << vechit.pos.X() << " " << vechit.pos.Y() << " " << vechit.pos.Z() << std::endl;
	  //std::cout << vechit.dir.X() << " " << vechit.dir.Y() << " " << vechit.dir.Z() << std::endl;
	  retval = true;
	}
      return retval;
    }

    void tracker1::fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir)
    {
      std::vector<double> x;
      std::vector<double> y;
      std::vector<double> z;
      for (size_t i=0; i<hlist.size();++i)
	{
	  x.push_back(hlist[i].X());
	  y.push_back(hlist[i].Y());
	  z.push_back(hlist[i].Z());
	}
      double slope_yx=0;
      double slope_zx=0;
      double intercept_yx=0;
      double intercept_zx=0;

      double slope_yz=0;
      double slope_xz=0;
      double intercept_yz=0;
      double intercept_xz=0;

      double slope_xy=0;
      double slope_zy=0;
      double intercept_xy=0;
      double intercept_zy=0;

      fitline(x,y,slope_xy,intercept_xy);
      fitline(x,z,slope_xz,intercept_xz);

      fitline(y,z,slope_yz,intercept_yz);
      fitline(y,x,slope_yx,intercept_yx);

      fitline(z,y,slope_zy,intercept_zy);
      fitline(z,x,slope_zx,intercept_zx);

      // pick the direction with the smallest sum of the absolute values of slopes to use to determine the line direction
      // in three-dimensional space

      double slopesumx = TMath::Abs(slope_xy) + TMath::Abs(slope_xz);
      double slopesumy = TMath::Abs(slope_yz) + TMath::Abs(slope_yx);
      double slopesumz = TMath::Abs(slope_zx) + TMath::Abs(slope_zy);

      if (slopesumx < slopesumy && slopesumx < slopesumz)
	{
	  dir.SetXYZ(1.0,slope_xy,slope_xz);
	  double avgx = 0;
	  for (size_t i=0; i<x.size(); ++i)
	    {
	      avgx += x[i];
	    }
	  avgx /= x.size();
          pos.SetXYZ(avgx, avgx*slope_xy + intercept_xy, avgx*slope_xz + intercept_xz);
	}
      else if (slopesumy < slopesumx && slopesumy < slopesumz)
	{
	  dir.SetXYZ(slope_yx,1.0,slope_yz);
	  double avgy = 0;
	  for (size_t i=0; i<y.size(); ++i)
	    {
	      avgy += y[i];
	    }
	  avgy /= y.size();
          pos.SetXYZ(avgy*slope_yx + intercept_yx, avgy, avgy*slope_yz + intercept_yz);
	}
      else
	{
	  dir.SetXYZ(slope_zx,slope_zy,1.0);
	  double avgz = 0;
	  for (size_t i=0; i<z.size(); ++i)
	    {
	      avgz += z[i];
	    }
	  avgz /= z.size();
          pos.SetXYZ(avgz*slope_zx + intercept_zx, avgz*slope_zy + intercept_zy, avgz);
	}
      dir *= 1.0/dir.Mag();

      // put in fit values for y and z
    }

    // fit with same weights on all points -- to think about: what are the uncertainties?

    void tracker1::fitline(std::vector<double> &x, std::vector<double> &y, double &slope, double &intercept)
    {
      size_t n = x.size();
      if (n < 2)
	{
	  throw cet::exception("tracker1: too few hits to fit a line in linefit");
	}
      double sumx = 0;
      double sumy = 0;
      double sumxx = 0;
      double sumxy = 0;

      for (size_t i=0; i<n; ++i)
	{
	  sumx += x[i];
	  sumy += y[i];
	  sumxx += TMath::Sq(x[i]);
	  sumxy += x[i]*y[i];
	}
      double denom = (n*sumxx) - TMath::Sq(sumx);
      if (denom == 0)
	{
	  slope = 1E6;   // is this right?
	  intercept = 0;
	}
      else
	{
	  slope = (n*sumxy - sumx*sumy)/denom;
	  intercept = (sumxx*sumy - sumx*sumxy)/denom;
	}
    }

    bool tracker1::vhclusmatch(std::vector<tracker1::vechit_t> &cluster, vechit_t &vh)
    {
      for (size_t ivh=0; ivh<cluster.size(); ++ivh)
	{
	  //std::cout << "Testing vh " << ivh << " in a cluster of size: " << cluster.size() << std::endl;

	  // require the two VH's directions to point along each other -- use dot product

	  if (TMath::Abs((vh.dir).Dot(cluster[ivh].dir)) < fVecHitMatchCos)
	    {
	      // std::cout << " Dot failure: " << TMath::Abs((vh.dir).Dot(cluster[ivh].dir)) << std::endl;
	      continue;
	    }

	  // require the positions to be within fVecHitMatchPos of each other

	  if ((vh.pos-cluster[ivh].pos).Mag() > fVecHitMatchPos)
	    {
	      //std::cout << " Pos failure: " << (vh.pos-cluster[ivh].pos).Mag() << std::endl;
	      continue;
	    }

	  // require the extrapolation of one VH's line to another VH's center to match up.  Do for
	  // both VH's.

	  if ( ((vh.pos-cluster[ivh].pos).Cross(vh.dir)).Mag() > fVecHitMatchPEX )
	    {
	      //std::cout << "PEX failure: " << ((vh.pos-cluster[ivh].pos).Cross(vh.dir)).Mag() << std::endl;
	      continue;
	    }
	  if ( ((vh.pos-cluster[ivh].pos).Cross(cluster[ivh].dir)).Mag() > fVecHitMatchPEX )
	    {
	      //std::cout << "PEX failure: " << ((vh.pos-cluster[ivh].pos).Cross(cluster[ivh].dir)).Mag() << std::endl;
	      continue;
	    }

	  //--------------------------
	  // compute a 2D eta

	  // normalized direction vector for the VH under test, just the components
	  // perpendicular to X

	  TVector3 vhdp(vh.dir);
	  vhdp.SetX(0);
	  float norm = vhdp.Mag();
	  if (norm > 0) vhdp *= (1.0/norm);

	  // same for the VH in the cluster under test

	  TVector3 vhcp(cluster[ivh].dir);
	  vhcp.SetX(0);
	  norm = vhcp.Mag();
	  if (norm > 0) vhcp *= (1.0/norm);

	  float relsign = 1.0;
	  if (vhdp.Dot(vhcp) < 0) relsign = -1;

	  TVector3 dcent = vh.pos-cluster[ivh].pos;
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

	  float vhpd = TMath::Sqrt( TMath::Sq(vh.dir.Y()) + TMath::Sq(vh.dir.Z()) );
	  float vhxd = TMath::Abs( vh.dir.X() );
	  float vhlambda = TMath::Pi()/2.0;
	  if (vhpd >0) vhlambda = TMath::ATan(vhxd/vhpd);

	  float cvhpd = TMath::Sqrt( TMath::Sq(cluster[ivh].dir.Y()) + TMath::Sq(cluster[ivh].dir.Z()) );
	  float cvhxd = TMath::Abs( cluster[ivh].dir.X() );
	  float cvhlambda = TMath::Pi()/2.0;
	  if (cvhpd >0) cvhlambda = TMath::ATan(cvhxd/cvhpd);

	  if ( TMath::Abs(vhlambda - cvhlambda) > fVecHitMatchLambda )
	    {
	      //std::cout << "dlambda  failure: " << vhlambda << " " << cvhlambda << std::endl;
	      continue;
	    }

          if ( vh.dir.Dot(cluster[ivh].dir) * vh.dir.X() * cluster[ivh].dir.X() < 0 &&
	       TMath::Abs(vh.dir.X()) > 0.01 && TMath::Abs(cluster[ivh].dir.X()) > 0.01)
	    {
	      //std::cout << "lambda sign failure" << std::endl;
	      continue;
	    }

	  //std::cout << " vh cluster match " << std::endl;
	  return true;
	}
      return false;
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar
