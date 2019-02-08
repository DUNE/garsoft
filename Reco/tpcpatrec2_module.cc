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
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Geometry/Geometry.h"

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

      std::string fHitLabel;        ///< label of module to get the hits from
      int fPrintLevel;              ///< debug printout:  0: none, 1: selected, 2: all
      float  fMaxVecHitLen;         ///< maximum vector hit length in patrec alg 2, in cm
      float  fVecHitRoad;           ///< max dist from a vector hit to a hit to assign it. for patrec alg 2.  in cm.
      float  fVecHitMatchCos;       ///< matching condition for pairs of vector hits cos angle between directions
      float  fVecHitMatchPos;       ///< matching condition for pairs of vector hits -- 3D distance (cm)
      float  fVecHitMatchPEX;       ///< matching condition for pairs of vector hits -- miss distance (cm)
      float  fVecHitMatchEta;       ///< matching condition for pairs of vector hits -- eta match (cm)
      float  fVecHitMatchLambda;    ///< matching condition for pairs of vector hits -- dLambda (radians)
      unsigned int    fVecHitMinHits;        ///< minimum number of hits on a vector hit for it to be considered
      unsigned int fInitialTPNHits; ///< number of hits to use for initial trackpar estimate, if present

      float  fHitRCut;              ///< only take hits within rcut of the center of the detector

      int initial_trackpar_estimate(const std::vector<Hit>  &hits,
				    std::vector<int> &hitlist,
				    float &curvature_init,
				    float &lambda_init,
				    float &phi_init,
				    float &xpos,
				    float &ypos,
				    float &zpos,
				    float &x_other_end);

      float capprox(float x1,float y1,
		    float x2,float y2,
		    float x3,float y3,
		    float &xc, float &yc);  ///< initial guess of curvature calculator -- from ALICE.  Also returns circle center

      typedef struct{
	TVector3 pos;
	TVector3 dir;
	std::vector<size_t> hitindex;
	float length;
      } vechit_t;

      bool vh_hitmatch(TVector3 &hpvec, int ihit, vechit_t &vechit, const std::vector<rec::Hit> &hits);
      void fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir);
      void fitline(std::vector<double> &x, std::vector<double> &y, double &lambda, double &intercept);
      bool vhclusmatch(std::vector<vechit_t> &cluster, vechit_t &vh);

      // rough estimate of track parameters

      int makepatrectrack(const std::vector<rec::Hit> &hits, std::vector<int> &thl, rec::TrackPar &trackpar);

    };


    tpcpatrec2::tpcpatrec2(fhicl::ParameterSet const& p) : EDProducer{p}  
      {
        fHitLabel          = p.get<std::string>("HitLabel","hit");
        fPrintLevel        = p.get<int>("PrintLevel",0);
	fMaxVecHitLen      = p.get<float>("MaxVecHitLen",10.0);
	fVecHitRoad        = p.get<float>("VecHitRoad",5.0);
	fVecHitMatchCos    = p.get<float>("VecHitMatchCos",0.9);
	fVecHitMatchPos    = p.get<float>("VecHitMatchPos",20.0);
	fVecHitMatchPEX    = p.get<float>("VecHitMatchPEX",5.0);
	fVecHitMatchEta    = p.get<float>("VecHitMatchEta",1.0);
	fVecHitMatchLambda = p.get<float>("VecHitMatchLambda",0.1);
	fVecHitMinHits     = p.get<unsigned int>("VecHitMinHits",3);
        fHitRCut             = p.get<float>("HitRCut",240);
        fInitialTPNHits    = p.get<int>("InitialTPNHits",100);

        art::InputTag hitTag(fHitLabel);
        consumes< std::vector<rec::Hit> >(hitTag);
        produces< std::vector<rec::Track> >();
        produces< std::vector<rec::VecHit> >();
        produces< art::Assns<rec::Hit, rec::VecHit> >();
        produces< art::Assns<rec::Hit, rec::Track> >();
        produces< art::Assns<rec::VecHit, rec::Track> >();
      }

    void tpcpatrec2::produce(art::Event& e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< std::vector<rec::VecHit> > vhCol(new std::vector<rec::VecHit>);
      std::unique_ptr< art::Assns<rec::Hit,rec::VecHit> > hitVHAssns(new ::art::Assns<rec::Hit,rec::VecHit>);
      std::unique_ptr< art::Assns<rec::VecHit,rec::Track> > vhTrkAssns(new ::art::Assns<rec::VecHit,rec::Track>);
      std::unique_ptr< art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      auto const trackPtrMaker = art::PtrMaker<rec::Track>(e);
      auto const vhPtrMaker = art::PtrMaker<rec::VecHit>(e);
      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, hitHandle.id());

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::ServiceHandle<geo::Geometry> geo;
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

      std::vector< std::vector<int> > hitlist;

      std::vector<vechit_t> vechits;
      std::vector<vechit_t> vhtmp;

      for (size_t ihitxs=0; ihitxs<hits.size(); ++ihitxs)
	{
	  int ihit = hsi[ihitxs];  // access hits sorted in X but keep original indices
	  const float *hpos = hits[ihit].Position();
	  TVector3 hpvec(hpos);
	  if ( ((hpos - tpccent).Cross(xhat)).Mag() > fHitRCut ) continue;  // skip hits if they are too far away from center as the
	  // last few pad rows may have distorted hits

	  bool matched=false;
	  for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
	    {
	      //std::cout << "testing hit: " << ihit << " with vector hit  " << ivh << std::endl;
	      if (vh_hitmatch(hpvec, ihit, vhtmp[ivh], hits ))  // updates vechit with this hit if matched
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
	      vh.length = 0;             // no length with just one hit
	      vh.hitindex.push_back(ihit);
	      vhtmp.push_back(vh);
	      //std::cout << "Created a new vector hit with one hit: " << hpos[0] << " " << hpos[1] << " " << hpos[2] << std::endl;
	    }
	}

      // trim the list of vechits down to only those with at least fVecHitMinHits

      for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
	{
	  if (vhtmp[ivh].hitindex.size() >= (size_t) fVecHitMinHits)
	    {
	      vechits.push_back(vhtmp[ivh]);
	    }
	}

      // stitch together vector hits into tracks
      // question -- do we need to iterate this, first looking for small-angle matches, and then
      // loosening up?  Also need to address tracks that get stitched across the primary vertex -- may need to split
      // these after other tracks have been found.

      std::vector< std::vector< vechit_t > > vhclusters;
      std::vector< std::vector< int > > vhcindex;

      for (size_t ivh = 0; ivh< vechits.size(); ++ ivh)
	{
	  if (fPrintLevel>1)
	    {
	      std::cout << " vhprint " << vechits[ivh].pos.X() << " " <<  vechits[ivh].pos.Y() << " " <<  vechits[ivh].pos.Z() << " " <<
		vechits[ivh].dir.X() << " " <<  vechits[ivh].dir.Y() << " " <<  vechits[ivh].dir.Z() << " " << vechits[ivh].hitindex.size() << std::endl;
	    }

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
	      std::vector<int> newci;
	      newci.push_back(ivh);
	      vhcindex.push_back(newci);
	    }
	  else if (clusmatchlist.size() == 1)
	    {
	      vhclusters[clusmatchlist[0]].push_back(vechits[ivh]);
	      vhcindex[clusmatchlist[0]].push_back(ivh);
	    }
	  else   // multiple matches -- merge clusters togetehr
	    {
	      for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		{
		  for (size_t ivh=0; ivh<vhclusters[clusmatchlist[icm]].size(); ++ivh)
		    {
		      vhclusters[clusmatchlist[0]].push_back(vhclusters[clusmatchlist[icm]][ivh]);
		      vhcindex[clusmatchlist[0]].push_back(vhcindex[clusmatchlist[icm]][ivh]);
		    }
		}
	      // remove the merged vh clusters, using the new indexes after removing earlier ones
	      for (size_t icm=1; icm<clusmatchlist.size(); ++icm)
		{
		  vhclusters.erase(vhclusters.begin() + (clusmatchlist[icm]-icm+1));
		  vhcindex.erase(vhcindex.begin() + (clusmatchlist[icm]-icm+1));
		}
	    }
	}

      // populate the hit list with hits from the vector hits in clusters.
      // sort them later

      std::vector<int> trkColIndex;
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
	  rec::TrackPar trackpar;
          if ( makepatrectrack(hits,hitlist[iclus],trackpar) == 0 )
	    {
	      trkCol->push_back(trackpar.CreateTrack());
	      trkColIndex.push_back(iclus);
	    }
	}

      // prepare data products for writing to the event

      for (size_t ivh=0; ivh<vechits.size(); ++ivh)
	{
	  vechit_t vh = vechits[ivh];
	  float pos[3] = {0,0,0};
	  float dir[3] = {0,0,0};
	  vh.pos.GetXYZ(pos);
	  vh.dir.GetXYZ(dir);
	  vhCol->emplace_back(pos,dir,vh.length);
	  auto const vhpointer = vhPtrMaker(ivh);
	  for (size_t ihit=0; ihit<vh.hitindex.size(); ++ihit)
	    {
	      auto const hitpointer = hitPtrMaker(vh.hitindex[ihit]);
	      hitVHAssns->addSingle(hitpointer,vhpointer);
	    }
	}

      for (size_t itrack=0; itrack<trkCol->size(); ++itrack)
	{
	  auto const trackpointer = trackPtrMaker(itrack);
	  for (size_t ihit=0; ihit<hitlist[trkColIndex[itrack]].size(); ++ihit)
	    {
	      auto const hitpointer = hitPtrMaker(hitlist[trkColIndex[itrack]][ihit]);
	      hitTrkAssns->addSingle(hitpointer,trackpointer);
	    }
	  for (size_t ivh=0; ivh<vhclusters[trkColIndex[itrack]].size(); ++ivh)
	    {
	      int ivh2 = vhcindex[trkColIndex[itrack]][ivh];
	      auto const vhpointer = vhPtrMaker(ivh2);
	      vhTrkAssns->addSingle(vhpointer,trackpointer);
	    }
	}

      e.put(std::move(trkCol));
      e.put(std::move(vhCol));
      e.put(std::move(hitVHAssns));
      e.put(std::move(vhTrkAssns));
      e.put(std::move(hitTrkAssns));

    }


    // see if a hit is consistent with a vector hit and add it if it is.
    // fit lines in y vs x and z vs x

    bool tpcpatrec2::vh_hitmatch(TVector3 &hpvec, int ihit, tpcpatrec2::vechit_t &vechit, const std::vector<rec::Hit> &hits)
    {
      bool retval = false;

      float dist = 1E10;
      float newlength = vechit.length;

      for (size_t iht=0; iht<vechit.hitindex.size(); ++iht)
	{
	  TVector3 ht(hits[vechit.hitindex[iht]].Position());
	  float d = (hpvec - ht).Mag();
	  if (d>fMaxVecHitLen) return retval;
	  dist = TMath::Min(dist,d);
	  newlength = TMath::Max(newlength,d);
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
	      hplist.emplace_back(hits[vechit.hitindex[i]].Position());
	    }
	  hplist.push_back(hpvec);
	  fitlinesdir(hplist,vechit.pos,vechit.dir);
	  vechit.hitindex.push_back(ihit);
	  vechit.length = newlength;  // we already calculated this above

	  //std::cout << "vechit now has " << hplist.size() << " hits" << std::endl;
	  //std::cout << vechit.pos.X() << " " << vechit.pos.Y() << " " << vechit.pos.Z() << std::endl;
	  //std::cout << vechit.dir.X() << " " << vechit.dir.Y() << " " << vechit.dir.Z() << std::endl;
	  retval = true;
	}
      return retval;
    }

    void tpcpatrec2::fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir)
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

    void tpcpatrec2::fitline(std::vector<double> &x, std::vector<double> &y, double &slope, double &intercept)
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

    bool tpcpatrec2::vhclusmatch(std::vector<tpcpatrec2::vechit_t> &cluster, vechit_t &vh)
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

    // rough estimate of track parameters -- both ends

    int tpcpatrec2::makepatrectrack(const std::vector<rec::Hit> &hits, std::vector<int> &thl, rec::TrackPar &trackpar)
    {
      // track parameters:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end

      // find candidate endpoints and sort hits

      float cmin[3];  // min x, y, and z coordinates over all hits
      float cmax[3];  // max x, y, and z coordinates over all hits
      size_t ihex[6];  // index of hit which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)

	 
      for (size_t ihit=0; ihit < thl.size(); ++ihit)
	{
	  for (int i=0; i<3; ++i)
	    {
	      float c = hits[thl[ihit]].Position()[i];
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
	  TVector3 poshc(hits[thl[ihex[i]]].Position());
	  for (size_t ihit=0; ihit<thl.size(); ++ihit)
	    {
	      TVector3 hp(hits[thl[ihit]].Position());
	      sumd += (poshc - hp).Mag();
	    }
	  if (sumd > sumdmax)
	    {
	      sumdmax = sumd;
	      imax = i;
	    }
	}

      //  Use this hit as a starting point -- find the closest hit to the last
      //  and add it to the newly sorted list hlf.  Change -- sort hits in order of how
      //  far they are from the first hit.  Prevents oscillations in position on sort order.

      float lengthforwards = 0;


      std::vector<int> hlf;
      hlf.push_back(thl[ihex[imax]]);
      TVector3 lpos(hits[hlf[0]].Position());
      TVector3 lastpos=lpos;
      TVector3 hpadd;

      for (size_t inh=1;inh<thl.size();++inh)
	{
	  float dmin=0;
	  float jmin=-1;
	  for (size_t jh=0;jh<thl.size();++jh)
	    {
	      bool found = false;
	      for (size_t kh=0;kh<hlf.size();++kh)
		{
		  if (hlf[kh] == thl[jh])
		    {
		      found = true;
		      break;
		    }
		}
	      if (found) continue;   // skip if we've already assigned this hit on this track
	      TVector3 hpos(hits[thl[jh]].Position());
	      float d=(hpos-lpos).Mag();
	      if (jmin == -1)
		{
		  jmin = jh;
		  dmin = d;
		  hpadd = hpos;
		}
	      else
		{
		  if (d<dmin)
		    {
		      jmin = jh;
		      dmin = d;
		      hpadd = hpos;
		    }
		}
	    }
	  //  std::cout << "dmin: " << dmin << std::endl;
	  hlf.push_back(thl[jmin]);
	  lengthforwards += (hpadd-lastpos).Mag();   // add up track length
	  lastpos = hpadd;
	}

      // replace our hit list with our newly sorted hit list.

      if (fPrintLevel>2)
	{
	  for (size_t ihit=0; ihit<thl.size(); ++ihit)
	    {
	      printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
		     thl[ihit],
		     hits[thl[ihit]].Position()[0],
		     hits[thl[ihit]].Position()[1],
		     hits[thl[ihit]].Position()[2],
		     hlf[ihit],
		     hits[hlf[ihit]].Position()[0],
		     hits[hlf[ihit]].Position()[1],
		     hits[hlf[ihit]].Position()[2]);
	    }
	}

      // now go backwards -- start at the end hit and use that as a starting point

      float lengthbackwards = 0;

      std::vector<int> hlb;
      hlb.push_back(hlf.back());
      TVector3 lpos2(hits[hlb[0]].Position());
      TVector3 lastpos2 = lpos2;

      for (size_t inh=1;inh<thl.size();++inh)
	{
	  float dmin=0;
	  float jmin=-1;
	  for (size_t jh=0;jh<thl.size();++jh)
	    {
	      bool found = false;
	      for (size_t kh=0;kh<hlb.size();++kh)
		{
		  if (hlb[kh] == thl[jh])
		    {
		      found = true;
		      break;
		    }
		}
	      if (found) continue;   // skip if we've already assigned this hit on this track
	      TVector3 hpos(hits[thl[jh]].Position());
	      float d=(hpos-lpos2).Mag();
	      if (jmin == -1)
		{
		  jmin = jh;
		  dmin = d;
		  hpadd = hpos;
		}
	      else
		{
		  if (d<dmin)
		    {
		      jmin = jh;
		      dmin = d;
		      hpadd = hpos;
		    }
		}
	    }
	  //  std::cout << "dmin: " << dmin << std::endl;
	  hlb.push_back(thl[jmin]);
	  lengthbackwards += (hpadd - lastpos2).Mag(); // add up length
	  lastpos2 = hpadd;
	}
      // replace our hit list with our newly sorted hit list.

      if (fPrintLevel>2)
	{
	  for (size_t ihit=0; ihit<thl.size(); ++ihit)
	    {
	      printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
		     thl[ihit],
		     hits[thl[ihit]].Position()[0],
		     hits[thl[ihit]].Position()[1],
		     hits[thl[ihit]].Position()[2],
		     hlb[ihit],
		     hits[hlb[ihit]].Position()[0],
		     hits[hlb[ihit]].Position()[1],
		     hits[hlb[ihit]].Position()[2]);
	    }
	}


      std::vector<float> tparbeg(6,0);
      float xother = 0;
      if ( initial_trackpar_estimate(hits, hlf, tparbeg[2], tparbeg[4], 
				     tparbeg[3], tparbeg[5], tparbeg[0], tparbeg[1], xother) != 0) 
	{
	  return 1;
	}

      std::vector<float> tparend(6,0);
      if ( initial_trackpar_estimate(hits, hlb, tparend[2], tparend[4], 
				     tparend[3], tparend[5], tparend[0], tparend[1], xother) != 0)
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

      trackpar.setNHits(thl.size());
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

  int tpcpatrec2::initial_trackpar_estimate(const std::vector<Hit>  &hits,
				    std::vector<int> &hitlist,
				    float &curvature_init,
				    float &lambda_init,
				    float &phi_init,
				    float &xpos,
				    float &ypos,
				    float &zpos,
				    float &x_other_end)
      {

      size_t nhits = hitlist.size();
      size_t firsthit = 0;
      size_t farhit = TMath::Min(nhits-1, (size_t) fInitialTPNHits);;
      size_t inthit = farhit/2;
      size_t lasthit = nhits-1;

      float trackbeg[3] = {hits[hitlist[firsthit]].Position()[0],
			   hits[hitlist[firsthit]].Position()[1],
			   hits[hitlist[firsthit]].Position()[2]};

      float tp1[3] = {hits[hitlist[inthit]].Position()[0],
		      hits[hitlist[inthit]].Position()[1],
		      hits[hitlist[inthit]].Position()[2]};

      float tp2[3] = {hits[hitlist[farhit]].Position()[0],
		      hits[hitlist[farhit]].Position()[1],
		      hits[hitlist[farhit]].Position()[2]};

      if (fPrintLevel>1)
	{
	  std::cout << "Hit Dump in initial_trackpar_estimate: " << std::endl;
	  for (size_t i=0;i<nhits;++i)
	    {
	      size_t ihf = i;
	      std::cout << i << " : " <<
		hits[hitlist[ihf]].Position()[0] << " " <<
		hits[hitlist[ihf]].Position()[1] << " " <<
		hits[hitlist[ihf]].Position()[2] << std::endl;
	    }
	}
      if (fPrintLevel>0)
	{
	  std::cout << "first hit: " << firsthit << ", inter hit: " << inthit << " " << " far hit: " << farhit << std::endl;
	  std::cout << "in the hit list: " << hitlist[firsthit] << " " << hitlist[inthit] << " " << hitlist[farhit] << std::endl;
	  std::cout << "First hit x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
	  std::cout << "Inter hit x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
	  std::cout << "Far   hit x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
	}

      xpos = trackbeg[0];
      ypos = trackbeg[1];
      zpos = trackbeg[2];
      x_other_end = hits[hitlist[lasthit]].Position()[0];


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

    float tpcpatrec2::capprox(float x1,float y1,
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


    DEFINE_ART_MODULE(tpcpatrec2)

  } // namespace rec
} // namespace gar
