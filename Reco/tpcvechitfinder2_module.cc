////////////////////////////////////////////////////////////////////////
// Class:       tpcvechitfinder2
// Plugin Type: producer (art v3_00_00)
// File:        tpcvechitfinder2_module.cc
//
// Generated at Tue Feb 12 15:44:27 2019 by Thomas Junk using cetskelgen
// A first step in GArSoft tracking pattern recognition
// finds short track segments from hits
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
// GArSoft Includes
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/VecHit.h"

#include "Geometry/Geometry.h"
#include "TMath.h"
#include "TVector3.h"

namespace gar {
  namespace rec {

    class tpcvechitfinder2 : public art::EDProducer {
    public:
      explicit tpcvechitfinder2(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpcvechitfinder2(tpcvechitfinder2 const&) = delete;
      tpcvechitfinder2(tpcvechitfinder2&&) = delete;
      tpcvechitfinder2& operator=(tpcvechitfinder2 const&) = delete;
      tpcvechitfinder2& operator=(tpcvechitfinder2&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

    // a struct used for convenience before making data products
    // not meant to be stored, but used by producers

      typedef struct{
	TVector3 pos;
	TVector3 dir;
	std::vector<size_t> hitindex;
	float length;
      } vechit_t;

      std::string fHitLabel;        ///< label of module to get the hits from
      int fPrintLevel;              ///< debug printout:  0: none, 1: selected, 2: all
      float  fMaxVecHitLen;         ///< maximum vector hit length in patrec alg 2, in cm
      float  fVecHitRoad;           ///< max dist from a vector hit to a hit to assign it. for patrec alg 2.  in cm.
      unsigned int    fVecHitMinHits;        ///< minimum number of hits on a vector hit for it to be considered
      float  fHitRCut;              ///< only take hits within rcut of the center of the detector

      bool vh_hitmatch(TVector3 &hpvec, int ihit, vechit_t &vechit, const std::vector<rec::Hit> &hits);
      void fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir);
      void fitline(std::vector<double> &x, std::vector<double> &y, double &lambda, double &intercept);

    };


    tpcvechitfinder2::tpcvechitfinder2(fhicl::ParameterSet const& p) : EDProducer{p}
    {
        fHitLabel          = p.get<std::string>("HitLabel","hit");
        fPrintLevel        = p.get<int>("PrintLevel",0);
	fMaxVecHitLen      = p.get<float>("MaxVecHitLen",10.0);
	fVecHitRoad        = p.get<float>("VecHitRoad",2.0);
	fVecHitMinHits     = p.get<unsigned int>("VecHitMinHits",3);
        fHitRCut             = p.get<float>("HitRCut",240);

	art::InputTag hitTag(fHitLabel);
        consumes< std::vector<rec::Hit> >(hitTag);
        produces< std::vector<rec::VecHit> >();
        produces< art::Assns<rec::Hit, rec::VecHit> >();

    }

    void tpcvechitfinder2::produce(art::Event& e)
    {
      std::unique_ptr< std::vector<rec::VecHit> > vhCol(new std::vector<rec::VecHit>);
      std::unique_ptr< art::Assns<rec::Hit,rec::VecHit> > hitVHAssns(new ::art::Assns<rec::Hit,rec::VecHit>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      auto const vhPtrMaker = art::PtrMaker<rec::VecHit>(e);
      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, hitHandle.id());

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


      e.put(std::move(vhCol));
      e.put(std::move(hitVHAssns));

    }


    // see if a hit is consistent with a vector hit and add it if it is.
    // fit lines in y vs x and z vs x

    bool tpcvechitfinder2::vh_hitmatch(TVector3 &hpvec, int ihit, vechit_t &vechit, const std::vector<rec::Hit> &hits)
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

    void tpcvechitfinder2::fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir)
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

    void tpcvechitfinder2::fitline(std::vector<double> &x, std::vector<double> &y, double &slope, double &intercept)
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

    DEFINE_ART_MODULE(tpcvechitfinder2)

  } // namespace rec
} // namespace gar
