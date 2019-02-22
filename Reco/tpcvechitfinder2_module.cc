////////////////////////////////////////////////////////////////////////
// Class:       tpcvechitfinder2
// Plugin Type: producer (art v3_00_00)
// File:        tpcvechitfinder2_module.cc
//
// Generated at Tue Feb 12 15:44:27 2019 by Thomas Junk using cetskelgen
// A first step in GArSoft tracking pattern recognition
// finds short track segments from TPCClusters
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
#include "ReconstructionDataProducts/TPCCluster.h"
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
	std::vector<size_t> TPCClusterindex;
	float length;
	float c2ndf;
      } vechit_t;

      std::string fTPCClusterLabel;        ///< label of module to get the TPCClusters from
      int fPrintLevel;                     ///< debug printout:  0: none, 1: selected, 2: all
      float  fMaxVecHitLen;                ///< maximum vector hit length in patrec alg 2, in cm
      float  fVecHitRoad;                  ///< max dist from a vector hit to a TPCCluster to assign it. for patrec alg 2.  in cm.
      unsigned int    fVecHitMinTPCClusters;        ///< minimum number of TPCClusters on a vector hit for it to be considered
      float  fTPCClusterRCut;              ///< only take TPCClusters within rcut of the center of the detector
      size_t fNPasses;              ///< number of passes to reassign TPCClusters from vector hits with poor chisquareds
      float  fC2Cut;                ///< "chisquared" per ndof cut to reassign TPCClusters

      bool vh_TPCClusmatch(int iTPCCluster, vechit_t &vechit, const std::vector<gar::rec::TPCCluster> &TPCClusters, vechit_t &proposedvh);
      void fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir, float &chi2ndf);
      void fitline(std::vector<double> &x, std::vector<double> &y, double &lambda, double &intercept, double &chi2ndf);

    };


    tpcvechitfinder2::tpcvechitfinder2(fhicl::ParameterSet const& p) : EDProducer{p}
    {
        fTPCClusterLabel   = p.get<std::string>("TPCClusterLabel","tpccluster");
        fPrintLevel        = p.get<int>("PrintLevel",0);
	fMaxVecHitLen      = p.get<float>("MaxVecHitLen",10.0);
	fVecHitRoad        = p.get<float>("VecHitRoad",2.0);
	fVecHitMinTPCClusters     = p.get<unsigned int>("VecHitMinTPCClusters",3);
        fTPCClusterRCut           = p.get<float>("TPCClusterRCut",240);
	fNPasses           = p.get<unsigned int>("NPasses",2);
	fC2Cut             = p.get<float>("C2Cut",0.5);

	art::InputTag TPCClusterTag(fTPCClusterLabel);
        consumes< std::vector<gar::rec::TPCCluster> >(TPCClusterTag);
        produces< std::vector<gar::rec::VecHit> >();
        produces< art::Assns<gar::rec::TPCCluster, gar::rec::VecHit> >();

    }

    void tpcvechitfinder2::produce(art::Event& e)
    {
      std::unique_ptr< std::vector<gar::rec::VecHit> > vhCol(new std::vector<gar::rec::VecHit>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::VecHit> > TPCClusterVHAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::VecHit>);

      auto TPCClusterHandle = e.getValidHandle< std::vector<TPCCluster> >(fTPCClusterLabel);
      auto const& TPCClusters = *TPCClusterHandle;

      auto const vhPtrMaker = art::PtrMaker<gar::rec::VecHit>(e);
      auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());

      art::ServiceHandle<geo::Geometry> geo;
      double xtpccent = geo->TPCXCent();
      double ytpccent = geo->TPCYCent();
      double ztpccent = geo->TPCZCent();
      TVector3 tpccent(xtpccent,ytpccent,ztpccent);
      TVector3 xhat(1,0,0);

      // make an array of TPCCluster indices sorted by TPCCluster X position
      std::vector<float> TPCClusterx;
      for (size_t i=0; i<TPCClusters.size(); ++i)
	{
	  TPCClusterx.push_back(TPCClusters[i].Position()[0]);
	}
      std::vector<int> hsi(TPCClusterx.size());
      TMath::Sort((int) TPCClusterx.size(),TPCClusterx.data(),hsi.data());

      std::vector<vechit_t> vechits;
      std::vector<vechit_t> vhtmp;

      // iterate this, breaking up vechits with bad chisquareds and reassigning their TPCClusters

      for (size_t ipass=0; ipass<fNPasses; ++ipass)
	{
	  for (size_t iTPCClusterxs=0; iTPCClusterxs<hsi.size(); ++iTPCClusterxs)
	    {
	      int iTPCCluster = hsi[iTPCClusterxs];  // access TPCClusters sorted in X but keep original indices
	      const float *hpos = TPCClusters[iTPCCluster].Position();
	      TVector3 hpvec(hpos);
	      if ( ((hpvec - tpccent).Cross(xhat)).Mag() > fTPCClusterRCut ) continue;  // skip TPCClusters if they are too far away from center as the
	      // last few pad rows may have distorted TPCClusters

	      bool matched=false;
	      vechit_t proposedvh;
	      vechit_t bestproposedvh;
	      size_t ibestmatch = 0;

	      for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
		{
		  //std::cout << "testing TPCCluster: " << iTPCCluster << " with vector hit  " << ivh << std::endl;
		  if (vh_TPCClusmatch(iTPCCluster, vhtmp[ivh], TPCClusters, proposedvh ))  // updates vechit with this TPCCluster if matched
		    {
		      if (!matched || proposedvh.c2ndf < bestproposedvh.c2ndf )
			{
			  matched = true;
			  ibestmatch = ivh;
			  bestproposedvh = proposedvh;
			}
		    }
		}
	      if (!matched)   // make a new vechit if we haven't found one yet
		{
		  vechit_t vh;
		  vh.pos.SetXYZ(hpos[0],hpos[1],hpos[2]);
		  vh.dir.SetXYZ(0,0,0);      // new vechit with just one TPCCluster; don't know the direction yet
		  vh.length = 0;             // no length with just one TPCCluster
		  vh.c2ndf = 0;
		  vh.TPCClusterindex.push_back(iTPCCluster);
		  vhtmp.push_back(vh);
		  //std::cout << "Created a new vector hit with one TPCCluster: " << hpos[0] << " " << hpos[1] << " " << hpos[2] << std::endl;
		}
	      else
		{
		  vhtmp[ibestmatch] = bestproposedvh;
		}
	    }

	  // prepare for the next pass -- find vector hits with poor chisquareds and remove them from the list.  Make a new
	  // list instead of constantly removing from the existing list.  Repurpose hsi to contain the list of TPCClusters we want to reassign
	  // on the next pass.  

	  hsi.clear();
	  std::vector<vechit_t> vhtmp2;
	  for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
	    {
	      if (vhtmp[ivh].c2ndf < fC2Cut)
		{
		  vhtmp2.push_back(vhtmp[ivh]);
		}
	      else
		{
		  for (size_t iTPCCluster=0; iTPCCluster<vhtmp[ivh].TPCClusterindex.size(); ++iTPCCluster)
		    {
		      hsi.push_back(vhtmp[ivh].TPCClusterindex[iTPCCluster]);
		    }
		}
	    }
	  vhtmp = vhtmp2;
	  //std::cout << "left with " << vhtmp.size() << " vec hits and " << hsi.size() << " TPCClusters to reassociate " << std::endl;

	}

      // trim the list of vechits down to only those with at least fVecHitMinTPCClusters

      for (size_t ivh=0; ivh<vhtmp.size(); ++ivh)
	{
	  if (vhtmp[ivh].TPCClusterindex.size() >= (size_t) fVecHitMinTPCClusters)
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
	  for (size_t iTPCCluster=0; iTPCCluster<vh.TPCClusterindex.size(); ++iTPCCluster)
	    {
	      auto const TPCClusterpointer = TPCClusterPtrMaker(vh.TPCClusterindex[iTPCCluster]);
	      TPCClusterVHAssns->addSingle(TPCClusterpointer,vhpointer);
	    }
	}


      e.put(std::move(vhCol));
      e.put(std::move(TPCClusterVHAssns));
    }

    // see if a TPCCluster is consistent with a vector hit and add it if it is.
    // fit lines in y vs x and z vs x

    bool tpcvechitfinder2::vh_TPCClusmatch(int iTPCCluster, vechit_t &vechit, const std::vector<gar::rec::TPCCluster> &TPCClusters, vechit_t &proposedvh)
    {
      bool retval = false;

      TVector3 hpvec(TPCClusters.at(iTPCCluster).Position());
      float dist = 1E10;
      float newlength = vechit.length;

      for (size_t iht=0; iht<vechit.TPCClusterindex.size(); ++iht)
	{
	  TVector3 ht(TPCClusters.at(vechit.TPCClusterindex.at(iht)).Position());
	  float d = (hpvec - ht).Mag();
	  if (d>fMaxVecHitLen) return retval;
	  dist = TMath::Min(dist,d);
	  newlength = TMath::Max(newlength,d);
	}

      if (vechit.TPCClusterindex.size() > 1)
	{
	  dist = ((hpvec - vechit.pos).Cross(vechit.dir)).Mag();
	  //std::cout << " Distance cross comparison: " << dist << std::endl;
	}

      if (dist < fVecHitRoad)  // add TPCCluster to vector hit if we have a match
	{
	  //std::cout << "matched a TPCCluster to a vh" << std::endl;
	  std::vector<TVector3> hplist;
	  for (size_t i=0; i< vechit.TPCClusterindex.size(); ++i)
	    {
	      hplist.emplace_back(TPCClusters[vechit.TPCClusterindex[i]].Position());
	    }
	  hplist.push_back(hpvec);

	  proposedvh.TPCClusterindex = vechit.TPCClusterindex;
	  fitlinesdir(hplist,proposedvh.pos,proposedvh.dir,proposedvh.c2ndf);
	  proposedvh.TPCClusterindex.push_back(iTPCCluster);
	  proposedvh.length = newlength;  // we already calculated this above

	  //std::cout << "vechit now has " << hplist.size() << " TPCClusters" << std::endl;
	  //std::cout << proposedvh.pos.X() << " " << proposedvh.pos.Y() << " " << proposedvh.pos.Z() << std::endl;
	  //std::cout << proposedvh.dir.X() << " " << proposedvh.dir.Y() << " " << proposedvh.dir.Z() << std::endl;
	  retval = true;
	}
      return retval;
    }

    void tpcvechitfinder2::fitlinesdir(std::vector<TVector3> &hlist, TVector3 &pos, TVector3 &dir, float &chi2ndf)
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

      double chi2ndf_xy=0;
      double chi2ndf_xz=0;
      double chi2ndf_yz=0;
      double chi2ndf_yx=0;
      double chi2ndf_zx=0;
      double chi2ndf_zy=0;

      fitline(x,y,slope_xy,intercept_xy,chi2ndf_xy);
      fitline(x,z,slope_xz,intercept_xz,chi2ndf_xz);

      fitline(y,z,slope_yz,intercept_yz,chi2ndf_yz);
      fitline(y,x,slope_yx,intercept_yx,chi2ndf_yx);

      fitline(z,y,slope_zy,intercept_zy,chi2ndf_zy);
      fitline(z,x,slope_zx,intercept_zx,chi2ndf_zx);

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
	  chi2ndf = chi2ndf_xy + chi2ndf_xz;  // not really a chisquared per ndf, but perhaps good enough to pick which VH to assign a TPCCluster to
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
	  chi2ndf = chi2ndf_yx + chi2ndf_yz;  // not really a chisquared per ndf, but perhaps good enough to pick which VH to assign a TPCCluster to
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
	  chi2ndf = chi2ndf_zx + chi2ndf_zy;  // not really a chisquared per ndf, but perhaps good enough to pick which VH to assign a TPCCluster to
	}
      dir *= 1.0/dir.Mag();

      // put in fit values for y and z
    }

    // fit with same weights on all points -- to think about: what are the uncertainties?

    void tpcvechitfinder2::fitline(std::vector<double> &x, std::vector<double> &y, double &slope, double &intercept, double &chi2ndf)
    {
      chi2ndf = 0;
      size_t n = x.size();
      if (n < 2)
	{
	  throw cet::exception("tpcvechitfinder2_module.cc: too few TPCClusters to fit a line in linefit");
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
      if (n > 2)
	{
           for (size_t i=0; i<n; ++i)
	     {
	       chi2ndf += TMath::Sq( y[i] - slope*x[i] - intercept );  // no errors for now.
	     }
           chi2ndf /= (n-2);
	}
      else
	{
	  chi2ndf = 1.0; // a choice here -- if we can add the TPCCluster to an existing VH or make a "perfect-chisquare" 2-TPCCluster VH?
	  // define a 2-TPCCluster VH to be not perfect.
	}
    }

    DEFINE_ART_MODULE(tpcvechitfinder2)

  } // namespace rec
} // namespace gar
