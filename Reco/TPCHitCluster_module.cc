////////////////////////////////////////////////////////////////////////
// Class:       TPCHitCluster
// Plugin Type: producer (art v3_00_00)
// File:        TPCHitCluster_module.cc
//
// Takes Hits and clusters them together
//
// Generated at Fri Feb 22 13:16:03 2019 by Thomas Junk using cetskelgen
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

#include "TMath.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/TPCCluster.h"

namespace gar {

  namespace rec {

    class TPCHitCluster : public art::EDProducer {
    public:
      explicit TPCHitCluster(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      TPCHitCluster(TPCHitCluster const&) = delete;
      TPCHitCluster(TPCHitCluster&&) = delete;
      TPCHitCluster& operator=(TPCHitCluster const&) = delete;
      TPCHitCluster& operator=(TPCHitCluster&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      std::string fHitLabel;        ///< label of module to get the hits from
      int fClusterHits;         ///< hit clustering algorithm number
      float fHitClusterDx;      ///< range in cm to look for hits to cluster in x
      float fHitClusterDyDz;    ///< range in cm to look for hits to cluster in y and z
      int fPrintLevel;              ///< debug printout:  0: none, 1: selected, 2: all

    };


    TPCHitCluster::TPCHitCluster(fhicl::ParameterSet const& p) 
      : EDProducer{p}
    {
      fHitLabel         = p.get<std::string>("HitLabel","hit");
      fPrintLevel       = p.get<int>("PrintLevel",0);
      fClusterHits      = p.get<int>("ClusterHits",1);
      fHitClusterDx     = p.get<float>("HitClusterDx",1.0);
      fHitClusterDyDz   = p.get<float>("HitClusterDyDz",5.0);

      art::InputTag hitTag(fHitLabel);
      consumes< std::vector<gar::rec::Hit> >(hitTag);
      produces< std::vector<gar::rec::TPCCluster> >();
      produces< art::Assns<gar::rec::Hit, gar::rec::TPCCluster> >();
    }

    void TPCHitCluster::produce(art::Event& e)
    {

      // input: hits

      auto hitHandle = e.getValidHandle< std::vector<gar::rec::Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      // output: TPCCluster and associations with hits

      std::unique_ptr<std::vector<gar::rec::TPCCluster> > TPCClusterCol (new std::vector<gar::rec::TPCCluster> );
      std::unique_ptr<art::Assns<gar::rec::Hit,gar::rec::TPCCluster> > hitClusterAssns(new ::art::Assns<gar::rec::Hit,gar::rec::TPCCluster>);

      auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e);
      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, hitHandle.id());

      size_t nhits = hits.size();
      std::vector<float> hptmp(nhits);
      std::vector<int> hsi(nhits);
      std::vector<int> used(nhits,0);

      // sort hits by position in x

      for (size_t ihit=0;ihit<nhits;++ihit)
	{
	  hptmp[ihit] = hits.at(ihit).Position()[0];
	}
      TMath::Sort( (int) nhits, hptmp.data(), hsi.data() );

      // loop through hits in x, clustering as we go, marking hits as used.

      for (size_t ihitx=0; ihitx< nhits; ++ihitx)
	{
	  size_t ihit = hsi[ihitx];  // unwound index to hit array
	  if (used[ihit]) continue;

	  const float *xyz = hits.at(ihit).Position();

	  // start a cluster with just this hit
          std::vector<size_t> hitsinclus;  // keep track of hit indices in this cluster so we can make associations
	  hitsinclus.push_back(ihit);
	  used[ihit] = 1;

	  double cpos[3] = {xyz[0], xyz[1], xyz[2]};
	  double csig = hits.at(ihit).Signal();
	  double cstime = hits.at(ihit).StartTime();
	  double cetime = hits.at(ihit).EndTime();
	  double crms = hits.at(ihit).RMS();
	  double ctime = hits.at(ihit).Time();

	  double xyzlow[3] =
	    {
	      xyz[0] - fHitClusterDx,
	      xyz[1] - fHitClusterDyDz,
	      xyz[2] - fHitClusterDyDz
	    };
	  if (xyzlow[0]<0 && xyz[0] >= 0) xyzlow[0] = 0;  // don't cluster across the cathode

	  double xyzhigh[3] =
	    {
	      xyz[0] + fHitClusterDx,
	      xyz[1] + fHitClusterDyDz,
	      xyz[2] + fHitClusterDyDz
	    };
	  if (xyzhigh[0]>0 && xyz[0] <= 0) xyzhigh[0] = 0;  // don't cluster across the cathode

	  for (size_t ix = ihitx+1; ix<nhits; ++ix) // look for candidate hits to cluster in with this one
	    {
	      size_t ihc = hsi[ix];  // candidate hit to add to the cluster if it's in range
	      const float *xyz2 = hits.at(ihc).Position();
	      if (xyz2[0] > xyzhigh[0] || xyz2[0] < xyzlow[0]) break;
	      if (xyz2[1] < xyzhigh[1] && xyz2[1] > xyzlow[1] && xyz2[2] < xyzhigh[2] && xyz2[2] > xyzlow[2] && (used[ihc] == 0))
		{
		  // add hit to cluster
		  hitsinclus.push_back(ihc);
		  used[ihc] = 1;
		  double signal = hits.at(ihit).Signal();
		  double totsig = csig + signal;
		  if (totsig > 0)
		    {
		      for (size_t idim=0; idim<3; ++idim)
			{
			  cpos[idim] = ( cpos[idim]*csig + xyz2[idim]*signal) / totsig;
			}
		      cstime = TMath::Min(cstime, (double) hits.at(ihit).StartTime());
		      cetime = TMath::Max(cetime, (double) hits.at(ihit).EndTime());

		      double htime = hits.at(ihit).Time();
		      double hrms = hits.at(ihit).RMS();
		      crms = TMath::Sqrt( (csig*crms*crms + signal*hrms*hrms)/totsig  +
					  (csig*signal)*TMath::Sq(htime-ctime)/TMath::Sq(totsig));
		      ctime = (ctime*csig + htime*signal) / totsig;
		      csig = totsig;
		    }
		}
	    }

	  float fcpos[3] = {0,0,0};
	  for (int i=0;i<3;++i)
	    {
	      fcpos[i] = cpos[i];
	    }
	  TPCClusterCol->emplace_back(csig,
				      fcpos,
				      cstime,
				      cetime,
				      ctime,
				      crms);
	  auto const TPCClusterPointer = TPCClusterPtrMaker(TPCClusterCol->size()-1);
	  for (size_t i=0; i<hitsinclus.size(); ++i)
	    {
	      auto const HitPointer = hitPtrMaker(hitsinclus.at(i));
	      hitClusterAssns->addSingle(HitPointer,TPCClusterPointer);
	    }

	}

      e.put(std::move(TPCClusterCol));
      e.put(std::move(hitClusterAssns));

    }

    DEFINE_ART_MODULE(TPCHitCluster)

  } // namespace rec
} // namespace gar
