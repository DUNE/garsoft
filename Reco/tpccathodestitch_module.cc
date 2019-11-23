////////////////////////////////////////////////////////////////////////
// Class:       tpccathodestitch
// Plugin Type: producer (art v3_00_00)
// File:        tpccathodestitch_module.cc
//
// Generated at Thu Nov 21 15:00:00 2019 by Thomas Junk modified from tpctrackfit2_module.cc
//  Finds pairs of tracks on either side of the cathode plane and merges them.
//  Inputs -- tracks fit by the track fitter, and associations with TPCClusters
//   Outputs -- a completely new set of tracks, trackioniz, and associations of tracks with TPCClusters and trackioniz.  
//   No change for tracks that have not been stitched, just copies.  Tracks that have been stitched correspond to fewer tracks on output
//    also fills in a time for the stitched track

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

// ROOT includes

#include "TMath.h"
#include "TVector3.h"

// GArSoft Includes
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "Geometry/Geometry.h"

#include "Geant4/G4ThreeVector.hh"

#include "nutools/MagneticField/MagneticField.h"

namespace gar {
  namespace rec {

    class tpccathodestitch : public art::EDProducer {
    public:
      explicit tpccathodestitch(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpccathodestitch(tpccathodestitch const&) = delete;
      tpccathodestitch(tpccathodestitch&&) = delete;
      tpccathodestitch& operator=(tpccathodestitch const&) = delete;
      tpccathodestitch& operator=(tpccathodestitch&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      // Declare member data here.

      std::string fInputTrackLabel;     ///< input tracks and associations
      int fPrintLevel;              ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      float fDistCut;               ///< cut in distance between best-matching points
      float fCTCut;                 ///< cut on cosine of angle matching
      float fMaxDX;                 ///< cut on maximum translation of dX on one side.
      float fMinDX;                 ///< cut on minimum translation of dX on one side.

      // returns true if trackA and trackB look like they were split across the cathode and need merging
      // fwdflag variables say whether the beginning endpoint is the one on the stitched end (1) or the endpoint
      // is the one on the stitched end (2).  Both flags are set to zero if the track is not stitchable
      // deltaX is the shift in X needed to be applied to track A, and -deltaX is to be applied to track B.

      bool cathodematch(const gar::rec::Track &trackA, const gar::rec::Track &trackB, int &fwdflagA, int &fwdflagB, float &deltaX);

    };


    tpccathodestitch::tpccathodestitch(fhicl::ParameterSet const& p) : EDProducer{p}  
    {
      fInputTrackLabel   = p.get<std::string>("InputTrackLabel","track");
      fPrintLevel        = p.get<int>("PrintLevel",0);
      fDistCut           = p.get<float>("DistCut",3);
      fCTCut             = p.get<float>("CTCut",0.99);
      fMaxDX             = p.get<float>("MaxDX",50.0);
      fMinDX             = p.get<float>("MinDX",-50.0);

      art::InputTag inputTrackTag(fInputTrackLabel);
      consumes< std::vector<gar::rec::Track> >(inputTrackTag);
      consumes< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >(inputTrackTag);
      consumes<std::vector<gar::rec::TrackIoniz>>(inputTrackTag);
      consumes<art::Assns<rec::TrackIoniz, rec::Track>>(inputTrackTag);

      // probably don't need the vector hits at this point if we have the TPCClusters
      //consumes< std::vector<gar::rec::VecHit> >(patrecTag);
      //consumes< art::Assns<gar::rec::VecHit, gar::rec::Track> >(patrecTag);

      produces< std::vector<gar::rec::Track> >();
      produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();
      produces<std::vector<gar::rec::TrackIoniz>>();
      produces<art::Assns<rec::TrackIoniz, rec::Track>>();
    }



    void tpccathodestitch::produce(art::Event& e)
    {

      // get the distance corresponding to one ADC tick so we can report the time in ticks

      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      auto clockService = gar::providerFrom<detinfo::DetectorClocksService>();
      float DriftVelocity = detProp->DriftVelocity(detProp->Efield(),detProp->Temperature());       // in cm per microsecond
      float distonetick = DriftVelocity * (clockService->TPCTick2Time(1) - clockService->TPCTick2Time(0)) ;

      // output collections

      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);
      std::unique_ptr< std::vector<rec::TrackIoniz> > ionCol(new std::vector<rec::TrackIoniz>);
      std::unique_ptr< art::Assns<rec::TrackIoniz,rec::Track> > ionTrkAssns(new ::art::Assns<rec::TrackIoniz,rec::Track>);

      // inputs

      auto inputTrackHandle = e.getValidHandle< std::vector<gar::rec::Track> >(fInputTrackLabel);
      auto const& inputTracks = *inputTrackHandle;

      const art::FindManyP<gar::rec::TPCCluster> TPCClustersFromInputTracks(inputTrackHandle,e,fInputTrackLabel);
      const art::FindManyP<gar::rec::TrackIoniz> TrackIonizFromInputTracks(inputTrackHandle,e,fInputTrackLabel);

      // we get the input ionization data products from the associations
      //auto inputIonHandle = e.getValidHandle< std::vector<gar::rec::TrackIoniz> >(fInputTrackLabel);
      //auto const& inputIoniz = *inputIonHandle;

      // for making output associations

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const ionizPtrMaker = art::PtrMaker<rec::TrackIoniz>(e);
      //auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());

      // vector of merged flags contains a list of which pairs of tracks are merged
      // A more general track stitcher may want to stitch kinked tracks together and thus be able to stitch
      // more than a pair of tracks, but this module so far just stitches across the cathode.

      size_t nti = inputTracks.size();
      std::vector<int> mergedflag(nti,-1); // save indexes here; negative means unmerged
      std::vector<int> fwdflag(nti,0);     // 1 for forwards, 2 for backwards
      std::vector<float> dx(nti,0);       // shift in X needed

      for (size_t itrack = 0; itrack < inputTracks.size(); ++itrack)
        {
	  if (mergedflag.at(itrack) >= 0) continue;  // this track is already part of an earlier merge
	  for (size_t jtrack = itrack+1; jtrack < inputTracks.size(); ++jtrack)
	    {
	      if (mergedflag.at(jtrack)>=0) continue;  // already merged, don't merge another one
	      if (cathodematch(inputTracks.at(itrack),inputTracks.at(jtrack),fwdflag.at(itrack),fwdflag.at(jtrack),dx.at(itrack)))
		{
		  std::cout << "found a merge.  dx= " << dx.at(itrack) << " flip flags: " << fwdflag.at(itrack) << " " << fwdflag.at(jtrack) <<  std::endl;
		  mergedflag.at(itrack) = jtrack;
		  mergedflag.at(jtrack) = itrack;
		  dx.at(jtrack) = -dx.at(itrack);
		  break;
		}
	    }
	  if (mergedflag.at(itrack) < 0)   // not merged, just copy the input track, trkioniz, and associations to the output
	    {
	      trkCol->push_back(inputTracks.at(itrack));

              auto const trackpointer = trackPtrMaker(trkCol->size()-1);
	      for (size_t iTPCCluster=0; iTPCCluster<TPCClustersFromInputTracks.at(itrack).size(); ++iTPCCluster)
		{
		  TPCClusterTrkAssns->addSingle(TPCClustersFromInputTracks.at(itrack).at(iTPCCluster),trackpointer);
		}
	      // this loop may be unnecessary as there is only one TrackIoniz per Track, but just in case
	      for (size_t iTrkIoniz=0; iTrkIoniz<TrackIonizFromInputTracks.at(itrack).size(); ++iTrkIoniz)
		{
		  ionCol->push_back(*TrackIonizFromInputTracks.at(itrack).at(iTrkIoniz));
                  auto const ionizpointer = ionizPtrMaker(ionCol->size()-1);
                  ionTrkAssns->addSingle(ionizpointer, trackpointer);
		}
	    }
	  else  // merge them -- adjust track parameters and put in the time
	    {
	      // assume directions are set so that itrack is the "beginning" and jtrack is the "end"
	      // flip the tracks around according to the fwdflag returned by the cathode match method

	      TrackPar tpi(inputTracks.at(itrack), fwdflag.at(itrack)==1);
	      int jtrack = mergedflag.at(itrack);
	      TrackPar tpj(inputTracks.at(jtrack), fwdflag.at(jtrack)==1);

	      tpi.setXBeg( tpi.getXBeg() + dx.at(itrack) );
	      tpi.setXEnd( tpi.getXEnd() + dx.at(itrack) );
	      tpj.setXBeg( tpi.getXBeg() + dx.at(jtrack) );
	      tpj.setXEnd( tpi.getXEnd() + dx.at(jtrack) );

	      // some checking due to the unsigned nature of the timestmap.  Assume in ticks.
	      ULong64_t ts = tpi.getTime();
	      int deltat = dx.at(itrack)/distonetick;
	      if ( (int) ts + deltat >= 0)
		{
		  ts += deltat;
		}

	      TrackPar tpm(tpi.getLengthForwards() + tpj.getLengthForwards(),
			   tpi.getLengthBackwards() + tpj.getLengthBackwards(),
			   tpi.getNTPCClusters() + tpj.getNTPCClusters(),
			   tpi.getXBeg(),
			   tpi.getTrackParametersBegin(),
			   tpi.getCovMatBeg(),
			   tpi.getChisqForwards() + tpj.getChisqForwards(),
			   tpj.getXEnd(),
			   tpj.getTrackParametersEnd(),
			   tpj.getCovMatEnd(),
			   tpi.getChisqBackwards() + tpj.getChisqBackwards(),
			   ts);
	      trkCol->push_back(tpm.CreateTrack());
              auto const trackpointer = trackPtrMaker(trkCol->size()-1);

	      for (size_t iTPCCluster=0; iTPCCluster<TPCClustersFromInputTracks.at(itrack).size(); ++iTPCCluster)
		{
		  TPCClusterTrkAssns->addSingle(TPCClustersFromInputTracks.at(itrack).at(iTPCCluster),trackpointer);
		}
	      for (size_t iTPCCluster=0; iTPCCluster<TPCClustersFromInputTracks.at(jtrack).size(); ++iTPCCluster)
		{
		  TPCClusterTrkAssns->addSingle(TPCClustersFromInputTracks.at(jtrack).at(iTPCCluster),trackpointer);
		}

	      // merge the track ionization data -- assume just one TrackIoniz per track, with forward and backward vectors
	      // reverse them if need be.

	      const auto itif = TrackIonizFromInputTracks.at(itrack).at(0)->getFWD_dSigdXs();
	      const auto itib = TrackIonizFromInputTracks.at(itrack).at(0)->getBAK_dSigdXs();
	      const auto jtif = TrackIonizFromInputTracks.at(jtrack).at(0)->getFWD_dSigdXs();
	      const auto jtib = TrackIonizFromInputTracks.at(jtrack).at(0)->getBAK_dSigdXs();

	      auto mtif = itif;
	      auto mtib = itib;
	      if (fwdflag.at(itrack) != 1)      //  flip itrack's ionization vector if need be
		{
		  mtif = itib;
		  mtib = itif;
		}
	      auto jmtif = jtif;
	      auto jmtib = jtib;
	      if (fwdflag.at(jtrack) != 1)
		{
		  jmtif = jtib;
		  jmtib = jtif;
		}

	      for (size_t i=0; i<jmtif.size(); ++i)
		{
		  mtif.push_back(jmtif.at(i));
		}
	      for (size_t i=0; i<jmtib.size(); ++i)
		{
		  mtib.push_back(jmtib.at(i));
		}
	      TrackIoniz tim;
	      tim.setData(mtif,mtib);
	      ionCol->push_back(tim);
              auto const ionizpointer = ionizPtrMaker(ionCol->size()-1);
              ionTrkAssns->addSingle(ionizpointer, trackpointer);
	    }
	}

      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterTrkAssns));
      e.put(std::move(ionCol));
      e.put(std::move(ionTrkAssns));
    }

    // returns true if trackA and trackB look like they were split across the cathode and need merging
    // fwdflag variables say whether the beginning endpoint is the one on the stitched end (1) or the endpoint
    // is the one on the stitched end (2).  Both flags are set to zero if the track is not stitchable
    // deltaX is the shift in X needed to be applied to track A, and -deltaX is to be applied to track B.

    bool tpccathodestitch::cathodematch(const gar::rec::Track &atrack, const gar::rec::Track &btrack, int &afw, int &bfw, float &deltaX)
    {
      afw = 0;  
      bfw = 0;
      bool matchable = false;

      // may not need this yet -- but magnetic field will be a function of position someday
      // art::ServiceHandle<mag::MagneticField> magFieldService;
      // G4ThreeVector zerovec(0,0,0);
      // G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::ServiceHandle<geo::Geometry> geo;
      double xtpccent = geo->TPCXCent();
      TVector3 xhat(1,0,0);
      TVector3 xcentv = xhat*xtpccent;

      std::vector<TVector3> apt;  // endpoint
      std::vector<TVector3> bpt;
      std::vector<TVector3> adir; // direction
      std::vector<TVector3> bdir;

      apt.emplace_back(atrack.Vertex());
      apt.emplace_back(atrack.End());
      adir.emplace_back(atrack.VtxDir());
      adir.emplace_back(atrack.EndDir());
      bpt.emplace_back(btrack.Vertex());
      bpt.emplace_back(btrack.End());
      bdir.emplace_back(btrack.VtxDir());
      bdir.emplace_back(btrack.EndDir());

      // center the postions on the center of the detector

      for (size_t i=0; i<apt.size(); ++i) apt.at(i) -= xcentv;
      for (size_t i=0; i<bpt.size(); ++i) bpt.at(i) -= xcentv;

      std::vector<int> fwflag = {2, 1};

      for (size_t ia=0; ia<apt.size(); ++ia)
	{
	  for (size_t ib=0; ib<bpt.size(); ++ib)
	    {
	      // directions should point along track pieces to stitch, so they should be back to back
	      TVector3 v=adir.at(ia);
	      TVector3 u=bdir.at(ib);

	      if (v.Dot(u) > -fCTCut) continue;

	      // simple comparisons of vertex positions.  Ignores the fact that we may be missing
	      // part of one of the tracks near the endpoint due to the gaps or patrec failure
	      //if ( ((apt.at(ia)-bpt.at(ib)).Cross(xhat)).Mag() > dYZCut) continue;
	      //if (TMath::Abs(apt.at(ia).X() + bpt.at(ib).X()) > fdXCut) continue;

	      // solved for X shift called tau

	      float tdenom = 1.0 - TMath::Sq(xhat.Dot(v));
	      if (tdenom == 0) continue;
	      TVector3 d = apt.at(ia) - bpt.at(ib);
	      float tau = (0.5/tdenom) * d.Dot( (v.Dot(xhat))*v - xhat );

	      if (tau > fMaxDX) continue;
	      if (tau < fMinDX) continue;

	      // solve for displacement along adir.

	      float gamma = -v.Dot(2*tau*xhat + d);
	      TVector3 chivec = d + gamma*v + 2.0*tau*xhat;
	      if (chivec.Mag() > fDistCut) continue;

	      // we have a match
	      afw = fwflag.at(ia);
	      bfw = fwflag.at(1-ib);
	      deltaX = tau;
	      matchable = true;

	      break;   // just find the first match -- possible to match the same tracks with different endpoints?
	    }
	  if (matchable) break;
	}
      return matchable;
    }



    DEFINE_ART_MODULE(tpccathodestitch)

  } // namespace rec
} // namespace gar
