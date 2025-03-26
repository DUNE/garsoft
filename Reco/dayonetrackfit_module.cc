////////////////////////////////////////////////////////////////////////
// Class:       dayonetrackfit
// Plugin Type: producer (art v3_00_00)
// File:        dayonetrackfit_module.cc
//
// copied from tpctrackfit2_module.cc and modified for the day one tracker
// need to step in Z and not X and compute track parameters nonlinearly
// second try, non-kalman: starting with the patrec track parameters
// and minimize chisquared swimming the particle from one station to the next
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <set>

// ROOT includes

#include "TMath.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "TVectorF.h"

// GArSoft Includes
#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/GArMagneticField.h"
#include "Geometry/GeometryGAr.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "RecoAlg/TrackPropagator.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "ReconstructionDataProducts/VecHit.h"

#include "Geant4/G4ThreeVector.hh"

#include "nug4/MagneticFieldServices/MagneticFieldService.h"

namespace gar {
  namespace rec {

    class dayonetrackfit : public art::EDProducer {
    public:
      explicit dayonetrackfit(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      dayonetrackfit(dayonetrackfit const&) = delete;
      dayonetrackfit(dayonetrackfit&&) = delete;
      dayonetrackfit& operator=(dayonetrackfit const&) = delete;
      dayonetrackfit& operator=(dayonetrackfit&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:
      // Declare member data here.

      std::string fPatRecLabel; ///< input patrec tracks and associations
      int fPrintLevel; ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      float fTPCClusterResolXY; ///< pad size in cm in YZ to determine step size
      float
        fTPCClusterResolZ; ///< drift direction contribution to determine step size (resolution of a TPCCluster)
      int fDumpTracks; ///< 0: do not print out tracks, 1: print out tracks
      //      float fRoadYZinFit;                  ///< cut in cm for dropping TPCClusters from tracks in fit
      //      float fMinIonizGapCut;               ///< Don't compute dEdx for this dx or larger

      int DayOneFitBothWays(const gar::rec::Track& patrectrack,
                            std::vector<gar::rec::TPCCluster>& TPCClusters,
                            TrackPar& trackpar,
                            TrackIoniz& trackions,
                            TrackTrajectory& tracktraj);

      int DayOneFit(std::vector<TPCCluster>& TPCClusters,
                    std::vector<int>& TPCClusterlist, // sort ordered list
                    std::vector<float>& trackparbeg,
                    float& chisquared,
                    float& length,
                    float* covmat, // 5x5 covariance matrix
                    std::vector<std::pair<float, float>>& dSigdXs,
                    std::vector<TVector3>& trajpts);

      void updatepar(TVectorF& parvec,
                     float zcur,
                     float xh,
                     float yh,
                     float zh,
                     TVectorF& predstep,
                     float& dlength);
    };

    dayonetrackfit::dayonetrackfit(fhicl::ParameterSet const& p) : EDProducer{p}
    {
      // Call appropriate produces<>() functions here.
      // Call appropriate consumes<>() for any products to be retrieved by this module.

      fPatRecLabel = p.get<std::string>("PatRecLabel", "patrec");
      fPrintLevel = p.get<int>("PrintLevel", 0);
      fDumpTracks = p.get<int>("DumpTracks", 0);
      fTPCClusterResolXY =
        p.get<float>("TPCClusterResolXY", 0.3); // TODO -- think about what this value is
      fTPCClusterResolZ = p.get<float>("TPCClusterResolZ", 0.3); // this is probably much better

      art::InputTag patrecTag(fPatRecLabel);
      consumes<std::vector<gar::rec::Track>>(patrecTag);
      consumes<art::Assns<gar::rec::TPCCluster, gar::rec::Track>>(patrecTag);

      produces<std::vector<gar::rec::Track>>();
      produces<art::Assns<gar::rec::TPCCluster, gar::rec::Track>>();
      produces<std::vector<gar::rec::TrackIoniz>>();
      produces<art::Assns<rec::TrackIoniz, rec::Track>>();
      produces<std::vector<gar::rec::TrackTrajectory>>();
      produces<art::Assns<rec::TrackTrajectory, rec::Track>>();
    }

    void dayonetrackfit::produce(art::Event& e)
    {
      // output collections

      std::unique_ptr<std::vector<gar::rec::Track>> trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr<art::Assns<gar::rec::TPCCluster, gar::rec::Track>> TPCClusterTrkAssns(
        new ::art::Assns<gar::rec::TPCCluster, gar::rec::Track>);
      std::unique_ptr<std::vector<rec::TrackIoniz>> ionCol(new std::vector<rec::TrackIoniz>);
      std::unique_ptr<art::Assns<rec::TrackIoniz, rec::Track>> ionTrkAssns(
        new ::art::Assns<rec::TrackIoniz, rec::Track>);
      std::unique_ptr<std::vector<rec::TrackTrajectory>> trajCol(
        new std::vector<rec::TrackTrajectory>);
      std::unique_ptr<art::Assns<rec::TrackTrajectory, rec::Track>> trajTrkAssns(
        new ::art::Assns<rec::TrackTrajectory, rec::Track>);

      // inputs

      auto patrecTrackHandle = e.getValidHandle<std::vector<gar::rec::Track>>(fPatRecLabel);
      auto const& patrecTracks = *patrecTrackHandle;

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const ionizPtrMaker = art::PtrMaker<rec::TrackIoniz>(e);
      auto const trajPtrMaker = art::PtrMaker<rec::TrackTrajectory>(e);
      //auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());

      auto const* magFieldService = gar::providerFrom<mag::MagneticFieldService>();
      G4ThreeVector zerovec(0, 0, 0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::ServiceHandle<geo::GeometryGAr> geo;
      //double xtpccent = geo->GetOriginX();
      //double ytpccent = geo->GetOriginY();
      //double ztpccent = geo->GetOriginZ();
      //TVector3 tpccent(xtpccent,ytpccent,ztpccent);
      TVector3 xhat(1, 0, 0);

      const art::FindManyP<gar::rec::TPCCluster> TPCClustersFromPatRecTracks(
        patrecTrackHandle, e, fPatRecLabel);

      for (size_t itrack = 0; itrack < patrecTracks.size(); ++itrack) {
        std::vector<gar::rec::TPCCluster> TPCClusters;
        for (size_t iTPCCluster = 0; iTPCCluster < TPCClustersFromPatRecTracks.at(itrack).size();
             ++iTPCCluster) {
          TPCClusters.push_back(*TPCClustersFromPatRecTracks.at(itrack).at(
            iTPCCluster)); // make our own local copy of TPCClusters.  Maybe we can skip this?
        }
        TrackPar trackparams;
        TrackIoniz trackions;
        TrackTrajectory tracktraj;
        if (DayOneFitBothWays(
              patrecTracks.at(itrack), TPCClusters, trackparams, trackions, tracktraj) == 0) {
          trkCol->push_back(trackparams.CreateTrack());
          ionCol->push_back(trackions);
          trajCol->push_back(tracktraj);
          auto const trackpointer = trackPtrMaker(trkCol->size() - 1);
          auto const ionizpointer = ionizPtrMaker(ionCol->size() - 1);
          auto const trajpointer = trajPtrMaker(trajCol->size() - 1);
          for (size_t iTPCCluster = 0; iTPCCluster < TPCClusters.size(); ++iTPCCluster) {
            TPCClusterTrkAssns->addSingle(TPCClustersFromPatRecTracks.at(itrack).at(iTPCCluster),
                                          trackpointer);
          }
          ionTrkAssns->addSingle(ionizpointer, trackpointer);
          trajTrkAssns->addSingle(trajpointer, trackpointer);
        }
      }

      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterTrkAssns));
      e.put(std::move(ionCol));
      e.put(std::move(ionTrkAssns));
      e.put(std::move(trajCol));
      e.put(std::move(trajTrkAssns));
    }

    int dayonetrackfit::DayOneFitBothWays(const gar::rec::Track& patrectrack,
                                          std::vector<gar::rec::TPCCluster>& TPCClusters,
                                          TrackPar& trackpar,
                                          TrackIoniz& trackions,
                                          TrackTrajectory& tracktraj)

    {
      // track parameters:  x is the independent variable in TPC track fits, but just report it here.
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end

      size_t nclus = TPCClusters.size();
      std::vector<float> TPCClusterz;
      for (size_t i = 0; i < nclus; ++i) {
        TPCClusterz.push_back(TPCClusters[i].Position()[2]);
      }
      std::vector<int> hsi(nclus); // "hit sort index" -- really TPC clusters
      TMath::Sort((int)nclus, TPCClusterz.data(), hsi.data());
      std::vector<int> hsib(nclus);
      for (size_t i = 0; i < nclus; ++i) {
        hsib[nclus - i - 1] = hsi[i];
      }

      std::vector<float> tparbeg(6); // fill with initial track parameters from patrec track
      for (size_t i = 0; i < 5; ++i) {
        tparbeg[i] = patrectrack.TrackParBeg()[i];
      }
      tparbeg[5] = patrectrack.Vertex()[0];
      float covmatbeg[25];
      patrectrack.CovMatBegSymmetric(covmatbeg);
      float chisqforwards = 0;
      float lengthforwards = 0;
      std::vector<std::pair<float, float>> dSigdXs_FWD;
      std::vector<TVector3> trajpts_FWD;

      int retcode = DayOneFit(TPCClusters,
                              hsi,
                              tparbeg,
                              chisqforwards,
                              lengthforwards,
                              covmatbeg,
                              dSigdXs_FWD,
                              trajpts_FWD);
      if (retcode != 0) return 1;

      std::vector<float> tparend(6); // fill with end track parameters from patrec track
      for (size_t i = 0; i < 5; ++i) {
        tparend[i] = patrectrack.TrackParEnd()[i];
      }
      tparend[5] = patrectrack.End()[0];
      float covmatend[25];
      patrectrack.CovMatEndSymmetric(covmatend);
      float chisqbackwards = 0;
      float lengthbackwards = 0;
      std::vector<std::pair<float, float>> dSigdXs_BAK;
      std::vector<TVector3> trajpts_BAK;

      retcode = DayOneFit(TPCClusters,
                          hsib,
                          tparend,
                          chisqbackwards,
                          lengthbackwards,
                          covmatend,
                          dSigdXs_BAK,
                          trajpts_BAK);
      if (retcode != 0) return 1;

      size_t nTPCClusters = TPCClusters.size();
      trackpar.setNTPCClusters(nTPCClusters);
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

      trackions.setData(dSigdXs_FWD, dSigdXs_BAK);
      tracktraj.setData(trajpts_FWD, trajpts_BAK);

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    // track parameters visible to caller

    // 0: y
    // 1: z
    // 2: curvature
    // 3: phi
    // 4: lambda
    // 5: x

    int dayonetrackfit::DayOneFit(std::vector<TPCCluster>& TPCClusters,
                                  std::vector<int>& TPCClusterlist, // sort ordered list
                                  std::vector<float>& trackparbeg,
                                  float& chisquared,
                                  float& length,
                                  float* /* covmat */, // 5x5 covariance matrix
                                  std::vector<std::pair<float, float>>& /* dSigdXs */,
                                  std::vector<TVector3>& trajpts)
    {

      // set some default values in case we return early

      size_t nclus = TPCClusterlist.size();
      chisquared = 0;
      length = 0;

      // first try, do nothing, but we should calculate a chisquared

      for (size_t i = 0; i < nclus; ++i) {
        const float* pos = TPCClusters.at(TPCClusterlist.at(i)).Position();
        TVector3 tvec(pos[0], pos[1], pos[2]);
        trajpts.push_back(tvec);
        float dist = 0;
        float tvpos[3] = {trackparbeg[5], trackparbeg[0], trackparbeg[1]};
        int retcode = util::TrackPropagator::DistXYZ(trackparbeg.data(), tvpos, pos, dist);
        if (retcode == 0) {
          chisquared +=
            TMath::Sq(dist / fTPCClusterResolXY); // todo -- break this up into X,Y vs. Z pieces.
        }
      }

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------

    // extrapolate a helix using the current helix parameters at parvec to a new hit position, given by z.
    // ouptut is predstep

    // 0: y
    // 1: x
    // 2: curvature
    // 3: phi
    // 4: lambda

    void dayonetrackfit::updatepar(TVectorF& parvec,
                                   float zcur,
                                   float xh,
                                   float /* yh */,
                                   float zh,
                                   TVectorF& predstep,
                                   float& dlength)
    {

      predstep = parvec;

      // radius and center of the circle

      float r = 1E4;
      if (parvec[2] != 0) { r = 1 / parvec[2]; }
      float yc = parvec[0] + r * TMath::Cos(parvec[3]);
      float zc = zcur - r * TMath::Sin(parvec[3]);

      // two solutions of arcsin, plus any number of 2pi.  We could be on
      // the boundary of a wraparound.

      float phip1 = TMath::ASin((zh - zc) / r);
      float phip2 = TMath::Pi() - phip1;

      // find the solution that's closest to the one suggested by x

      float tanlambda = TMath::Tan(parvec[4]);
      if (tanlambda == 0) tanlambda = 1E-4;
      float phip3 = parvec[3] + (xh - parvec[1]) / (r * tanlambda);

      float dphip1n = (phip1 - phip3) / TMath::TwoPi();
      float niphi1 = nearbyint(dphip1n);
      float dphip1a = TMath::Abs(dphip1n - niphi1);

      float dphip2n = (phip2 - phip3) / TMath::TwoPi();
      float niphi2 = nearbyint(dphip2n);
      float dphip2a = TMath::Abs(dphip2n - niphi2);

      float phiex = phip1 + niphi1 * TMath::TwoPi();
      if (dphip1a > dphip2a) { phiex = phip2 + niphi2 * TMath::TwoPi(); }

      // update predstep to values at this z.  Leave curvature and lambda
      // unchanged

      predstep[0] = yc - r * TMath::Cos(phiex);                      // y
      predstep[1] = parvec[1] + r * tanlambda * (phiex - parvec[3]); // x
      predstep[3] = phiex;

      dlength = TMath::Abs(r * (phiex - parvec[3]) / TMath::Cos(parvec[4]));
    }

    DEFINE_ART_MODULE(dayonetrackfit)

  } // namespace rec
} // namespace gar
