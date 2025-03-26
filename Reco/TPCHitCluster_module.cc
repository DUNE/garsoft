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
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "Geometry/GeometryGAr.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "TMath.h"

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
      std::string fHitLabel; ///< label of module to get the hits from
      int fClusterHits;      ///< hit clustering algorithm number
      float fHitClusterDx;   ///< range in cm to look for hits to cluster in x
      float fHitClusterDyDz; ///< range in cm to look for hits to cluster in y and z

      float fHitClusterDr_I;  ///< range in R to look for hits to cluster for IROC, in cm
      float fHitClusterDfi_I; ///< range in R*dPhi to look for hits to cluster for IROC, in cm
      float fHitClusterDrIO;  ///< range in R to look for hits to cluster for IOROC, in cm
      float fHitClusterDfiIO; ///< range in R*dPhi to look for hits to cluster for IROC, in cm
      float fHitClusterDrOO;  ///< range in R to look for hits to cluster for OOROC, in cm
      float fHitClusterDfiOO; ///< range in R*dPhi to look for hits to cluster for IROC, in cm

      int fPrintLevel; ///< debug printout:  0: none, 1: selected, 2: all
    };

    //==============================================================================
    //==============================================================================
    //==============================================================================
    TPCHitCluster::TPCHitCluster(fhicl::ParameterSet const& p) : EDProducer{p}
    {
      fHitLabel = p.get<std::string>("HitLabel", "hit");
      fPrintLevel = p.get<int>("PrintLevel", 0);
      fClusterHits = p.get<int>("ClusterHits", 1);
      fHitClusterDx = p.get<float>("HitClusterDx", 1.00);
      fHitClusterDyDz = p.get<float>("HitClusterDyDz", 1.20);
      fHitClusterDr_I = p.get<float>("HitClusterDr_I", 1.50);
      fHitClusterDfi_I = p.get<float>("HitClusterDfi_I", 1.06);
      fHitClusterDrIO = p.get<float>("HitClusterDrIO", 2.00);
      fHitClusterDfiIO = p.get<float>("HitClusterDfiIO", 1.55);
      fHitClusterDrOO = p.get<float>("HitClusterDrOO", 3.00);
      fHitClusterDfiOO = p.get<float>("HitClusterDfiO)", 1.73);

      art::InputTag hitTag(fHitLabel);
      consumes<std::vector<gar::rec::Hit>>(hitTag);
      produces<std::vector<gar::rec::TPCCluster>>();
      produces<art::Assns<gar::rec::Hit, gar::rec::TPCCluster>>();
    }

    //==============================================================================
    //==============================================================================
    //==============================================================================
    void TPCHitCluster::produce(art::Event& e)
    {
      art::ServiceHandle<geo::GeometryGAr> euclid;
      float xtpccent = euclid->TPCXCent();
      float ytpccent = euclid->TPCYCent();
      float ztpccent = euclid->TPCZCent();

      // input: hits

      auto hitHandle = e.getValidHandle<std::vector<gar::rec::Hit>>(fHitLabel);
      auto const& hits = *hitHandle;

      // output: TPCCluster and associations with hits

      std::unique_ptr<std::vector<gar::rec::TPCCluster>> TPCClusterCol(
        new std::vector<gar::rec::TPCCluster>);
      std::unique_ptr<art::Assns<gar::rec::Hit, gar::rec::TPCCluster>> hitClusterAssns(
        new ::art::Assns<gar::rec::Hit, gar::rec::TPCCluster>);

      auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e);
      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e, hitHandle.id());

      size_t nhits = hits.size();
      std::vector<float> hptmp(nhits);
      std::vector<int> hsi(nhits);
      std::vector<int> used(nhits, 0);

      // sort hits by position in x

      for (size_t ihit = 0; ihit < nhits; ++ihit) {
        hptmp[ihit] = hits.at(ihit).Position()[0];
      }
      TMath::Sort((int)nhits, hptmp.data(), hsi.data());

      // loop through hits in x, clustering as we go, marking hits as used.

      for (size_t ihitx = 0; ihitx < nhits; ++ihitx) {
        size_t ihit = hsi[ihitx]; // unwound index to hit array
        if (used[ihit]) continue;

        const float* xyz_fromhit = hits.at(ihit).Position();
        float xyz[3] = {
          xyz_fromhit[0] - xtpccent, xyz_fromhit[1] - ytpccent, xyz_fromhit[2] - ztpccent};

        // start a cluster with just this hit
        std::vector<size_t>
          hitsinclus; // keep track of hit indices in this cluster so we can make associations
        hitsinclus.push_back(ihit);
        used[ihit] = 1;

        int side = 0; // detector side of the seed hit.
        auto hitchan = hits.at(ihit).Channel();
        float chanpos[3] = {0, 0, 0};
        euclid->ChannelToPosition(hitchan, chanpos);
        if (chanpos[0] > xtpccent) { side = 1; }
        else {
          side = -1;
        }

        // Clusters can be displaced in time.  We may want to save which end of the
        // chamber the data came from however and not cluster together charge
        // that comes from different sides.

        if (fPrintLevel > 1) {
          std::cout << "Started a new TPC cluster.  Pos= " << xyz[0] << " " << xyz[1] << " "
                    << xyz[2] << std::endl;
          std::cout << "Side: " << side << std::endl;
        }

        // Determine cluster position in endcap for cuts to accumulate more hits
        float rHit = std::hypot(xyz[1], xyz[2]);
        bool In_CROC = rHit < euclid->GetIROCInnerRadius();
        bool In_IROC = euclid->GetIROCInnerRadius() < rHit && rHit < euclid->GetIROCOuterRadius();
        bool InIOROC =
          euclid->GetOROCInnerRadius() < rHit && rHit < euclid->GetOROCPadHeightChangeRadius();
        bool InOOROC =
          euclid->GetOROCPadHeightChangeRadius() < rHit && rHit < euclid->GetOROCOuterRadius();

        for (size_t ix = ihitx + 1; ix < nhits;
             ++ix) // look for candidate hits to cluster in with this one
        {
          size_t ihc = hsi[ix]; // candidate hit to add to the cluster if it's in range
          if (used[ihc]) continue;

          int sidetest = 0;
          auto hitchantest = hits.at(ihc).Channel();
          float chanpostest[3] = {0, 0, 0};
          euclid->ChannelToPosition(hitchantest, chanpostest);
          if (chanpostest[0] > xtpccent) { sidetest = 1; }
          else {
            sidetest = -1;
          }
          if (sidetest != side)
            continue; // don't cluster hits that were detected on opposite TPC endplates.

          const float* xyz2_fromhit = hits.at(ihc).Position();
          float xyz2[3] = {
            xyz2_fromhit[0] - xtpccent, xyz2_fromhit[1] - ytpccent, xyz2_fromhit[2] - ztpccent};

          if (fPrintLevel > 1) {
            std::cout << " Testing a hit: " << xyz2[0] << " " << xyz2[1] << " " << xyz2[2] << " "
                      << sidetest << " " << used[ihc] << " " << ihc << " " << ix << std::endl;
          }

          // Check if hit ix is close to ihit, first in drift direction, then
          // in transverse plane depending on which ROC has ihit
          if (fabs(xyz2[0] - xyz[0]) > fHitClusterDx) break;

          if (In_CROC) {
            if (fabs(xyz2[1] - xyz[1]) > fHitClusterDyDz) continue;
            if (fabs(xyz2[2] - xyz[2]) > fHitClusterDyDz) continue;
          }
          // Project out the radial and rotational components of xyz2 -xyz.
          float RhatY = xyz[1] / std::hypot(xyz[1], xyz[2]);
          float RhatZ = xyz[2] / std::hypot(xyz[1], xyz[2]);
          float clustDist = std::hypot(xyz2[1] - xyz[1], xyz2[2] - xyz[2]);
          float dR = TMath::Abs((xyz2[1] - xyz[1]) * RhatY + (xyz2[2] - xyz[2]) * RhatZ);
          float Rdfi = sqrt(clustDist * clustDist - dR * dR);
          if (In_IROC) {
            if (dR > fHitClusterDr_I) continue;
            if (Rdfi > fHitClusterDfi_I) continue;
          }
          if (InIOROC) {
            if (dR > fHitClusterDrIO) continue;
            if (Rdfi > fHitClusterDfiIO) continue;
          }
          if (InOOROC) {
            if (dR > fHitClusterDrOO) continue;
            if (Rdfi > fHitClusterDfiOO) continue;
          }

          // add hit to cluster
          hitsinclus.push_back(ihc);
          used[ihc] = 1;
          if (fPrintLevel > 1) {
            std::cout << "Added hit with pos: " << xyz2[0] << " " << xyz2[1] << " " << xyz2[2]
                      << std::endl;
          }
        }

        // calculate cluster charge, centroid, min and max times, and covariance

        double cpos[3] = {0, 0, 0}; // charge-weighted cluster position
        double csig = 0;            // cluster signal total
        double cstime = 0;          // cluster start time
        double cetime = 0;          // cluster end time
        double crms = 0;            // cluster RMS in X
        double ctime = 0;           // cluster time centroid
        double cov[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

        for (size_t ix = 0; ix < hitsinclus.size(); ++ix) {
          size_t ihx = hitsinclus.at(ix);
          double hsig = hits.at(ihx).Signal();
          csig += hsig;
          for (size_t idim = 0; idim < 3; ++idim) {
            cpos[idim] += hsig * hits.at(ihx).Position()[idim];
          }
          ctime += hsig * hits.at(ihx).Time();
          if (ix == 0) {
            cstime = hits.at(ihx).StartTime();
            cetime = hits.at(ihx).EndTime();
          }
          else {
            cstime = TMath::Min(cstime, (double)hits.at(ihx).StartTime());
            cetime = TMath::Max(cetime, (double)hits.at(ihx).EndTime());
          }
        }
        if (csig != 0) {
          ctime /= csig;
          for (size_t idim = 0; idim < 3; ++idim) {
            cpos[idim] /= csig;
          }

          // calculate covariance.  Do a numerical sum over hit widths assuming Gaussians.

          for (size_t ix = 0; ix < hitsinclus.size(); ++ix) {
            size_t ihx = hitsinclus.at(ix);
            double cfrac = hits.at(ihx).Signal() / csig;
            const float* hpos = hits.at(ihx).Position();
            for (size_t idim = 0; idim < 3; ++idim) {
              for (size_t jdim = 0; jdim < 3; ++jdim) {
                cov[idim][jdim] += cfrac * (hpos[idim] - cpos[idim]) * (hpos[jdim] - cpos[jdim]);
              }
            }
            cov[0][0] += cfrac * TMath::Sq(hits.at(ihx).RMS());
          }
        }
        crms =
          TMath::Sqrt(cov[0][0]); // to mean what it meant before -- RMS along the drift direction

        float fcpos[3] = {0, 0, 0};
        for (int i = 0; i < 3; ++i) {
          fcpos[i] = cpos[i];
        }
        float fccov[6] = {(float)cov[0][0],
                          (float)cov[1][0],
                          (float)cov[2][0],
                          (float)cov[1][1],
                          (float)cov[1][2],
                          (float)cov[2][2]};

        TPCClusterCol->emplace_back(csig, fcpos, cstime, cetime, ctime, crms, fccov);

        if (fPrintLevel > 0) {
          std::cout << "Made a TPC Cluster pos: " << fcpos[0] << " " << fcpos[1] << " " << fcpos[2]
                    << "  signal: " << csig << "  nHits: " << hitsinclus.size() << std::endl;
        }
        auto const TPCClusterPointer = TPCClusterPtrMaker(TPCClusterCol->size() - 1);
        for (size_t i = 0; i < hitsinclus.size(); ++i) {
          auto const HitPointer = hitPtrMaker(hitsinclus.at(i));
          hitClusterAssns->addSingle(HitPointer, TPCClusterPointer);
        }
      }

      e.put(std::move(TPCClusterCol));
      e.put(std::move(hitClusterAssns));
    }

    DEFINE_ART_MODULE(TPCHitCluster)

  } // namespace rec
} // namespace gar
