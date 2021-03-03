////////////////////////////////////////////////////////////////////////
// Class:       CaloClustering
// Plugin Type: producer (art v2_11_02)
// File:        CaloClustering_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "RecoAlg/KNNClusterAlg.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include <memory>

namespace gar {
    namespace rec {

        namespace alg{
            class KNNClusterFinderAlg;
        }

        class CaloClustering : public art::EDProducer {
        public:
            explicit CaloClustering(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            CaloClustering(CaloClustering const &) = delete;
            CaloClustering(CaloClustering &&) = delete;
            CaloClustering & operator = (CaloClustering const &) = delete;
            CaloClustering & operator = (CaloClustering &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        private:

            // Declare member data here.
            void CollectTracks(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::Track> > &trkVector);
            void CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

            std::string fTrackLabel;  ///< label to find the reco tracks
            std::string fCaloHitLabel;  ///< label to find the right reco calo hits
            std::string fECALInstanceName; ///< product instance name for the ECAL

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry

            std::unique_ptr<rec::alg::KNNClusterAlg> fClusterAlgo; //Cluster algorithm
        };


        CaloClustering::CaloClustering(fhicl::ParameterSet const & p) : EDProducer{p}
        {
            fTrackLabel = p.get<std::string>("TrackLabel", "track");
            fCaloHitLabel = p.get<std::string>("CaloHitLabel", "calohit");
            fECALInstanceName =  p.get<std::string >("ECALInstanceName", "");

            fGeo     = gar::providerFrom<geo::GeometryGAr>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            //configure the cluster algorithm
            auto fClusterAlgoPars = p.get<fhicl::ParameterSet>("ClusterAlgPars");
            fClusterAlgo = std::make_unique<rec::alg::KNNClusterAlg>(fClusterAlgoPars);

            art::InputTag ecaltag(fCaloHitLabel, fECALInstanceName);
            consumes< std::vector<gar::rec::CaloHit> >(ecaltag);
            produces< std::vector<gar::rec::Cluster> >();
            produces< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> >();
            if (fClusterAlgo->usesTracks()) produces< art::Assns<gar::rec::Cluster, gar::rec::Track> >();
        }

        void CaloClustering::produce(art::Event & e)
        {
            //maps of the hit/track pointer to art ptr pointers
            std::unordered_map< const gar::rec::CaloHit*, art::Ptr<gar::rec::CaloHit> > hitMaptoArtPtr;
            std::unordered_map< const gar::rec::Track*, art::Ptr<gar::rec::Track> > trkMaptoArtPtr;

            //Collect the tracks to be passed to the algo.  As of Jun 2019, the algo doesn't need tracks.
            //But maybe someday it will!
            std::vector< art::Ptr<gar::rec::Track> > artTrk;
            this->CollectTracks(e, fTrackLabel, artTrk);

            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
            this->CollectHits(e, fCaloHitLabel, fECALInstanceName, artHits);

            //Prepare the hits for clustering (tag isolated hits and possible mip hits)
            //trkMaptoArtPtr will be empty if ( !fClusterAlgo->usesTracks() )
            fClusterAlgo->PrepareAlgo(artTrk, artHits, trkMaptoArtPtr, hitMaptoArtPtr);
            //Perform the clustering
            fClusterAlgo->DoClustering();
            //Get back the cluster results
            std::vector< gar::rec::Cluster* > ClusterVec = fClusterAlgo->GetFoundClusters();

            // make an art::PtrVector of the clusters
            std::unique_ptr< std::vector<gar::rec::Cluster> > ClusterCol(new std::vector<gar::rec::Cluster>);
            std::unique_ptr< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> > ClusterHitAssns(new art::Assns<gar::rec::Cluster, gar::rec::CaloHit>);
            std::unique_ptr< art::Assns<gar::rec::Cluster, gar::rec::Track> > ClusterTrackAssns(new art::Assns<gar::rec::Cluster, gar::rec::Track>);

            art::PtrMaker<gar::rec::Cluster> makeClusterPtr(e);

            MF_LOG_DEBUG("CaloClustering_module")
            << "Found " << ClusterVec.size() << " Clusters";

            //Copy the clusters to the collection
            for(auto const &it : ClusterVec)
            {
                MF_LOG_DEBUG("CaloClustering_module")
                << "Cluster has " << it->CalorimeterHits().size() << " calo hits";

                ClusterCol->emplace_back(*it);

                art::Ptr<gar::rec::Cluster> clusterPtr = makeClusterPtr(ClusterCol->size() - 1);

                //get list of hits associated to the cluster
                const std::vector< gar::rec::CaloHit* > hitVec = it->CalorimeterHits();
                for (size_t i=0; i<hitVec.size(); ++i)
                {
                    const gar::rec::CaloHit* pCaloHit = hitVec[i];
                    //Need to find the corresponding art ptr in the map
                    if(hitMaptoArtPtr.find(pCaloHit) != hitMaptoArtPtr.end())
                    {
                        // associate the hits to this cluster
                        art::Ptr<gar::rec::CaloHit> hitpointer = hitMaptoArtPtr[pCaloHit];
                        ClusterHitAssns->addSingle(clusterPtr, hitpointer);
                    }
                }

                if (fClusterAlgo->usesTracks()) {
				    //get list of track associated to the cluster
                    const std::vector< gar::rec::Track* > trkVec = it->Tracks();
                    for(const gar::rec::Track *const pTrack : trkVec)
                    {
                        //Need to find the corresponding art ptr in the map
                        if(trkMaptoArtPtr.find(pTrack) != trkMaptoArtPtr.end())
                        {
                            // associate the hits to this cluster
                            art::Ptr<gar::rec::Track> trkpointer = trkMaptoArtPtr[pTrack];
                            ClusterTrackAssns->addSingle(clusterPtr, trkpointer);
                        }
                    }
                }
            }

            e.put(std::move(ClusterCol));
            e.put(std::move(ClusterHitAssns));
            if (fClusterAlgo->usesTracks()) e.put(std::move(ClusterTrackAssns));

            return;
        }

        void CaloClustering::CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
        {
            art::Handle< std::vector<gar::rec::CaloHit> > theHits;
            evt.getByLabel(label, instance, theHits);

            if (!theHits.isValid())
            return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        void CaloClustering::CollectTracks(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::Track> > &trkVector)
        {
            art::Handle< std::vector<gar::rec::Track> > theTracks;
            evt.getByLabel(label, theTracks);

            if (!theTracks.isValid())
            return;

            for (unsigned int i = 0; i < theTracks->size(); ++i)
            {
                const art::Ptr<gar::rec::Track> track(theTracks, i);
                trkVector.push_back(track);
            }
        }

        DEFINE_ART_MODULE(CaloClustering)

    } // namespace rec
} // namespace gar
