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
#include "RecoAlg/Cluster.h"

#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "RecoAlg/KNNClusterFinderAlg.h"

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
            void CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

            std::string fTrackLabel;  ///< label to find the reco tracks
            std::string fCaloHitLabel;  ///< label to find the right reco calo hits
            int fVerbosity;

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry

            std::unique_ptr<rec::alg::KNNClusterFinderAlg> fClusterAlgo; //Cluster algorithm
        };


        CaloClustering::CaloClustering(fhicl::ParameterSet const & p)
        {
            fTrackLabel = p.get<std::string>("TrackLabel", "track");
            fCaloHitLabel = p.get<std::string>("CaloHitLabel", "calohit");
            fVerbosity = p.get<int>("Verbosity", 0);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            //configure the cluster algorithm
            auto fClusterAlgoPars = p.get<fhicl::ParameterSet>("ClusterAlgPars");
            fClusterAlgo = std::make_unique<rec::alg::KNNClusterFinderAlg>(fClusterAlgoPars);

            produces< std::vector<gar::rec::Cluster> >();
            produces< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> >();
            produces< art::Assns<gar::rec::Cluster, gar::rec::Track> >();
        }

        void CaloClustering::produce(art::Event & e)
        {
            //maps of the hit/track pointer to art ptr pointers
            std::unordered_map< const gar::rec::CaloHit*, art::Ptr<gar::rec::CaloHit> > hitMaptoArtPtr;
            std::unordered_map< const gar::rec::Track*, art::Ptr<gar::rec::Track> > trkMaptoArtPtr;

            //Collect the tracks to be passed to the algo
            std::vector< art::Ptr<gar::rec::Track> > artTrk;
            this->CollectTracks(e, fTrackLabel, artTrk);

            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
            this->CollectHits(e, fCaloHitLabel, artHits);

            //Prepare the hits for clustering (tag isolated hits and possible mip hits)
            fClusterAlgo->PrepareAlgo(artTrk, artHits, trkMaptoArtPtr, hitMaptoArtPtr);
            //Perform the clustering
            fClusterAlgo->DoClustering();
            //Get back the cluster results
            std::vector< gar::rec::alg::Cluster* > ClusterVec = fClusterAlgo->GetFoundClusters();

            // make an art::PtrVector of the clusters
            std::unique_ptr< std::vector<gar::rec::Cluster> > ClusterCol(new std::vector<gar::rec::Cluster>);
            std::unique_ptr< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> > ClusterHitAssns(new art::Assns<gar::rec::Cluster, gar::rec::CaloHit>);
            std::unique_ptr< art::Assns<gar::rec::Cluster, gar::rec::Track> > ClusterTrackAssns(new art::Assns<gar::rec::Cluster, gar::rec::Track>);

            art::PtrMaker<gar::rec::Cluster> makeClusterPtr(e);

            if (fVerbosity>0) std::cout << "Found " << ClusterVec.size() << " Clusters" << std::endl;
            //Copy the clusters to the collection
            for(auto it : ClusterVec)
            {
                gar::rec::Cluster clus(it->getEnergy(), it->getDirection(), it->getInnerLayer(), it->getOuterLayer(), it->getNCaloHits(), it->getEnergyWeightedCenterOfGravity(), it->getEigenVectors(), it->getParticleId());

                if (fVerbosity>0) std::cout << "Cluster has " << clus.NCaloHits() << " calo hits" << std::endl;
                ClusterCol->push_back(clus);

                art::Ptr<gar::rec::Cluster> clusterPtr = makeClusterPtr(ClusterCol->size() - 1);

                //get list of hits associated to the cluster
                const auto hitList = it->getOrderedCaloHitList();
                for(auto it2 = hitList.begin(); it2!= hitList.end(); ++it2)
                {
                    std::list< const gar::rec::CaloHit* > *const pCaloHitList(it2->second);
                    for(const gar::rec::CaloHit *const pCaloHit : *pCaloHitList)
                    {
                        //Need to find the corresponding art ptr in the map
                        if(hitMaptoArtPtr.find(pCaloHit) != hitMaptoArtPtr.end())
                        {
                            // associate the hits to this cluster
                            art::Ptr<gar::rec::CaloHit> hitpointer = hitMaptoArtPtr[pCaloHit];
                            ClusterHitAssns->addSingle(clusterPtr, hitpointer);
                        }
                    }
                }

                //get list of track associated to the cluster
                const std::list< const gar::rec::Track* > trkList = it->getTrackList();
                for(const gar::rec::Track *const pTrack : trkList)
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

            e.put(std::move(ClusterCol));
            e.put(std::move(ClusterHitAssns));
            e.put(std::move(ClusterTrackAssns));

            return;
        }

        void CaloClustering::CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
        {
            art::Handle< std::vector<gar::rec::CaloHit> > theHits;
            evt.getByLabel(label, theHits);

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
