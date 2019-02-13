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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
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
            void CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

            std::string fCaloHitLabel;  ///< label to find the right reco calo hits

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry

            std::unique_ptr<rec::alg::KNNClusterFinderAlg> fClusterAlgo; //Cluster algorithm
        };


        CaloClustering::CaloClustering(fhicl::ParameterSet const & p)
        {
            fCaloHitLabel = p.get<std::string>("CaloHitLabel", "calohit");

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            //configure the cluster algorithm
            auto fClusterAlgoPars = p.get<fhicl::ParameterSet>("ClusterAlgPars");
            fClusterAlgo = std::make_unique<rec::alg::KNNClusterFinderAlg>(fClusterAlgoPars);

            produces< std::vector<rec::Cluster> >();
        }

        void CaloClustering::produce(art::Event & e)
        {
            std::unique_ptr<std::vector<Cluster> > ClusterCol (new std::vector<Cluster> );

            std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
            this->CollectHits(e, fCaloHitLabel, artHits);

            //Prepare the hits for clustering (tag isolated hits and possible mip hits)
            fClusterAlgo->PrepareCaloHits(artHits);
            fClusterAlgo->DoClustering();
            std::vector< gar::rec::Cluster* > ClusterVec = fClusterAlgo->GetFoundClusters();

            mf::LogDebug("CaloClustering") << "Collected " << ClusterVec.size() << " Clusters" << std::endl;
            mf::LogDebug("CaloClustering") << "Number of hits before clustering " << artHits.size() << std::endl;

            //Copy the clusters to the collection
            int nhitclus = 0;
            int i = 1;
            for(auto &it : ClusterVec)
            {
                gar::rec::Cluster* clus = it;

                //check number of hits
                nhitclus += clus->NCaloHits();
                std::cout << "CaloClustering::produce(): Cluster " << i << " has " << clus->NCaloHits() << std::endl;

                ClusterCol->push_back(*clus);
                i++;
            }

            mf::LogDebug("CaloClustering") << "Number of hits after clustering " << nhitclus << std::endl;

            e.put(std::move(ClusterCol));

            return;
        }

        void CaloClustering::CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
        {
            art::Handle< std::vector<gar::rec::CaloHit> > theHits;
            evt.getByLabel(label, theHits);

            if (!theHits.isValid())
            {
                mf::LogDebug("CaloClustering") << "  Failed to find hits... " << std::endl;
                return;
            }
            else
            {
                mf::LogDebug("CaloClustering") << "  Found: " << theHits->size() << " Hits " << std::endl;
            }

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        DEFINE_ART_MODULE(CaloClustering)

    } // namespace rec
} // namespace gar
