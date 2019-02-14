//
//  KNNClusterFinderAlg.h
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/include/LCClustering/ConeClusteringAlgorithm.h

#ifndef GAR_RECOALG_KNNClusterFinderAlg_h
#define GAR_RECOALG_KNNClusterFinderAlg_h

#include "art/Framework/Core/ModuleMacros.h"

#include "Geometry/GeometryCore.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "RecoAlg/Cluster.h"
#include "ReconstructionDataProducts/Track.h"

#include "RecoAlg/KDTreeAlgo.h"

#include <list>
#include <map>
#include <unordered_map>
#include <unordered_set>

namespace fhicl{
    class ParameterSet;
}

namespace gar{
    namespace rec{
        namespace alg{

            typedef std::list<const gar::rec::CaloHit*> CaloHitList;
            typedef std::list<const gar::rec::Track *> TrackList;
            typedef std::map<unsigned int, CaloHitList*> OrderedCaloHitList;
            typedef std::vector<const gar::rec::CaloHit*> CaloHitVector;
            typedef std::vector<gar::rec::alg::Cluster*> ClusterVector;
            typedef std::list<gar::rec::alg::Cluster *> ClusterList;
            typedef std::unordered_set<gar::rec::alg::Cluster *> ClusterSet;

            class SortingHelper
            {
            public:

                static bool SortClustersByNHits(const gar::rec::alg::Cluster *const pLhs, const gar::rec::alg::Cluster *const pRhs);
            };

            class KNNClusterFinderAlg {

            public:

                KNNClusterFinderAlg(fhicl::ParameterSet const& pset);

                virtual ~KNNClusterFinderAlg();

                void reconfigure(fhicl::ParameterSet const& pset);

                void PrepareAlgo(const std::vector< art::Ptr<gar::rec::Track> > &trkVector, const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector, std::unordered_map< const gar::rec::Track*, art::Ptr<gar::rec::Track> > &trkMaptoArtPtr, std::unordered_map< const gar::rec::CaloHit*, art::Ptr<gar::rec::CaloHit> > &hitMaptoArtPtr);

                void DoClustering();

                ClusterVector GetFoundClusters() { return clusterVector; }

            private:
                typedef KDTreeLinkerAlgo<const gar::rec::Track*, 3> TrackKDTree;
                typedef KDTreeNodeInfoT<const gar::rec::Track*, 3> TrackKDNode;

                typedef KDTreeLinkerAlgo<const gar::rec::CaloHit*, 4> HitKDTree;
                typedef KDTreeNodeInfoT<const gar::rec::CaloHit*, 4> HitKDNode;

                void ClearLists();

                void OrderCaloHitList(const CaloHitList& rhs);

                void InitializeKDTrees(const TrackList *const pTrackList, const CaloHitList *const pCaloHitList);

                void SeedClustersWithTracks(const TrackList *const pTrackList, ClusterVector &clusterVector);

                void FindHitsInPreviousLayers(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector &clusterVector);

                void FindHitsInSameLayer(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector &clusterVector);

                void GetGenericDistanceToHit(const gar::rec::alg::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, const unsigned int searchLayer, float &genericDistance) const;

                void GetDistanceToHitInSameLayer(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList, float &distance) const;

                void GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList,
                const CLHEP::Hep3Vector &clusterDirection, float &distance) const;

                void GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CLHEP::Hep3Vector &clusterPosition, const CLHEP::Hep3Vector &clusterDirection, float &distance) const;

                void GetDistanceToTrackSeed(const gar::rec::alg::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, unsigned int searchLayer, float &distance) const;

                void GetDistanceToTrackSeed(const gar::rec::alg::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, float &distance) const;

                void RemoveEmptyClusters(ClusterVector& clusterVector);

                gar::geo::GeometryCore const* fGeo;        ///< geometry information

                std::string fClusterAlgName;

                unsigned int m_clusterSeedStrategy;
                float m_maxTrackSeedSeparation2;
                float m_genericDistanceCut;
                unsigned int m_layersToStepBack;

                float m_coneApproachMaxSeparation2;
                float m_maxClusterDirProjection;
                float m_minClusterDirProjection;
                float m_tanConeAngle;
                float m_additionalPadWidths;
                float m_sameLayerPadWidths;

                bool m_shouldUseTrackSeed;
                bool m_shouldFollowInitialDirection;
                unsigned int m_trackSeedCutOffLayer;
                unsigned int m_firstLayer;
                unsigned int m_maxLayersToTrackSeed;
                float m_trackPathWidth;
                unsigned int m_maxLayersToTrackLikeHit;

                CaloHitList m_CaloHitList;
                TrackList m_TrackList;

                OrderedCaloHitList m_OrderedCaloHitList;
                ClusterVector clusterVector;

                std::vector<TrackKDNode> m_trackNodes;
                TrackKDTree m_tracksKdTree;

                std::vector<HitKDNode> m_hitNodes;
                HitKDTree m_hitsKdTree;

                std::unordered_map<const gar::rec::CaloHit*, gar::rec::alg::Cluster*> m_hitsToClusters;
                std::unordered_map<const gar::rec::Track*, gar::rec::alg::Cluster*> m_tracksToClusters;
            };

        } // namespace alg
    } // namespace rec
} //namespace gar

#endif /* GAR_RECO_KNNClusterFinderAlg_h */
