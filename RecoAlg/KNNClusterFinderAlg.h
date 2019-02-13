//
//  KNNClusterFinderAlg.h
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/include/LCClustering/ConeClusteringAlgorithm.h

#ifndef GAR_RECO_KNNClusterFinderAlg_h
#define GAR_RECO_KNNClusterFinderAlg_h

#include "art/Framework/Core/ModuleMacros.h"

#include "Geometry/GeometryCore.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

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
            typedef std::map<unsigned int, CaloHitList*> OrderedCaloHitList;
            typedef std::vector<const gar::rec::CaloHit*> CaloHitVector;
            typedef std::vector<gar::rec::Cluster*> ClusterVector;
            typedef std::list<gar::rec::Cluster *> ClusterList;
            typedef std::unordered_set<gar::rec::Cluster *> ClusterSet;

            class SortingHelper
            {
            public:

                static bool SortClustersByNHits(const gar::rec::Cluster *const pLhs, const gar::rec::Cluster *const pRhs);
            };

            class KNNClusterFinderAlg {

            public:

                KNNClusterFinderAlg(fhicl::ParameterSet const& pset);

                virtual ~KNNClusterFinderAlg();

                void reconfigure(fhicl::ParameterSet const& pset);

                void PrepareCaloHits(const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

                void DoClustering();

                ClusterVector GetFoundClusters() { return clusterVector; }

            private:
                typedef KDTreeLinkerAlgo<const gar::rec::CaloHit*, 4> HitKDTree;
                typedef KDTreeNodeInfoT<const gar::rec::CaloHit*, 4> HitKDNode;

                void ClearLists();

                void OrderCaloHitList(const CaloHitList& rhs);

                void InitializeKDTrees(const CaloHitList *const pCaloHitList);

                void FindHitsInPreviousLayers(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector &clusterVector);

                void FindHitsInSameLayer(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector &clusterVector);

                void GetGenericDistanceToHit(const gar::rec::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, const unsigned int searchLayer, float &genericDistance) const;

                void GetDistanceToHitInSameLayer(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList, float &distance) const;

                void GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList,
                const CLHEP::Hep3Vector &clusterDirection, float &distance) const;

                void GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CLHEP::Hep3Vector &clusterPosition, const CLHEP::Hep3Vector &clusterDirection, float &distance) const;

                void RemoveEmptyClusters(ClusterVector& clusterVector);

                gar::geo::GeometryCore const* fGeo;        ///< geometry information

                std::string fClusterAlgName;
                bool fVerbose;
                float fGeVtoMIP;

                float m_maxTrackSeedSeparation2;
                float m_genericDistanceCut;
                float m_layersToStepBack;

                float m_coneApproachMaxSeparation2;
                float m_maxClusterDirProjection;
                float m_minClusterDirProjection;
                float m_tanConeAngle;
                float m_additionalPadWidths;
                float m_sameLayerPadWidths;

                CaloHitList m_CaloHitList;
                OrderedCaloHitList m_OrderedCaloHitList;
                ClusterVector clusterVector;

                std::vector<HitKDNode> m_hitNodes;
                HitKDTree m_hitsKdTree;

                std::unordered_map<const gar::rec::CaloHit*, gar::rec::Cluster*> m_hitsToClusters;
            };

        } // namespace alg
    } // namespace rec
} //namespace gar

#endif /* GAR_RECO_KNNClusterFinderAlg_h */
