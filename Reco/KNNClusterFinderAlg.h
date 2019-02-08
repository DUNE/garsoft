//
//  KNNClusterFinderAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_RECO_KNNClusterFinderAlg_h
#define GAR_RECO_KNNClusterFinderAlg_h

#include "Geometry/GeometryCore.h"

#include "ReconstructionDataProducts/CaloHit.h"

#include <list>
#include <map>
#include <unordered_map>

namespace fhicl{
    class ParameterSet;
}

namespace util {
    template<typename, unsigned int> class KDTreeLinkerAlgo;
    template<typename, unsigned int> class KDTreeNodeInfoT;
}

namespace gar{
    namespace rec{

        enum HitRegion {
            BARREL,
            ENDCAP
        };

        enum HitType {
            ECAL,
            TRACKER
        };

        enum CellGeometry {
            TILE,
            STRIP
        };

        class InternalCaloHit {

        public:

            InternalCaloHit(HitRegion region, HitType type, unsigned int layer, CellGeometry cellgeom, CLHEP::Hep3Vector positionVector, CLHEP::Hep3Vector expectedDirection, const void* pParentAddress, float inputEnergy, float time, float cellThickness, float mipEquivalentEnergy, float cellsizeX, float cellsizeY);

            virtual ~InternalCaloHit() {}

            const HitRegion GetHitRegion() const { return m_hitRegion; }

            const HitType GetHitType() const { return m_hitType; }

            const unsigned int GetLayer() const { return m_layer; }

            const CellGeometry GetCellGeometry() const { return m_cellGeometry; }

            const CLHEP::Hep3Vector& GetPositionVector() const { return m_positionVector; }

            const CLHEP::Hep3Vector& GetExpectedDirection() const { return m_expectedDirection; }

            const void* GetParentAddress() const { return m_pParentAddress; }

            const float GetEnergy() const { return m_inputEnergy; }

            const float GetTime() const { return m_time; }

            const float GetCellThickness() const { return m_cellThickness; }

            const float GetMIPEquivalentEnergy() const { return m_mipEquivalentEnergy; }

            bool IsIsolated() const { return m_isIsolated; }

            bool IsPossibleMIP() const { return m_isPossibleMip; }

            void SetIsIsolated(bool IsIsolated) const { m_isIsolated = IsIsolated; }

            void SetIsPossibleMIP(bool IsPossibleMIP) const { m_isPossibleMip = IsPossibleMIP; }

            const float GetCellLengthScale() const { return m_cellLengthScale; }

            const float GetCellSizeX() const { return m_cellSizeX; }

            const float GetCellSizeY() const { return m_cellSizeY; }

        private:

            float CalculateCellLengthScale()
            {
                if (TILE == this->GetCellGeometry())
                {
                    return std::sqrt(this->GetCellSizeX() * this->GetCellSizeX());
                }
                else if (STRIP == this->GetCellGeometry())
                {
                    return std::sqrt(this->GetCellSizeX() * this->GetCellSizeY());
                }
                else
                {
                    return 1.;
                }
            }

            HitRegion m_hitRegion;
            HitType m_hitType;
            unsigned int m_layer;
            CellGeometry m_cellGeometry;
            CLHEP::Hep3Vector m_positionVector;
            CLHEP::Hep3Vector m_expectedDirection;
            const void* m_pParentAddress;
            float m_inputEnergy;
            float m_time;
            float m_cellThickness;
            float m_mipEquivalentEnergy;
            mutable bool m_isIsolated;
            mutable bool m_isPossibleMip;
            float m_cellLengthScale;
            float m_cellSizeX;
            float m_cellSizeY;
        };

        class InternalCluster {
        public:
            InternalCluster();
            virtual ~InternalCluster() {}

            unsigned int GetNCaloHits() const { return m_nCaloHits; }

        private:

            unsigned int m_nCaloHits;
        };

        class ClusterFitResult {
        public:
            ClusterFitResult();
            virtual ~ClusterFitResult() {}
        };

        typedef std::list< const InternalCaloHit* > CaloHitList;
        typedef std::map< unsigned int, CaloHitList* > OrderedCaloHitList;
        typedef std::vector< const InternalCaloHit* > CaloHitVector;
        typedef std::vector< const InternalCluster* > ClusterVector;

        class KNNClusterFinderAlg {

        public:

            typedef std::unordered_map<const InternalCluster*, const ClusterFitResult> ClusterFitResultMap;
            typedef util::KDTreeLinkerAlgo<const InternalCaloHit*, 4> HitKDTree4D;
            typedef util::KDTreeNodeInfoT<const InternalCaloHit*, 4> HitKDNode4D;

            KNNClusterFinderAlg(fhicl::ParameterSet const& pset);
            virtual ~KNNClusterFinderAlg();

            void reconfigure(fhicl::ParameterSet const& pset);

            void FindClusters(std::vector<CaloHit> *allHits);

        private:

            void ClearandResizeVectors();
            void PrepareCaloHits(std::vector<CaloHit> *allHits);
            void DoClustering();

            void InitializeKDTree(const CaloHitList *pCaloHitList);
            void OrderCaloHitList(const CaloHitList& rhs);
            void CalculateCaloHitProperties(const InternalCaloHit *const pCaloHit, const OrderedCaloHitList &orderedCaloHitList);
            unsigned int IsolationCountNearbyHits(unsigned int searchLayer, const InternalCaloHit *const pCaloHit);
            unsigned int MipCountNearbyHits(unsigned int searchLayer, const InternalCaloHit *const pCaloHit);

            void GetCurrentClusterFitResults(const ClusterVector &clusterVector, ClusterFitResultMap &clusterFitResultMap) const;

            gar::geo::GeometryCore const* fGeo;        ///< geometry information

            std::string fClusterAlgName;
            bool fVerbose;
            float feCalToMip;
            float feCalMipThreshold;

            //Hit preparation Algo params
            float m_isolationNLayers;
            float m_isolationMaxNearbyHits;
            float m_mipLikeMipCut;
            float m_mipMaxNearbyHits;
            float m_isolationCutDistance;
            float m_isolationSearchSafetyFactor;
            float m_isolationCaloHitMaxSeparation2;
            float m_mipNCellsForNearbyHit;
            float m_caloHitMaxSeparation2;

            //Clustering Algo params
            bool m_shouldUseIsolatedHits;
            float m_nLayersSpannedForFit;
            float m_nLayersToFit;
            float m_nLayersToFitLowMipCut;
            float m_nLayersToFitLowMipMultiplier;
            float m_fitSuccessDotProductCut1;
            float m_fitSuccessChi2Cut1;
            float m_fitSuccessDotProductCut2;
            float m_fitSuccessChi2Cut2;
            float m_nLayersSpannedForApproxFit;

            CaloHitList m_CaloHitList;
            OrderedCaloHitList m_OrderedCaloHitList;

            std::vector<HitKDNode4D>   *m_hitNodes4D;           ///< nodes for the KD tree (used for filling)
            HitKDTree4D                *m_hitsKdTree4D;         ///< the kd-tree itself, 4D in x,y,z,pseudolayer

            std::unordered_map<const InternalCaloHit*, const InternalCluster*> m_hitsToClusters;
        };

    } // namespace rec
} //namespace gar

#endif /* GAR_RECO_KNNClusterFinderAlg_h */
