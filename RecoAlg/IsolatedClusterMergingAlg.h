//
//  IsolatedClusterMergingAlg.h
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/include/LCTopologicalAssociation/IsolatedHitMergingAlgorithm.h

#ifndef GAR_RECOALG_IsolatedClusterMergingAlg_h
#define GAR_RECOALG_IsolatedClusterMergingAlg_h

#include "RecoAlg/Cluster.h"
#include "ReconstructionDataProducts/CaloHit.h"

namespace fhicl{
    class ParameterSet;
}

namespace gar{
  namespace rec{
    namespace alg{

      typedef std::list<gar::rec::alg::Cluster*> ClusterList;
      typedef std::vector<gar::rec::alg::Cluster*> ClusterVector;

      class IsolatedClusterMergingAlg {

      public:

        IsolatedClusterMergingAlg(fhicl::ParameterSet const& pset);

        virtual ~IsolatedClusterMergingAlg();

        void reconfigure(fhicl::ParameterSet const& pset);

        void MergeIsolatedClusters(ClusterVector &clusVector);

      private:

        float GetDistanceToHit(const gar::rec::alg::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit) const;

        void DeleteClusterFromList(gar::rec::alg::Cluster *pCluster, ClusterVector &clusterVector);

        std::string fClusterAlgName;
        unsigned int m_minHitsInCluster;
        float m_maxRecombinationDistance;

      };

    }
  }
}

#endif
