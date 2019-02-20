//
//  IsolatedClusterMergingAlg.cxx
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/src/LCTopologicalAssociation/IsolatedHitMergingAlgorithm.cc
//

#include "fhiclcpp/ParameterSet.h"

#include "RecoAlg/IsolatedClusterMergingAlg.h"
#include "RecoAlg/SortingHelper.h"

#include <algorithm>

namespace gar {
    namespace rec{
        namespace alg{

            //----------------------------------------------------------------------------
            IsolatedClusterMergingAlg::IsolatedClusterMergingAlg(fhicl::ParameterSet const& pset)
            {
                this->reconfigure(pset);

                return;
            }

            //----------------------------------------------------------------------------
            IsolatedClusterMergingAlg::~IsolatedClusterMergingAlg()
            {
                return;
            }

            //----------------------------------------------------------------------------
            void IsolatedClusterMergingAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                fClusterAlgName = pset.get<std::string>("ClusterAlgName");

                m_minHitsInCluster = pset.get<unsigned int>("MinHitsInCluster", 4);
                m_maxRecombinationDistance = pset.get<float>("MaxRecombinationDistance", 250.f);

                return;
            }

            //----------------------------------------------------------------------------
            void IsolatedClusterMergingAlg::MergeIsolatedClusters(ClusterVector &clusVector)
            {
                //Copy the vector locally
                ClusterVector clusterVector(clusVector.begin(), clusVector.end());
                //Sort by layers
                std::sort(clusterVector.begin(), clusterVector.end(), SortingHelper::SortClustersByInnerLayer);

                //find "small" clusters, below threshold number of calo hits, delete them and associate hits with closest cluster
                for (ClusterVector::iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
                {
                    gar::rec::alg::Cluster *pClusterToDelete = *iterI;

                    //avoid deleted clusters
                    if (nullptr == pClusterToDelete)
                    continue;

                    const unsigned int nCaloHits(pClusterToDelete->getNCaloHits());

                    if (nCaloHits > m_minHitsInCluster)
                    continue;

                    //Get the list of hits in that cluster
                    OrderedCaloHitList orderedCaloHitList = pClusterToDelete->getOrderedCaloHitList();

                    //Delete the cluster from original list
                    this->DeleteClusterFromList(pClusterToDelete, clusVector);

                    // Redistribute hits that used to be in cluster I amongst other clusters
                    for (OrderedCaloHitList::const_iterator it = orderedCaloHitList.begin(); it != orderedCaloHitList.end(); ++it)
                    {
                        CaloHitList *const pCaloHitList(it->second);

                        //Loop over the list of calo hits of the deleted cluster
                        for(const gar::rec::CaloHit *const pCaloHit : *pCaloHitList)
                        {
                            gar::rec::alg::Cluster *pBestHostCluster(nullptr);
                            float bestHostClusterEnergy(0.);
                            float minDistance(m_maxRecombinationDistance);

                            // Find the most appropriate cluster for this newly-available hit
                            for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
                            {
                                gar::rec::alg::Cluster *const pNewHostCluster = *iterJ;

                                if (nullptr == pNewHostCluster)
                                continue;

                                if (pNewHostCluster->getNCaloHits() < nCaloHits)
                                continue;

                                const float distance(this->GetDistanceToHit(pNewHostCluster, pCaloHit));
                                const float hostClusterEnergy(pNewHostCluster->getEnergy());

                                // In event of equidistant host candidates, choose highest energy cluster
                                if ((distance < minDistance) || ((distance == minDistance) && (hostClusterEnergy > bestHostClusterEnergy)))
                                {
                                    minDistance = distance;
                                    pBestHostCluster = pNewHostCluster;
                                    bestHostClusterEnergy = hostClusterEnergy;
                                }
                            }

                            if (nullptr != pBestHostCluster)
                            {
                                //Add the hit to the found cluster
                                pBestHostCluster->AddToCluster(pCaloHit);
                            }
                        }
                    }
                }

                return;
            }

            //----------------------------------------------------------------------------
            float IsolatedClusterMergingAlg::GetDistanceToHit(const gar::rec::alg::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit) const
            {
                float minDistanceSquared(std::numeric_limits<float>::max());
                const CLHEP::Hep3Vector &hitPosition(pCaloHit->GetPositionVector());
                const OrderedCaloHitList &orderedCaloHitList(pCluster->getOrderedCaloHitList());

                for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
                {
                    const float distanceSquared((pCluster->GetCentroid(iter->first) - hitPosition).mag2());

                    if (distanceSquared < minDistanceSquared)
                    minDistanceSquared = distanceSquared;
                }

                if (minDistanceSquared < std::numeric_limits<float>::max())
                return std::sqrt(minDistanceSquared);

                return std::numeric_limits<float>::max();
            }

            //----------------------------------------------------------------------------
            void IsolatedClusterMergingAlg::DeleteClusterFromList(gar::rec::alg::Cluster *pCluster, ClusterVector &clusterVector)
            {
                ClusterVector::iterator iter = std::find(clusterVector.begin(), clusterVector.end(), pCluster);

                if(clusterVector.end() == iter)
                return;

                //Delete the cluster and set it to nullptr
                clusterVector.erase(iter);

                return;
            }

        }
    }
}
