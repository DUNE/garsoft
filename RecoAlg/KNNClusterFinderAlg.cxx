//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"

#include "RecoAlg/KNNClusterFinderAlg.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar {
    namespace rec{
        namespace alg{

            //----------------------------------------------------------------------------
            KNNClusterFinderAlg::KNNClusterFinderAlg(fhicl::ParameterSet const& pset)
            {
                fGeo = gar::providerFrom<geo::Geometry>();
                this->reconfigure(pset);

                ClearLists();

                return;
            }

            //----------------------------------------------------------------------------
            KNNClusterFinderAlg::~KNNClusterFinderAlg()
            {
                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                fClusterAlgName = pset.get<std::string>("ClusterAlgName");
                fVerbose = pset.get<bool>("Verbose", false);
                fGeVtoMIP = pset.get<float>("GeVtoMIP", 1.f);

                m_maxTrackSeedSeparation2 = 250.f * 250.f;
                m_genericDistanceCut = 1.f;
                m_maxClusterDirProjection = 200.f;
                m_layersToStepBack = 3;

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::ClearLists()
            {
                m_CaloHitList.clear();
                m_OrderedCaloHitList.clear();
                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::PrepareCaloHits(const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
            {
                std::cout << "KNNClusterFinderAlg::PrepareCaloHits()" << std::endl;

                //Loop over all hits
                for (std::vector< art::Ptr<gar::rec::CaloHit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
                {
                    const art::Ptr<gar::rec::CaloHit> hitPtr = *iter;
                    const gar::rec::CaloHit *hit = hitPtr.get();
                    m_CaloHitList.push_back(hit);
                }

                this->OrderCaloHitList(m_CaloHitList);

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::OrderCaloHitList(const CaloHitList& rhs)
            {
                std::cout << "KNNClusterFinderAlg::OrderCaloHitList()" << std::endl;

                for (const gar::rec::CaloHit* pCaloHit : rhs)
                {

                    long long int cellID = pCaloHit->CellID();
                    unsigned int layer = fGeo->getIDbyCellID(cellID, "layer");

                    OrderedCaloHitList::iterator iter = m_OrderedCaloHitList.find(layer);

                    if (m_OrderedCaloHitList.end() == iter)
                    {
                        CaloHitList *const pCaloHitList = new std::list< const gar::rec::CaloHit* >;
                        pCaloHitList->push_back(pCaloHit);

                        if (!(m_OrderedCaloHitList.insert(OrderedCaloHitList::value_type(layer, pCaloHitList)).second))
                        {
                            delete pCaloHitList;
                            return;
                        }
                    }
                    else
                    {
                        if (iter->second->end() != std::find(iter->second->begin(), iter->second->end(), pCaloHit))
                        return;

                        iter->second->push_back(pCaloHit);
                    }
                }
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::DoClustering()
            {
                std::cout << "KNNClusterFinderAlg::DoClustering()" << std::endl;

                //Initialise the KD Tree
                this->InitializeKDTrees(&m_CaloHitList);

                ClusterVector clusterVector;

                //Here perform the clustering on the ordered list of calo hits
                m_hitsToClusters.clear();
                for (OrderedCaloHitList::const_iterator iter = m_OrderedCaloHitList.begin(), iterEnd = m_OrderedCaloHitList.end(); iter != iterEnd; ++iter)
                {
                    const unsigned int layer(iter->first);
                    CaloHitVector relevantCaloHits;

                    for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
                    {
                        //Can be used to filter hits to perform the clustering
                        const CaloHit *const pCaloHit = *hitIter;
                        relevantCaloHits.push_back(pCaloHit);
                    }

                    //Try to find hits in the previous layers
                    this->FindHitsInPreviousLayers(layer, relevantCaloHits, clusterVector);

                    //Try to find hits in the same layer
                    this->FindHitsInSameLayer(layer, relevantCaloHits, clusterVector);
                }

                m_hitsKdTree.clear();
                m_hitsToClusters.clear();

                std::cout << "KNNClusterFinderAlg::DoClustering()" << " Created " << clusterVector.size() << " Cluster(s)" << std::endl;

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::InitializeKDTrees(const CaloHitList *const pCaloHitList)
            {
                m_hitsKdTree.clear();
                m_hitNodes.clear();
                KDTreeTesseract hitsBoundingRegion = fill_and_bound_4d_kd_tree(*pCaloHitList, m_hitNodes);
                m_hitsKdTree.build(m_hitNodes, hitsBoundingRegion);
                m_hitNodes.clear();
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::FindHitsInPreviousLayers(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector & /*clusterVector*/)
            {
                const float maxTrackSeedSeparation = std::sqrt(m_maxTrackSeedSeparation2);

                std::vector<HitKDNode> found_hits;
                ClusterSet nearby_clusters;

                for (const gar::rec::CaloHit *const pCaloHit : relevantCaloHits)
                {
                    gar::rec::Cluster *pBestCluster = nullptr;
                    float bestClusterEnergy(0.f);
                    float smallestGenericDistance(m_genericDistanceCut);
                    const float largestAllowedDistanceForSearch = std::max(maxTrackSeedSeparation, m_maxClusterDirProjection);

                    const unsigned int layersToStepBack(m_layersToStepBack);

                    // Associate with existing clusters in stepBack layers. If stepBackLayer == pseudoLayer, will examine track projections
                    for (unsigned int stepBackLayer = 1; (stepBackLayer <= layersToStepBack) && (stepBackLayer <= layer); ++stepBackLayer)
                    {
                        const unsigned int searchLayer(layer - stepBackLayer);

                        // now search for hits-in-clusters that would also satisfy the criteria
                        KDTreeTesseract searchRegionHits = build_4d_kd_search_region(pCaloHit, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch, largestAllowedDistanceForSearch, searchLayer);
                        m_hitsKdTree.search(searchRegionHits, found_hits);

                        for (auto &hit : found_hits)
                        {
                            auto assc_cluster = m_hitsToClusters.find(hit.data);
                            if (assc_cluster != m_hitsToClusters.end())
                            {
                                nearby_clusters.insert(assc_cluster->second);
                            }
                        }

                        found_hits.clear();

                        ClusterList nearbyClusterList(nearby_clusters.begin(), nearby_clusters.end());
                        // Need tp sort the list?
                        nearby_clusters.clear();

                        // Instead of using the full cluster list we use only those clusters that are found to be nearby according to the KD-tree
                        // ---- This can be optimized further for sure. (for instance having a match by KD-tree qualifies a ton of the loops later
                        // See if hit should be associated with any existing clusters
                        for (gar::rec::Cluster *pCluster : nearbyClusterList)
                        {
                            float genericDistance(std::numeric_limits<float>::max());
                            const float clusterEnergy(pCluster->Energy());

                            this->GetGenericDistanceToHit(pCluster, pCaloHit, searchLayer, genericDistance);

                            if ((genericDistance < smallestGenericDistance) || ((genericDistance == smallestGenericDistance) && (clusterEnergy > bestClusterEnergy)))
                            {
                                pBestCluster = pCluster;
                                bestClusterEnergy = clusterEnergy;
                                smallestGenericDistance = genericDistance;
                            }
                        }

                        // Add best hit found after completing examination of a stepback layer
                        if (nullptr != pBestCluster)
                        {
                            pBestCluster->AddToCluster(pCaloHit);
                            m_hitsToClusters.emplace(pCaloHit, pBestCluster);
                            break;
                        }

                    }//end step back loop
                }//end list loop

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::FindHitsInSameLayer(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector & clusterVector)
            {
                //keep a list of available hits with the most energetic available hit at the back
                std::list<unsigned> available_hits_in_layer;
                for (unsigned i = 0; i < relevantCaloHits.size(); ++i )
                available_hits_in_layer.push_back(i);

                //tactical cache hits -> hits
                std::unordered_multimap<const gar::rec::CaloHit*, const gar::rec::CaloHit*> hitsToHitsLocal;

                //pull out result caches to help keep memory locality
                std::vector<HitKDNode> found_hits;

                ClusterSet nearby_clusters;

                while (!available_hits_in_layer.empty())
                {
                    bool clustersModified = true;

                    while (clustersModified)
                    {
                        clustersModified = false;

                        for (auto iter = available_hits_in_layer.begin(); iter != available_hits_in_layer.end();)
                        {
                            // this his is assured to be usable by the lines above and algorithm course
                            const gar::rec::CaloHit *const pCaloHit = relevantCaloHits[*iter];

                            const float pad_search_width = 100.f;
                            const float hit_search_width = pad_search_width;

                            Cluster *pBestCluster = nullptr;
                            float bestClusterEnergy(0.f);
                            float smallestGenericDistance(m_genericDistanceCut);

                            // now search for hits-in-clusters that would also satisfy the criteria
                            auto hits_assc_cache = hitsToHitsLocal.find(pCaloHit);
                            if (hits_assc_cache != hitsToHitsLocal.end())
                            {
                                auto range = hitsToHitsLocal.equal_range(pCaloHit);
                                for (auto itr = range.first; itr != range.second; ++itr)
                                {
                                    auto assc_cluster = m_hitsToClusters.find(itr->second);
                                    if( assc_cluster != m_hitsToClusters.end() )
                                    {
                                        nearby_clusters.insert(assc_cluster->second);
                                    }
                                }
                            }
                            else
                            {
                                KDTreeTesseract searchRegionHits = build_4d_kd_search_region(pCaloHit, hit_search_width, hit_search_width, hit_search_width, layer);
                                m_hitsKdTree.search(searchRegionHits, found_hits);

                                for (auto &hit : found_hits)
                                {
                                    hitsToHitsLocal.emplace(pCaloHit, hit.data);
                                    auto assc_cluster = m_hitsToClusters.find(hit.data);
                                    if (assc_cluster != m_hitsToClusters.end())
                                    {
                                        nearby_clusters.insert(assc_cluster->second);
                                    }
                                }

                                found_hits.clear();
                            }

                            ClusterList nearbyClusterList(nearby_clusters.begin(), nearby_clusters.end());
                            //Sort the clusterlist?
                            nearby_clusters.clear();

                            // See if hit should be associated with any existing clusters
                            for (ClusterList::iterator clusterIter = nearbyClusterList.begin(), clusterIterEnd = nearbyClusterList.end(); clusterIter != clusterIterEnd; ++clusterIter)
                            {
                                gar::rec::Cluster *pCluster = *clusterIter;
                                float genericDistance(std::numeric_limits<float>::max());
                                const float clusterEnergy(pCluster->Energy());

                                this->GetGenericDistanceToHit(pCluster, pCaloHit, layer, genericDistance);

                                if ((genericDistance < smallestGenericDistance) || ((genericDistance == smallestGenericDistance) && (clusterEnergy > bestClusterEnergy)))
                                {
                                    pBestCluster = pCluster;
                                    bestClusterEnergy = clusterEnergy;
                                    smallestGenericDistance = genericDistance;
                                }
                            }

                            if (nullptr != pBestCluster)
                            {
                                pBestCluster->AddToCluster(pCaloHit);
                                m_hitsToClusters.emplace(pCaloHit, pBestCluster);
                                // remove this hit and advance to the next
                                iter = available_hits_in_layer.erase(iter);
                                clustersModified = true;
                            }
                            else
                            {
                                // otherwise advance the iterator
                                ++iter;
                            }
                        }
                    } // end clustersModified = true

                    // If there is no cluster within the search radius, seed a new cluster with this hit
                    if (!available_hits_in_layer.empty())
                    {
                        unsigned index = *(available_hits_in_layer.begin());
                        available_hits_in_layer.pop_front();

                        const gar::rec::CaloHit *const  pCaloHit = relevantCaloHits[index];

                        // hit is assured to be valid
                        Cluster *pCluster = new gar::rec::Cluster();
                        pCluster->AddToCluster(pCaloHit);
                        clusterVector.push_back(pCluster);
                        m_hitsToClusters.emplace(pCaloHit, pCluster);
                    }
                }// end while (!available_hits_in_layer.empty())

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::GetGenericDistanceToHit(const gar::rec::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, const unsigned int searchLayer, float &genericDistance) const
            {
                float initialDirectionDistance(std::numeric_limits<float>::max());
                float currentDirectionDistance(std::numeric_limits<float>::max());

                const CaloHitList pClusterCaloHitList = pCluster->CaloHitList();

                if (searchLayer == pCaloHit->GetLayer())
                {
                    return this->GetDistanceToHitInSameLayer(pCaloHit, &pClusterCaloHitList, genericDistance);
                }

                //Use initial direction
                this->GetConeApproachDistanceToHit(pCaloHit, &pClusterCaloHitList, pCluster->GetInitialDirection(), initialDirectionDistance);

                //Use current direction
                if(pCluster->CaloHitList().size() > 1)
                this->GetConeApproachDistanceToHit(pCaloHit, pClusterCaloHitList, pCluster->GetDirection(), currentDirectionDistance)

                const float smallestDistance( std::min(initialDirectionDistance, currentDirectionDistance) );

                if (smallestDistance < genericDistance)
                {
                    genericDistance = smallestDistance;

                    if (std::numeric_limits<float>::max() != genericDistance)
                    return;
                }

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::GetDistanceToHitInSameLayer(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList, float &distance) const
            {
                const float dCut(100.f);
                const CLHEP::Hep3Vector &hitPosition(pCaloHit->GetPositionVector());

                bool hitFound(false);
                float smallestDistanceSquared(std::numeric_limits<float>::max());
                const float rDCutSquared(1.f / (dCut * dCut));

                for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
                {
                    const gar::rec::CaloHit *const pHitInCluster = *iter;
                    const CLHEP::Hep3Vector &hitInClusterPosition(pHitInCluster->GetPositionVector());
                    const float separationSquared((hitPosition - hitInClusterPosition).mag2());
                    const float hitDistanceSquared(separationSquared * rDCutSquared);

                    if (hitDistanceSquared < smallestDistanceSquared)
                    {
                        smallestDistanceSquared = hitDistanceSquared;
                        hitFound = true;
                    }
                }

                if (!hitFound)
                return;

                distance = std::sqrt(smallestDistanceSquared);
                return;
            }

        } // namespace alg
    } // namespace rec
} // namespace gar
