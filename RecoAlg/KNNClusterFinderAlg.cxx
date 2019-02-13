//
//  KNNClusterFinderAlg.cxx
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/src/LCClustering/ConeClusteringAlgorithm.cc
//

#include "fhiclcpp/ParameterSet.h"

#include "RecoAlg/KNNClusterFinderAlg.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

#include <algorithm>

namespace gar {
    namespace rec{
        namespace alg{

            bool SortingHelper::SortClustersByNHits(const gar::rec::Cluster *const pLhs, const gar::rec::Cluster *const pRhs)
            {
                // NHits
                const unsigned int nCaloHitsLhs(pLhs->NCaloHits()), nCaloHitsRhs(pRhs->NCaloHits());

                if (nCaloHitsLhs != nCaloHitsRhs)
                return (nCaloHitsLhs > nCaloHitsRhs);

                // Energy
                const float energyLhs(pLhs->Energy()), energyRhs(pRhs->Energy());

                if (std::fabs(energyLhs - energyRhs) > std::numeric_limits<float>::epsilon())
                return (energyLhs > energyRhs);

                // Final attempt to distinguish
                if ((nCaloHitsLhs > 0) && (nCaloHitsRhs > 0))
                {
                    const gar::rec::CaloHit *const pFirstHitLhs((pLhs->getOrderedCaloHitList().begin())->second->front());
                    const gar::rec::CaloHit *const pFirstHitRhs((pRhs->getOrderedCaloHitList().begin())->second->front());

                    return (*pFirstHitLhs < *pFirstHitRhs);
                }

                throw cet::exception("SortingHelper::SortClustersByNHits") << "Can't sort cluster by nHits";
            }

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

                //Algorithm parameters
                m_clusterSeedStrategy = pset.get<unsigned int>("ClusterSeedStrategy", 2);
                m_maxTrackSeedSeparation2 = pset.get<float>("MaxTrackSeedSeparation2", 25.f * 25.f);
                m_genericDistanceCut = pset.get<float>("GenericDistanceCut", 1.f);
                m_layersToStepBack = pset.get<unsigned int>("LayersToStepBack", 3);

                m_coneApproachMaxSeparation2 = pset.get<float>("ConeApproachMaxSeparation2", 100.f * 100.f);
                m_minClusterDirProjection = pset.get<float>("MinClusterDirProjection", -1.f);
                m_maxClusterDirProjection = pset.get<float>("MaxClusterDirProjection", 20.f);
                m_tanConeAngle = pset.get<float>("TanConeAngle", 0.5f);
                m_additionalPadWidths = pset.get<float>("AdditionalPadWidths", 2.5f);
                m_sameLayerPadWidths = pset.get<float>("SameLayerPadWidths", 2.8f);

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::ClearLists()
            {
                clusterVector.clear();
                m_CaloHitList.clear();
                m_TrackList.clear();
                m_OrderedCaloHitList.clear();
                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::PrepareAlgo(const std::vector< art::Ptr<gar::rec::Track> > &trkVector, const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
            {
                std::cout << "KNNClusterFinderAlg::PrepareAlgo()" << std::endl;

                //Clear the lists
                ClearLists();

                //Loop over all tracks
                for (std::vector< art::Ptr<gar::rec::Track> >::const_iterator iter = trkVector.begin(), iterEnd = trkVector.end(); iter != iterEnd; ++iter)
                {
                    const art::Ptr<gar::rec::Track> trkPtr = *iter;
                    const gar::rec::Track *track = trkPtr.get();
                    m_TrackList.push_back(track);
                }

                //Loop over all hits
                for (std::vector< art::Ptr<gar::rec::CaloHit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
                {
                    const art::Ptr<gar::rec::CaloHit> hitPtr = *iter;
                    const gar::rec::CaloHit *hit = hitPtr.get();
                    m_CaloHitList.push_back(hit);
                }

                //Order the hit list by layers
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

                //Here perform the clustering on the ordered list of calo hits
                clusterVector.clear();
                m_hitsToClusters.clear();
                m_tracksToClusters.clear();

                this->SeedClustersWithTracks(&m_TrackList, clusterVector);

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

                this->RemoveEmptyClusters(clusterVector);

                //Sort the cluster by number of hits
                std::sort(clusterVector.begin(), clusterVector.end(), SortingHelper::SortClustersByNHits);

                m_hitsKdTree.clear();
                m_hitsToClusters.clear();
                m_tracksToClusters.clear();

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::InitializeKDTrees(const CaloHitList *const pCaloHitList)
            {
                std::cout << "KNNClusterFinderAlg::InitializeKDTrees()" << std::endl;

                m_hitsKdTree.clear();
                m_hitNodes.clear();
                KDTreeTesseract hitsBoundingRegion = fill_and_bound_4d_kd_tree(*pCaloHitList, m_hitNodes);
                m_hitsKdTree.build(m_hitNodes, hitsBoundingRegion);
                m_hitNodes.clear();
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::SeedClustersWithTracks(const TrackList *const pTrackList, ClusterVector &clusterVector)
            {
                std::cout << "KNNClusterFinderAlg::SeedClustersWithTracks()" << std::endl;

                if (0 == m_clusterSeedStrategy)
                return;

                if (nullptr == pTrackList)
                return;

                for (TrackList::const_iterator iter = pTrackList->begin(), iterEnd = pTrackList->end(); iter != iterEnd; ++iter)
                {
                    const Track *const pTrack = *iter;

                    bool useTrack(false);

                    if (2 == m_clusterSeedStrategy)
                    useTrack = true;

                    if (useTrack)
                    {
                        Cluster *pCluster = new gar::rec::Cluster();
                        pCluster->AddToCluster(pTrack);
                        clusterVector.push_back(pCluster);
                        m_tracksToClusters.emplace(pTrack,pCluster);
                    }
                }

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

                            const float hit_search_width = m_sameLayerPadWidths * pCaloHit->GetCellLengthScale();

                            Cluster *pBestCluster = nullptr;
                            float bestClusterEnergy(0.f);
                            float smallestGenericDistance(m_genericDistanceCut);

                            //search for hits-in-clusters that would satisfy the criteria
                            auto hits_assc_cache = hitsToHitsLocal.find(pCaloHit);
                            //if the cache is not empty
                            if (hits_assc_cache != hitsToHitsLocal.end())
                            {
                                //get the range
                                auto range = hitsToHitsLocal.equal_range(pCaloHit);
                                for (auto itr = range.first; itr != range.second; ++itr)
                                {
                                    //get the cluster associated to the hits in the equal range
                                    auto assc_cluster = m_hitsToClusters.find(itr->second);
                                    if( assc_cluster != m_hitsToClusters.end() )
                                    nearby_clusters.insert(assc_cluster->second);
                                }
                            }
                            else
                            {
                                //Search in the region hit_search_width * hit_search_width in that layer
                                KDTreeTesseract searchRegionHits = build_4d_kd_search_region(pCaloHit, hit_search_width, hit_search_width, hit_search_width, layer);
                                m_hitsKdTree.search(searchRegionHits, found_hits);

                                //loop over the found hits
                                for (auto &hit : found_hits)
                                {
                                    hitsToHitsLocal.emplace(pCaloHit, hit.data);
                                    auto assc_cluster = m_hitsToClusters.find(hit.data);
                                    //check if the hit is already associated to a cluster
                                    if (assc_cluster != m_hitsToClusters.end())
                                    nearby_clusters.insert(assc_cluster->second);
                                }

                                found_hits.clear();
                            }

                            //copy the list of nearby clusters
                            ClusterList nearbyClusterList(nearby_clusters.begin(), nearby_clusters.end());
                            //sort the list
                            nearbyClusterList.sort(SortingHelper::SortClustersByNHits);
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
                                //Add the hit to the cluster
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
                        }// end loop over available_hits_in_layer
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
            void KNNClusterFinderAlg::FindHitsInPreviousLayers(unsigned int layer, const CaloHitVector &relevantCaloHits, ClusterVector & /*clusterVector*/)
            {
                const float maxTrackSeedSeparation = std::sqrt(m_maxTrackSeedSeparation2);

                std::vector<HitKDNode> found_hits;
                ClusterSet nearby_clusters;

                for (const gar::rec::CaloHit *const pCaloHit : relevantCaloHits)
                {
                    const float additionalPadWidths = m_additionalPadWidths * pCaloHit->GetCellLengthScale();

                    gar::rec::Cluster *pBestCluster = nullptr;
                    float bestClusterEnergy(0.f);
                    float smallestGenericDistance(m_genericDistanceCut);
                    const float largestAllowedDistanceForSearch = std::max(maxTrackSeedSeparation, m_maxClusterDirProjection + additionalPadWidths);

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
                        nearbyClusterList.sort(SortingHelper::SortClustersByNHits);
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
            void KNNClusterFinderAlg::GetGenericDistanceToHit(const gar::rec::Cluster *const pCluster, const gar::rec::CaloHit *const pCaloHit, const unsigned int searchLayer, float &genericDistance) const
            {
                // Check that cluster is occupied in the searchlayer and is reasonably compatible with calo hit
                OrderedCaloHitList::const_iterator clusterHitListIter = pCluster->getOrderedCaloHitList().find(searchLayer);

                if (pCluster->getOrderedCaloHitList().end() == clusterHitListIter)
                return;

                float initialDirectionDistance(std::numeric_limits<float>::max());
                float currentDirectionDistance(std::numeric_limits<float>::max());

                const CaloHitList *pClusterCaloHitList = clusterHitListIter->second;

                if (searchLayer == pCaloHit->GetLayer())
                {
                    return this->GetDistanceToHitInSameLayer(pCaloHit, pClusterCaloHitList, genericDistance);
                }

                //Use initial direction
                this->GetConeApproachDistanceToHit(pCaloHit, pClusterCaloHitList, pCluster->InitialDirection(), initialDirectionDistance);

                //Use current direction
                if(pCluster->NCaloHits() > 1)
                this->GetConeApproachDistanceToHit(pCaloHit, pClusterCaloHitList, pCluster->Direction(), currentDirectionDistance);

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
                const float dCut = m_sameLayerPadWidths * pCaloHit->GetCellLengthScale();
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

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CaloHitList *const pCaloHitList,
            const CLHEP::Hep3Vector &clusterDirection, float &distance) const
            {
                bool hitFound(false);
                float smallestDistance(std::numeric_limits<float>::max());

                for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
                {
                    const gar::rec::CaloHit *const pHitInCluster = *iter;
                    float hitDistance(std::numeric_limits<float>::max());

                    this->GetConeApproachDistanceToHit(pCaloHit, pHitInCluster->GetPositionVector(), clusterDirection, hitDistance);

                    if (hitDistance < smallestDistance)
                    {
                        smallestDistance = hitDistance;
                        hitFound = true;
                    }
                }

                if (!hitFound)
                return;

                distance = smallestDistance;
                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterFinderAlg::GetConeApproachDistanceToHit(const gar::rec::CaloHit *const pCaloHit, const CLHEP::Hep3Vector &clusterPosition, const CLHEP::Hep3Vector &clusterDirection, float &distance) const
            {
                const CLHEP::Hep3Vector &hitPosition(pCaloHit->GetPositionVector());
                const CLHEP::Hep3Vector positionDifference(hitPosition - clusterPosition);

                if (positionDifference.mag2() > m_coneApproachMaxSeparation2)
                return;

                const float dAlong(clusterDirection.dot(positionDifference));

                if ((dAlong < m_maxClusterDirProjection) && (dAlong > m_minClusterDirProjection))
                {
                    const float dCut = (std::fabs(dAlong) * m_tanConeAngle) + (m_additionalPadWidths * pCaloHit->GetCellLengthScale());

                    if (dCut < std::numeric_limits<float>::epsilon())
                    return;

                    const float dPerp (clusterDirection.cross(positionDifference).mag());

                    distance = dPerp / dCut;
                    return;
                }

                return;
            }

            void KNNClusterFinderAlg::RemoveEmptyClusters(ClusterVector& clusterVector)
            {
                for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
                {
                    if (0 == (*iter)->NCaloHits())
                    {
                        //Erase the cluster
                        clusterVector.erase(iter);
                    }
                }
            }

        } // namespace alg
    } // namespace rec
} // namespace gar
