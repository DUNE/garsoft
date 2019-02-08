//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"

#include "Reco/KNNClusterFinderAlg.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

#include "Utilities/KDTreeLinkerAlgoT.h"

namespace gar {
    namespace rec{

        //----------------------------------------------------------------------------
        InternalCaloHit::InternalCaloHit(HitRegion region, HitType type, unsigned int layer, CellGeometry cellgeom, CLHEP::Hep3Vector positionVector, CLHEP::Hep3Vector expectedDirection, const void* pParentAddress, float inputEnergy, float time, float cellThickness, float mipEquivalentEnergy, float cellsizeX, float cellsizeY)
        : m_hitRegion(region), m_hitType(type), m_layer(layer), m_cellGeometry(cellgeom), m_positionVector(positionVector), m_expectedDirection(expectedDirection), m_pParentAddress(pParentAddress), m_inputEnergy(inputEnergy), m_time(time), m_cellThickness(cellThickness), m_mipEquivalentEnergy(mipEquivalentEnergy), m_isIsolated(false), m_isPossibleMip(false), m_cellSizeX(cellsizeX), m_cellSizeY(cellsizeY)
        {
            m_cellLengthScale = this->CalculateCellLengthScale();
        }

        //----------------------------------------------------------------------------
        KNNClusterFinderAlg::KNNClusterFinderAlg(fhicl::ParameterSet const& pset)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            this->reconfigure(pset);

            ClearandResizeVectors();

            return;
        }

        //----------------------------------------------------------------------------
        KNNClusterFinderAlg::~KNNClusterFinderAlg()
        {
            delete m_hitNodes4D;
            delete m_hitsKdTree4D;
            return;
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::reconfigure(fhicl::ParameterSet const& pset)
        {
            fClusterAlgName = pset.get<std::string>("ClusterAlgName");
            fVerbose = pset.get<bool>("Verbose");
            feCalToMip = pset.get<float>("eCalToMip");
            feCalMipThreshold = pset.get<float>("eCalMipThreshold");

            //Algo params for the preparation of the calo hits
            m_caloHitMaxSeparation2 = pset.get<float>("caloHitMaxSeparation2", 100.f * 100.f);
            m_isolationCaloHitMaxSeparation2 = pset.get<float>("isolationCaloHitMaxSeparation2", 1000.f * 1000.f);
            m_isolationNLayers = pset.get<float>("isolationNLayers", 2);
            m_isolationCutDistance = pset.get<float>("isolationCutDistance", 200.f * 200.f);
            m_isolationSearchSafetyFactor = pset.get<float>("isolationSearchSafetyFactor", 2.f);
            m_isolationMaxNearbyHits = pset.get<float>("isolationMaxNearbyHits", 2);
            m_mipLikeMipCut = pset.get<float>("mipLikeMipCut", 5.f);
            m_mipNCellsForNearbyHit = pset.get<float>("mipNCellsForNearbyHit", 2);
            m_mipMaxNearbyHits = pset.get<float>("mipMaxNearbyHits", 1);

            m_hitNodes4D = new std::vector<HitKDNode4D>;
            m_hitsKdTree4D = new HitKDTree4D;

            //Cone Algo params for the clustering
            m_shouldUseIsolatedHits = pset.get<bool>("shouldUseIsolatedHits", false);
            m_nLayersSpannedForFit = pset.get<float>("nLayersSpannedForFit", 6);
            m_nLayersToFit = pset.get<float>("nLayersToFit", 8);
            m_nLayersToFitLowMipCut = pset.get<float>("nLayersToFitLowMipCut", 0.5f);
            m_nLayersToFitLowMipMultiplier = pset.get<float>("nLayersToFitLowMipMultiplier", 2);
            m_fitSuccessDotProductCut1 = pset.get<float>("fitSuccessDotProductCut1", 0.75f);
            m_fitSuccessChi2Cut1 = pset.get<float>("fitSuccessChi2Cut1", 5.0f);
            m_fitSuccessDotProductCut2 = pset.get<float>("fitSuccessDotProductCut2", 0.50f);
            m_fitSuccessChi2Cut2 = pset.get<float>("fitSuccessChi2Cut2", 2.5f);
            m_nLayersSpannedForApproxFit = pset.get<float>("nLayersSpannedForApproxFit", 10);

            return;
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::ClearandResizeVectors()
        {
            m_CaloHitList.clear();
            m_OrderedCaloHitList.clear();
            return;
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::FindClusters(std::vector<CaloHit> *allHits)
        {
            std::cout << "KNNClusterFinderAlg::FindClusters()" << std::endl;

            if( allHits->size() == 0 )
            {
                if (fVerbose) std::cout << "No hits received! exiting" << std::endl;
                return;
            }

            ClearandResizeVectors();

            //Prepare calo hits before clustering
            this->PrepareCaloHits(allHits);

            //Perform the clustering
            this->DoClustering();
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::PrepareCaloHits(std::vector<CaloHit> *allHits)
        {
            std::cout << "KNNClusterFinderAlg::PrepareCaloHits()" << std::endl;

            //Loop over all hits
            for(auto &iterCaloHit : *allHits)
            {
                unsigned int det_id = fGeo->getIDbyCellID(iterCaloHit.CellID(), "system");
                float absorberCorrection(1.);

                HitRegion region = (1 == det_id ) ? BARREL : ENDCAP;
                HitType type = ECAL;
                unsigned int layer = fGeo->getIDbyCellID(iterCaloHit.CellID(), "layer");
                CellGeometry cellGeometry = fGeo->isTile(det_id, layer) ? TILE : STRIP;
                const float* position(iterCaloHit.Position());
                const CLHEP::Hep3Vector positionVector(position[0], position[1], position[2]);
                float thickness = (TILE == cellGeometry) ? 10. : 5.;
                feCalToMip = (TILE == cellGeometry) ? feCalToMip*2 : feCalToMip;

                float mipEquivalentEnergy = iterCaloHit.Energy() * feCalToMip * absorberCorrection;

                if (mipEquivalentEnergy < feCalMipThreshold)
                continue;

                //TO DO
                float cellsizeX = fGeo->isTile(det_id, layer) ? 25 : 4;
                float cellsizeY = fGeo->isTile(det_id, layer) ? 25 : 4;

                const InternalCaloHit *internalCaloHit = new InternalCaloHit(region, type, layer, cellGeometry, positionVector, positionVector.unit(), (void*)&iterCaloHit, iterCaloHit.Energy(), iterCaloHit.Time(), thickness, mipEquivalentEnergy, cellsizeX, cellsizeY);

                m_CaloHitList.push_back(internalCaloHit);
            }

            //Initialize the KD Tree with the list of CaloHits
            this->InitializeKDTree(&m_CaloHitList);

            //Ordering the hit list by layers
            this->OrderCaloHitList(m_CaloHitList);

            for (OrderedCaloHitList::const_iterator iter = m_OrderedCaloHitList.begin(), iterEnd = m_OrderedCaloHitList.end(); iter != iterEnd; ++iter)
            {
                for (CaloHitList::iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
                {
                    //Check for isolated hits and possible MIP-like hits
                    this->CalculateCaloHitProperties(*hitIter, m_OrderedCaloHitList);
                }
            }
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::DoClustering()
        {
            std::cout << "KNNClusterFinderAlg::DoClustering()" << std::endl;

            ClusterVector clusterVector;
            m_hitsToClusters.clear();

            for (OrderedCaloHitList::const_iterator iter = m_OrderedCaloHitList.begin(), iterEnd = m_OrderedCaloHitList.end(); iter != iterEnd; ++iter)
            {
                // const unsigned int layer(iter->first);
                CaloHitVector relevantCaloHits;

                for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
                {
                    const InternalCaloHit *const pCaloHit = *hitIter;

                    if ( (m_shouldUseIsolatedHits || !pCaloHit->IsIsolated()) )
                    relevantCaloHits.push_back(pCaloHit);
                }

                ClusterFitResultMap clusterFitResultMap;
                this->GetCurrentClusterFitResults(clusterVector, clusterFitResultMap);
                // this->FindHitsInPreviousLayers(layer, relevantCaloHits, clusterFitResultMap, clusterVector);
                // this->FindHitsInSameLayer(layer, relevantCaloHits, clusterFitResultMap, clusterVector);
            }

            // this->RemoveEmptyClusters(clusterVector);

            //reset our kd trees and maps if everything turned out well
            m_hitsToClusters.clear();
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::OrderCaloHitList(const CaloHitList& rhs)
        {
            std::cout << "KNNClusterFinderAlg::OrderCaloHitList()" << std::endl;

            for (const InternalCaloHit *pCaloHit : rhs)
            {
                OrderedCaloHitList::iterator iter = m_OrderedCaloHitList.find(pCaloHit->GetLayer());

                if (m_OrderedCaloHitList.end() == iter)
                {
                    CaloHitList *const pCaloHitList = new std::list<const InternalCaloHit*>;
                    pCaloHitList->push_back(pCaloHit);

                    if (!(m_OrderedCaloHitList.insert(OrderedCaloHitList::value_type(pCaloHit->GetLayer(), pCaloHitList)).second))
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
        void KNNClusterFinderAlg::InitializeKDTree(const CaloHitList *pCaloHitList)
        {
            std::cout << "KNNClusterFinderAlg::InitializeKDTree()" << std::endl;

            m_hitsKdTree4D->clear();
            m_hitNodes4D->clear();
            util::KDTreeTesseract hitsBoundingRegion4D = util::fill_and_bound_4d_kd_tree(*pCaloHitList, *m_hitNodes4D);
            m_hitsKdTree4D->build(*m_hitNodes4D, hitsBoundingRegion4D);
            m_hitNodes4D->clear();
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::CalculateCaloHitProperties(const InternalCaloHit *const pCaloHit, const OrderedCaloHitList &orderedCaloHitList)
        {
            // Calculate number of adjacent pseudolayers to examine
            const unsigned int layer( pCaloHit->GetLayer() );
            const unsigned int isolationMaxLayer(layer + m_isolationNLayers);
            const unsigned int isolationMinLayer((layer < m_isolationNLayers) ? 0 : layer - m_isolationNLayers);

            // Initialize variables
            bool isIsolated = true;
            unsigned int isolationNearbyHits = 0;

            // Loop over adjacent pseudolayers
            for (unsigned int iLayer = isolationMinLayer; iLayer <= isolationMaxLayer; ++iLayer)
            {
                OrderedCaloHitList::const_iterator adjacentLayerIter = orderedCaloHitList.find(iLayer);

                if (orderedCaloHitList.end() == adjacentLayerIter)
                continue;

                // IsIsolated flag
                if (isIsolated && (isolationMinLayer <= iLayer) && (isolationMaxLayer >= iLayer))
                {
                    isolationNearbyHits += this->IsolationCountNearbyHits(iLayer, pCaloHit);
                    isIsolated = isolationNearbyHits < m_isolationMaxNearbyHits;
                }

                // Possible mip flag
                if (layer == iLayer)
                {
                    const CLHEP::Hep3Vector &positionVector(pCaloHit->GetPositionVector());

                    const float x(positionVector.x());
                    const float y(positionVector.y());

                    const float angularCorrection( (BARREL == pCaloHit->GetHitRegion()) ? positionVector.mag() / std::sqrt(x * x + y * y) : positionVector.mag() / std::fabs(positionVector.z()) );

                    if ((pCaloHit->GetMIPEquivalentEnergy() <= (m_mipLikeMipCut * angularCorrection)) && (m_mipMaxNearbyHits >= this->MipCountNearbyHits(iLayer, pCaloHit)))
                    pCaloHit->SetIsPossibleMIP(true);
                }
            }

            if (isIsolated)
            pCaloHit->SetIsIsolated(true);
        }

        //----------------------------------------------------------------------------
        unsigned int KNNClusterFinderAlg::IsolationCountNearbyHits(unsigned int searchLayer, const InternalCaloHit *const pCaloHit)
        {
            const CLHEP::Hep3Vector &positionVector(pCaloHit->GetPositionVector());
            const float positionMagnitudeSquared(positionVector.mag2());
            const float isolationCutDistanceSquared(m_isolationCutDistance);

            unsigned int nearbyHitsFound = 0;

            // construct the kd tree search
            CaloHitList nearby_hits;
            const float searchDistance(m_isolationSearchSafetyFactor * std::sqrt(isolationCutDistanceSquared));
            util::KDTreeTesseract searchRegionHits = util::build_4d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance, searchLayer);

            std::vector<HitKDNode4D> found;
            m_hitsKdTree4D->search(searchRegionHits, found);

            for (const auto &hit : found)
            {
                nearby_hits.push_back(hit.data);
            }

            for (CaloHitList::const_iterator iter = nearby_hits.begin(), iterEnd = nearby_hits.end(); iter != iterEnd; ++iter)
            {
                if (pCaloHit == *iter)
                continue;

                const CLHEP::Hep3Vector positionDifference(positionVector - (*iter)->GetPositionVector());
                const CLHEP::Hep3Vector crossProduct(positionVector.cross(positionDifference));

                if (positionDifference.mag2() > m_isolationCaloHitMaxSeparation2)
                continue;

                if ((crossProduct.mag2() / positionMagnitudeSquared) < isolationCutDistanceSquared)
                ++nearbyHitsFound;
            }

            return nearbyHitsFound;
        }

        //----------------------------------------------------------------------------
        unsigned int KNNClusterFinderAlg::MipCountNearbyHits(unsigned int searchLayer, const InternalCaloHit *const pCaloHit)
        {
            const float mipNCellsForNearbyHit(m_mipNCellsForNearbyHit + 0.5f);

            unsigned int nearbyHitsFound = 0;
            const CLHEP::Hep3Vector &positionVector(pCaloHit->GetPositionVector());
            const bool isHitInBarrelRegion(pCaloHit->GetHitRegion() == BARREL);

            // construct the kd tree search
            CaloHitList nearby_hits;
            const float searchDistance(std::sqrt(m_caloHitMaxSeparation2));
            util::KDTreeTesseract searchRegionHits = util::build_4d_kd_search_region(pCaloHit, searchDistance, searchDistance, searchDistance, searchLayer);

            std::vector<HitKDNode4D> found;
            m_hitsKdTree4D->search(searchRegionHits, found);

            for (const auto &hit : found)
            {
                nearby_hits.push_back(hit.data);
            }

            for (CaloHitList::const_iterator iter = nearby_hits.begin(), iterEnd = nearby_hits.end(); iter != iterEnd; ++iter)
            {
                if (pCaloHit == *iter)
                continue;

                const CLHEP::Hep3Vector positionDifference(positionVector - (*iter)->GetPositionVector());

                if (positionDifference.mag2() > m_caloHitMaxSeparation2)
                continue;

                const float cellLengthScale(pCaloHit->GetCellLengthScale());

                if (isHitInBarrelRegion)
                {
                    const float dX(std::fabs(positionDifference.x()));
                    const float dY(std::fabs(positionDifference.y()));
                    const float dZ(std::fabs(positionDifference.z()));
                    const float dPhi(std::sqrt(dX * dX + dY * dY));

                    if ((dZ < (mipNCellsForNearbyHit * cellLengthScale)) && (dPhi < (mipNCellsForNearbyHit * cellLengthScale)))
                    ++nearbyHitsFound;
                }
                else
                {
                    const float dX(std::fabs(positionDifference.x()));
                    const float dY(std::fabs(positionDifference.y()));

                    if ((dX < (mipNCellsForNearbyHit * cellLengthScale)) && (dY < (mipNCellsForNearbyHit * cellLengthScale)))
                    ++nearbyHitsFound;
                }
            }

            return nearbyHitsFound;
        }

        //----------------------------------------------------------------------------
        void KNNClusterFinderAlg::GetCurrentClusterFitResults(const ClusterVector &clusterVector, ClusterFitResultMap &clusterFitResultMap) const
        {
            std::cout << "KNNClusterFinderAlg::GetCurrentClusterFitResults()" << std::endl;

            if (!clusterFitResultMap.empty())
            return;

            for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
            {
                const InternalCluster *const pCluster = *iter;
                ClusterFitResult clusterFitResult;

                if (pCluster->GetNCaloHits() > 1)
                {
                    const unsigned int innerLayer(pCluster->GetInnerPseudoLayer());
                    const unsigned int outerLayer(pCluster->GetOuterPseudoLayer());
                    const unsigned int nLayersSpanned(outerLayer - innerLayer);

                    if (nLayersSpanned > m_nLayersSpannedForFit)
                    {
                        unsigned int nLayersToFit(m_nLayersToFit);

                        if (pCluster->GetMipFraction() - m_nLayersToFitLowMipCut < std::numeric_limits<float>::epsilon())
                        nLayersToFit *= m_nLayersToFitLowMipMultiplier;

                        const unsigned int startLayer( (nLayersSpanned > nLayersToFit) ? (outerLayer - nLayersToFit) : innerLayer);
                        (void) ClusterFitHelper::FitLayerCentroids(pCluster, startLayer, outerLayer, clusterFitResult);

                        if (clusterFitResult.IsFitSuccessful())
                        {
                            const float dotProduct(clusterFitResult.GetDirection().GetDotProduct(pCluster->GetInitialDirection()));
                            const float chi2(clusterFitResult.GetChi2());

                            if (((dotProduct < m_fitSuccessDotProductCut1) && (chi2 > m_fitSuccessChi2Cut1)) || ((dotProduct < m_fitSuccessDotProductCut2) && (chi2 > m_fitSuccessChi2Cut2)) )
                            {
                                clusterFitResult.SetSuccessFlag(false);
                            }
                        }
                    }
                    else if (nLayersSpanned > m_nLayersSpannedForApproxFit)
                    {
                        const CLHEP::Hep3Vector centroidChange(pCluster->GetCentroid(outerLayer) - pCluster->GetCentroid(innerLayer));
                        clusterFitResult.Reset();
                        clusterFitResult.SetDirection(centroidChange.unit());
                        clusterFitResult.SetSuccessFlag(true);
                    }
                }

                if (!clusterFitResultMap.insert(ClusterFitResultMap::value_type(pCluster, clusterFitResult)).second)
                return;
            }

            return;
        }


    } // namespace rec
} // namespace gar
