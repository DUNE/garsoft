//
//  SortingHelper.cxx
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/src/LCHelpers/SortingHelper.cc
//

#include "RecoAlg/SortingHelper.h"

#include "ReconstructionDataProducts/CaloHit.h"

#include <cmath>
#include <limits>

namespace gar {
    namespace rec{
        namespace alg{

            //----------------------------------------------------------------------------
            bool SortingHelper::SortClustersByNHits(const gar::rec::alg::Cluster *const pLhs, const gar::rec::alg::Cluster *const pRhs)
            {
                // NHits
                const unsigned int nCaloHitsLhs(pLhs->getNCaloHits()), nCaloHitsRhs(pRhs->getNCaloHits());

                if (nCaloHitsLhs != nCaloHitsRhs)
                return (nCaloHitsLhs > nCaloHitsRhs);

                // Energy
                const float energyLhs(pLhs->getEnergy()), energyRhs(pRhs->getEnergy());

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
            bool SortingHelper::SortClustersByInnerLayer(const gar::rec::alg::Cluster *const pLhs, const gar::rec::alg::Cluster *const pRhs)
            {
                const unsigned int nCaloHitsLhs(pLhs->getNCaloHits()), nCaloHitsRhs(pRhs->getNCaloHits());

                if ((nCaloHitsLhs > 0) && (nCaloHitsRhs > 0))
                {
                    const unsigned int innerLayerLhs(pLhs->getInnerLayer()), innerLayerRhs(pRhs->getInnerLayer());

                    if (innerLayerLhs != innerLayerRhs)
                    return (innerLayerLhs < innerLayerRhs);
                }

                return SortingHelper::SortClustersByNHits(pLhs, pRhs);
            }

        }
    }
}
