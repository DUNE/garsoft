//
//  SortingHelper.h
//
//  Created by Eldwan Brianne on 10.02.2019
//  Based on Pandora Clustering Algo
//  https://github.com/PandoraPFA/LCContent/blob/master/include/LCHelpers/SortingHelper.h
//

#ifndef GAR_RECOALG_SortingHelper_h
#define GAR_RECOALG_SortingHelper_h

#include "RecoAlg/Cluster.h"

namespace gar{
  namespace rec{
    namespace alg{

      class SortingHelper
      {
      public:

        static bool SortClustersByNHits(const gar::rec::alg::Cluster *const pLhs, const gar::rec::alg::Cluster *const pRhs);

        static bool SortClustersByInnerLayer(const gar::rec::alg::Cluster *const pLhs, const gar::rec::alg::Cluster *const pRhs);
      };

    }
  }
}

#endif
