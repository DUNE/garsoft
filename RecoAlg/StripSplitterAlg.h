//
//  StripSplitterAlg.h
//
//  Created by Eldwan Brianne on 18.04.2019
//  Based on DDStripSplitter
//  https://github.com/iLCSoft/MarlinReco/blob/master/Clustering/hybridEcalSplitter/include/DDStripSplitter.h

#ifndef GAR_RECOALG_StripSplitterAlg_h
#define GAR_RECOALG_StripSplitterAlg_h

#include "art/Framework/Core/ModuleMacros.h"
#include "Geometry/GeometryCore.h"

#include "ReconstructionDataProducts/CaloHit.h"

#include "TVector3.h"

namespace fhicl{
    class ParameterSet;
}

namespace gar{
    namespace rec{
        namespace alg{

            class StripSplitterAlg {
            public:

                StripSplitterAlg(fhicl::ParameterSet const& pset);

                virtual ~StripSplitterAlg();

                void reconfigure(fhicl::ParameterSet const& pset);

                void ClearLists();

                void PrepareAlgo(const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

                void DoStripSplitting();

                std::vector <gar::rec::CaloHit*> getSplitHits() { return splitStripHits; }

                std::vector <gar::rec::CaloHit*> getUnSplitHits() { return unSplitStripHits; }

            protected:

                std::vector <gar::rec::CaloHit*> getVirtualHits(gar::rec::CaloHit* hit);

            private:

                gar::geo::GeometryCore const* fGeo; ///< geometry information
                TGeoManager* fGeoManager;

                std::string fSSAAlgName;
                unsigned int m_Verbose;
                int fInnerSymmetry;
                double fStripWidth;
                double fStripLength;
                int fnVirtual;

                std::vector<gar::rec::CaloHit*> m_CaloHitVecOdd;
                std::vector<gar::rec::CaloHit*> m_CaloHitVecEven;

                std::vector <gar::rec::CaloHit*> unSplitStripHits;
                std::vector <gar::rec::CaloHit*> splitStripHits;
            };

        } // namespace alg
    } // namespace rec
} //namespace gar

#endif /* GAR_RECOALG_KNNClusterAlg_h */
