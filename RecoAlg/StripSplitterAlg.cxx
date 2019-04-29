//
//  StripSplitterAlg.cxx
//
//  Created by Eldwan Brianne on 18.04.2019
//

#include "fhiclcpp/ParameterSet.h"

#include "RecoAlg/StripSplitterAlg.h"

#include "TGeoNode.h"
#include "TGeoManager.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/LocalTransformation.h"

#include <algorithm>

namespace gar {
    namespace rec{
        namespace alg{

            //----------------------------------------------------------------------------
            StripSplitterAlg::StripSplitterAlg(fhicl::ParameterSet const& pset)
            {
                ClearLists();

                fGeo = gar::providerFrom<geo::Geometry>();
                fGeoManager = fGeo->ROOTGeoManager();
                this->reconfigure(pset);

                fInnerSymmetry = fGeo->GetECALInnerSymmetry();
                fStripWidth = fGeo->getStripWidth();
                fStripLength = -999;
                fnVirtual = 1;

                return;
            }

            //----------------------------------------------------------------------------
            StripSplitterAlg::~StripSplitterAlg()
            {
                return;
            }

            //----------------------------------------------------------------------------
            void StripSplitterAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                fSSAAlgName = pset.get<std::string>("SSAAlgName");
                //Algorithm parameters
                m_Verbose = pset.get<unsigned int>("Verbose", 0);

                return;
            }

            //----------------------------------------------------------------------------
            void StripSplitterAlg::ClearLists()
            {
                m_CaloHitVecOdd.clear();
                m_CaloHitVecEven.clear();

                unSplitStripHits.clear();
                splitStripHits.clear();

                return;
            }

            //----------------------------------------------------------------------------
            void StripSplitterAlg::PrepareAlgo(const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
            {
                std::cout << "StripSplitterAlg::PrepareAlgo()" << std::endl;

                //Clear the lists
                ClearLists();

                //Loop over all hits
                for (std::vector< art::Ptr<gar::rec::CaloHit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
                {
                    art::Ptr<gar::rec::CaloHit> hitPtr = *iter;
                    const gar::rec::CaloHit *hit = hitPtr.get();

                    if(fGeo->isTile(hit->CellID())) continue;

                    //check the layer
                    unsigned int layer = hit->GetLayer();
                    if(layer%2 == 0)
                    m_CaloHitVecEven.push_back( const_cast<gar::rec::CaloHit *>(hit) );
                    else{
                        m_CaloHitVecOdd.push_back( const_cast<gar::rec::CaloHit *>(hit) );
                    }
                }

                return;
            }

            //----------------------------------------------------------------------------
            void StripSplitterAlg::DoStripSplitting()
            {
                if(m_Verbose)
                std::cout << "StripSplitterAlg::DoStripSplitting()" << std::endl;

                return;
            }

            //----------------------------------------------------------------------------
            std::vector <gar::rec::CaloHit*> StripSplitterAlg::getVirtualHits(gar::rec::CaloHit* hit)
            {
                if(m_Verbose)
                std::cout << "StripSplitterAlg::getVirtualHits()" << std::endl;

                std::vector <gar::rec::CaloHit*> newhits;
                return newhits;
            }

        }// end namespace alg
    }// end namespace rec
}// end namespace gar
