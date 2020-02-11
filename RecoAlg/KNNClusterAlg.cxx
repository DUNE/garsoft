//
//  KNNClusterAlg.cxx
//
//  Created by Eldwan Brianne on 10.02.2019
//

#include "fhiclcpp/ParameterSet.h"

#include "RecoAlg/KNNClusterAlg.h"
#include "RecoAlg/NNClusters.h"

#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include <algorithm>

namespace gar {
    namespace rec{
        namespace alg{

            //----------------------------------------------------------------------------
            KNNClusterAlg::KNNClusterAlg(fhicl::ParameterSet const& pset)
            {
                fGeo = gar::providerFrom<geo::Geometry>();
                this->reconfigure(pset);
                ClearLists();
                return;
            }
            //----------------------------------------------------------------------------
            KNNClusterAlg::~KNNClusterAlg()
            {
                return;
            }
            //----------------------------------------------------------------------------
            void KNNClusterAlg::reconfigure(fhicl::ParameterSet const& pset)
            {
                fClusterAlgName = pset.get<std::string>("ClusterAlgName");

                //Algorithm parameters
                m_EnergyCut = pset.get<float>("EnergyCut", 0.);
                m_DistanceCut = pset.get<float>("DistanceCut", 5.f);
                m_Verbose = pset.get<unsigned int>("Verbose", 0);

                return;
            }
            //----------------------------------------------------------------------------
            void KNNClusterAlg::ClearLists()
            {
                clusterVector.clear();
                m_CaloHitVec.clear();
                m_TrackVec.clear();

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterAlg::PrepareAlgo(const std::vector< art::Ptr<gar::rec::Track> > &trkVector, const std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector, std::unordered_map< const gar::rec::Track*, art::Ptr<gar::rec::Track> > &trkMaptoArtPtr, std::unordered_map< const gar::rec::CaloHit*, art::Ptr<gar::rec::CaloHit> > &hitMaptoArtPtr)
            {
                //Clear the lists
                ClearLists();

                //Loop over all tracks
                for (std::vector< art::Ptr<gar::rec::Track> >::const_iterator iter = trkVector.begin(), iterEnd = trkVector.end(); iter != iterEnd; ++iter)
                {
                    art::Ptr<gar::rec::Track> trkPtr = *iter;
                    const gar::rec::Track *track = trkPtr.get();
                    m_TrackVec.push_back( const_cast<gar::rec::Track *>(track) );
                    //update the track art ptr map
                    trkMaptoArtPtr.insert( std::make_pair(track, trkPtr) );
                }
                //Loop over all hits
                for (std::vector< art::Ptr<gar::rec::CaloHit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
                {
                    art::Ptr<gar::rec::CaloHit> hitPtr = *iter;
                    const gar::rec::CaloHit *hit = hitPtr.get();
                    m_CaloHitVec.push_back( const_cast<gar::rec::CaloHit *>(hit) );
                    //update the hit art ptr map
                    hitMaptoArtPtr.insert( std::make_pair(hit, hitPtr) );
                }

                return;
            }

            //----------------------------------------------------------------------------
            void KNNClusterAlg::DoClustering()
            {
                GenericHitVec< gar::rec::CaloHit > h;

                GenericClusterVec< gar::rec::CaloHit > cl;

                EnergyCut< gar::rec::CaloHit > eCut( m_EnergyCut );

                XIndex< gar::rec::CaloHit, 280 > xIndex( -350 , 350. );

                NNDistance< gar::rec::CaloHit, float> dist( m_DistanceCut );

                AlgCluster< gar::rec::CaloHit > converter( fGeo->GetECALEndcapStartX() );

                int nHit = 0;
                nHit += m_CaloHitVec.size();

                // create a vector of generic hits from the collection applying an energy cut
                addToGenericHitVec( h , m_CaloHitVec , eCut , xIndex ) ;

                // cluster the hits with a nearest neighbour condition
                cluster( h.begin() , h.end() , std::back_inserter( cl )  , &dist ) ;

                if(m_Verbose > 0){
                    std::cout << "  passing " << h.size() << " of " << nHit
                    << "  hits to clustering (E_cut: " << m_EnergyCut << ") "
                    << "  found  " << cl.size() << " clusters " << std::endl ;
                }

                // create Clusters from the clustered GenericHits
                std::transform( cl.begin(), cl.end(), std::back_inserter( clusterVector ), converter ) ;

                return;
            }

        }// end namespace alg
    }// end namespace rec
}// end namespace gar
