#include <string>
#include <algorithm>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "Geometry/GeometryGAr.h"
#include "MCCheater/BackTracker.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "nugen/EventGeneratorBase/GENIE/EVGBAssociationUtil.h"

#include "RecoAlg/ClusterShapes.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace gar {
    namespace rec {

        class eveLoc {
        public:
            eveLoc(int id)
            : eveID(id)
            {}

            ~eveLoc() {}

            friend bool operator < (eveLoc const& a, eveLoc const& b) {
                return a.eveID < b.eveID;
            }

            int GetEveID() const { return eveID; }

        private:
            int                    eveID;
        };

        class CaloClusterCheater : public art::EDProducer {
        public:
            explicit CaloClusterCheater(fhicl::ParameterSet const& pset);

            // Plugins should not be copied or assigned.
            CaloClusterCheater(CaloClusterCheater const &) = delete;
            CaloClusterCheater(CaloClusterCheater &&) = delete;
            CaloClusterCheater & operator = (CaloClusterCheater const &) = delete;
            CaloClusterCheater & operator = (CaloClusterCheater &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        private:

            std::string  fCaloHitModuleLabel; ///< label for module creating rec::CaloHit objects
            std::string  fInstanceName; ///< product instance name

            unsigned int fMinHits;            ///< minimum number of hits to make a cluster
            bool fIgnoreNeutrons;             ///< flag to ignore creating neutron clusters

            gar::geo::GeometryCore const* fGeo;        ///< geometry information

            void CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);
            gar::rec::Cluster CreateCluster(const simb::MCParticle *part, std::vector< art::Ptr<gar::rec::CaloHit> > hitVec);
        };

        //-----------------------------------------------------------------------------------

        CaloClusterCheater::CaloClusterCheater(fhicl::ParameterSet const& pset) :
            EDProducer{pset}
            {
                fCaloHitModuleLabel    = pset.get< std::string  >("CaloHitModuleLabel", "calohit"  );
                fInstanceName          =  pset.get<std::string >("InstanceLabelName", "");

                fMinHits               = pset.get< unsigned int >("MinHits",            1          );
                fIgnoreNeutrons        = pset.get< bool >("IgnoreNeutrons",             true      );

                fGeo = gar::providerFrom<geo::GeometryGAr>();

                art::InputTag tag(fCaloHitModuleLabel, fInstanceName);
                consumes< std::vector<gar::rec::CaloHit> >(tag);
                produces< std::vector<gar::rec::Cluster> >(fInstanceName);
                produces< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> >(fInstanceName);
            }

            //-----------------------------------------------------------------------------------
            void CaloClusterCheater::produce(art::Event & e)
            {
                cheat::BackTrackerCore const* bt = gar::providerFrom<cheat::BackTracker>();

                //Collect the hits to be passed to the algo
                std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
                this->CollectHits(e, fCaloHitModuleLabel, fInstanceName, artHits);

                std::map< eveLoc, std::vector< art::Ptr<gar::rec::CaloHit> > > eveCaloHitMap;

                // loop over all hits and fill in the map
                for( auto const& itr : artHits ) {

                    std::vector<gar::cheat::CalIDE> eveides = bt->CaloHitToCalIDEs(itr);
                    // loop over all eveides for this hit
                    for(size_t ieve = 0; ieve < eveides.size(); ieve++) {
                        // don't worry about eve particles that contribute less than 10% of the
                        // energy in the current hit
                        if( eveides[ieve].energyFrac < 0.1) continue;

                        eveLoc el(eveides[ieve].trackID);
                        eveCaloHitMap[el].push_back(itr);

                    } // end loop over eve IDs for this hit
                }// end loop over hits

                // loop over the map and make clusters
                std::unique_ptr< std::vector<gar::rec::Cluster> > clustercol(new std::vector<gar::rec::Cluster>);
                std::unique_ptr< art::Assns<gar::rec::Cluster, gar::rec::CaloHit> > assn(new art::Assns<gar::rec::Cluster, gar::rec::CaloHit>);

                for(auto hitMapItr : eveCaloHitMap) {

                    MF_LOG_DEBUG("CaloClusterCheater")
                    << "Found cluster of " << hitMapItr.second.size() << " hits"
                    << " with eveid " << hitMapItr.first.GetEveID();

                    if(hitMapItr.second.size() < fMinHits) continue;

                    const simb::MCParticle *part = bt->TrackIDToParticle(hitMapItr.first.GetEveID());
                    if(!part) {
                        MF_LOG_WARNING("CaloClusterCheater")
                        << "Cannot find MCParticle for eveid: " << hitMapItr.first.GetEveID();
                        continue;
                    }

                    if(fIgnoreNeutrons && part->PdgCode() == 2112) continue;

                    gar::rec::Cluster cluster = this->CreateCluster(part, hitMapItr.second);

                    clustercol->emplace_back(cluster);
                    // association the hits to this cluster
                    evgb::util::CreateAssn(*this, e, *clustercol, hitMapItr.second, *assn, clustercol->size() - 1);
                    MF_LOG_DEBUG("CaloClusterCheater") << "adding cluster: \n" << clustercol->back() << "\nto collection.";

                } // end loop over the map

                e.put(std::move(clustercol), fInstanceName);
                e.put(std::move(assn), fInstanceName);
            }

            //-----------------------------------------------------------------------------------
            void CaloClusterCheater::CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
            {
		art::InputTag itag(label,instance);
		auto theHits = evt.getHandle< std::vector<gar::rec::CaloHit> >(itag);
                if (!theHits) return;

                for (unsigned int i = 0; i < theHits->size(); ++i)
                {
                    const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                    hitVector.push_back(hit);
                }
            }

            //-----------------------------------------------------------------------------------
            gar::rec::Cluster CaloClusterCheater::CreateCluster(const simb::MCParticle *part, std::vector< art::Ptr<gar::rec::CaloHit> > hitVec)
            {
                gar::rec::Cluster clu;

                unsigned n = hitVec.size() ;
                unsigned i = 0 ;
                float rmin = 9999.;
                float rmax = 0.;
                std::vector<float> tmin;
                std::vector<float> tmax;

                std::vector<float> a, t, x, y, z;
                a.resize(n);
                t.resize(n);
                x.resize(n);
                y.resize(n);
                z.resize(n);

                for( auto const &it : hitVec ) {
                    a[i] = it->Energy() ;
                    t[i] = it->Time().first ;
                    x[i] = it->Position()[0] ;
                    y[i] = it->Position()[1] ;
                    z[i] = it->Position()[2] ;

                    // clu.addHit( it , a[i] ) ;

                    //get time of first layer and last layer to compute the difference
                    float r = std::sqrt(y[i]*y[i] + z[i]*z[i]);
                    if( std::fabs(x[i]) < fGeo->GetECALEndcapStartX() ){
                        //case in barrel
                        if(rmin >= r)
                        {
                            rmin = r;
                            tmin.push_back(t[i]);
                        }
                        if(rmax <= r)
                        {
                            rmax = r;
                            tmax.push_back(t[i]);
                        }
                    }
                    else{
                        //case endcap
                        if( rmin >= std::fabs(x[i]) )
                        {
                            rmin = std::fabs(x[i]);
                            tmin.push_back(t[i]);
                        }
                        if( rmax <= std::fabs(x[i]) )
                        {
                            rmax = std::fabs(x[i]);
                            tmax.push_back(t[i]);
                        }
                    }

                    ++i ;
                }

                //get minimum time of the first layer
                auto time_first = std::min_element( tmin.begin(), tmin.end() );
                assert(time_first != tmin.end());
                //get minimum time of the last layer
                auto time_last = std::min_element( tmax.begin(), tmax.end() );
                assert(time_last != tmax.end());

                util::ClusterShapes cs(n, a, t, x, y, z) ;

                clu.setEnergy( cs.getTotalAmplitude()  ) ;
                clu.setTime( cs.getAverageTime(), *time_first - *time_last ) ;
                clu.setPosition( cs.getCenterOfGravity() ) ;

                // direction of cluster's PCA
                float* d = cs.getEigenVecInertia() ;
                clu.setEigenVectors( d );

                //Main eigenvector
                CLHEP::Hep3Vector v( d[0], d[1], d[2] ) ;

                clu.setITheta( v.theta() )  ;
                clu.setIPhi( v.phi() ) ;

                float param[6] ;

                // param[0] = cs.getElipsoid_r1() ;
                param[0] = cs.getElipsoid_r_forw();
                param[1] = cs.getElipsoid_r_back();
                param[2] = cs.getElipsoid_r2() ;
                param[3] = cs.getElipsoid_r3() ;
                param[4] = cs.getElipsoid_vol() ;
                param[5] = cs.getWidth() ;

                clu.setShape( param ) ;

                clu.setParticleID( part->PdgCode() );

                return clu;
            }

        } //namespace rec
    } //namespace gar

    namespace gar {
        namespace rec {
            DEFINE_ART_MODULE(CaloClusterCheater)
        }
    }
