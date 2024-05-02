////////////////////////////////////////////////////////////////////////
// Class:       TPCECALAssociation2
// Plugin Type: producer (art v2_11_02)
// File:        TPCECALAssociation2_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "MCCheater/BackTracker.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "RecoAlg/TrackPropagator.h"

#include "art_root_io/TFileService.h"

#include <Math/RootFinder.h>
#include <Math/WrappedFunction.h>
#include <Math/Integrator.h>
#include <functional>
#include <unordered_map>

namespace gar {
    namespace rec {

        class TrackAssn {
        public:
            TrackAssn(int id);
            virtual ~TrackAssn();

            void set_chisq_data(gar::rec::TrackEnd iEnd, int clusid, float chisq, float chisq_cut);
            std::vector<std::pair<int, float>> get_chisq_data(gar::rec::TrackEnd iEnd) const;
            void set_direction(gar::rec::TrackEnd iEnd);
            void check_direction();
            std::vector<std::pair<int, float>> get_chisq_candidate() const;
            int get_nchisq_candidate() const;
            gar::rec::TrackEnd get_end_candidate() const;

        private:
            int trkid;
            gar::rec::TrackEnd candidate_end;
            std::vector<std::pair<int, float>> chisq_fwd;
            int n_chisq_fwd=0;
            float chisq_fwd_max=-1.0;
            std::vector<std::pair<int, float>> chisq_bak;
            int n_chisq_bak=0;
            float chisq_bak_max=-1.0;

        };

        TrackAssn::TrackAssn(int id)
        {
            this->trkid = id;
            return;
        }

        TrackAssn::~TrackAssn()
        {
            return;
        }

        void TrackAssn::set_chisq_data(gar::rec::TrackEnd iEnd, int clusid, float chisq, float chisq_cut) {
            if ((chisq>0)&&(chisq<=chisq_cut)) {
                std::pair<int, float> chisq_data_point = std::make_pair(clusid, chisq);
                if (iEnd==TrackEndBeg) {
                    chisq_fwd.push_back(chisq_data_point);
                    n_chisq_fwd++;
                    if(chisq > chisq_fwd_max) chisq_fwd_max = chisq;
                } else {
                    chisq_bak.push_back(chisq_data_point);
                    n_chisq_bak++;
                    if(chisq > chisq_bak_max) chisq_bak_max = chisq;
                }
            }
        }

        std::vector<std::pair<int, float>> TrackAssn::get_chisq_data(gar::rec::TrackEnd iEnd) const {
            if (iEnd==TrackEndBeg) {
                return chisq_fwd;
            } else {
                return chisq_bak;
            } 
        }

        void TrackAssn::set_direction(gar::rec::TrackEnd iEnd) {
            // Set association direction by hand
            candidate_end = iEnd;
        }

        void TrackAssn::check_direction() {
            // Check what direction has the most associations and select it
            if (n_chisq_fwd > n_chisq_bak) {
                candidate_end = TrackEndBeg;
            } else if (n_chisq_fwd < n_chisq_bak) {
                candidate_end = TrackEndEnd;
            } else {
                // What happens if both FWD and BAK have the same number of candidates?
                // Keep the one with lowest maximum value
                if (chisq_fwd_max > chisq_bak_max) {
                    candidate_end = TrackEndEnd;
                } else {
                    candidate_end = TrackEndBeg;
                }
            }
        }

        std::vector<std::pair<int, float>> TrackAssn::get_chisq_candidate() const {
            if (candidate_end==TrackEndBeg) {
                return chisq_fwd;
            } else {
                return chisq_bak;
            } 
        }

        int TrackAssn::get_nchisq_candidate() const {
            if (candidate_end==TrackEndBeg) {
                return n_chisq_fwd;
            } else {
                return n_chisq_bak;
            } 
        }

        gar::rec::TrackEnd TrackAssn::get_end_candidate() const {
            return candidate_end;
        }

        struct Association {
            int TrkID;
            gar::rec::TrackEnd TrkEnd;
            int ClusID;
            float Chisq;
        };

        class TPCECALAssociation2 : public art::EDProducer {
        public:
            explicit TPCECALAssociation2(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            TPCECALAssociation2(TPCECALAssociation2 const &) = delete;
            TPCECALAssociation2(TPCECALAssociation2 &&) = delete;
            TPCECALAssociation2 & operator = (TPCECALAssociation2 const &) = delete;
            TPCECALAssociation2 & operator = (TPCECALAssociation2 &&) = delete;

            // Required functions.
            void beginJob() override;
            void produce(art::Event & e) override;

            // Auxiliary functions
            float get_global_ecal_t0(art::Handle<std::vector<gar::rec::Cluster>> ClusterHandle);
            float helix_circle_intersections(float phi, float y0, float z0, float R, float phi0, float r);
            void helix(float phi, float *trackPar, float *trackEnd, float *projected, float t0 = 0.0);

        private:

            // Declare member data here.
            std::string fTrackLabel;    ///< label to find the reco tracks
            std::string fClusterLabel;  ///< label to find the right reco caloclusters
            std::string fInstanceName;
            int fVerbosity;
            bool fGlobalTimeCorrection;
            bool fClusterTimeCorrection;
            bool fCheatDirection;

            float fDriftVelocity;

            float fMaxAngle;   ///< max angle for propagation in units of pi
            float fChisqCut;   ///< 
            
            cheat::BackTrackerCore* BackTrack;

            const geo::GeometryCore*            fGeo;        ///< pointer to the geometry

            // Position of TPC from geometry service; 1 S Boston Ave.
            float ItsInTulsa[3];

        };



        TPCECALAssociation2::TPCECALAssociation2(fhicl::ParameterSet const & p) : EDProducer{p} {

            fTrackLabel            = p.get<std::string>("TrackLabel", "track");
            fClusterLabel          = p.get<std::string>("ClusterLabel","calocluster");
            fInstanceName          = p.get<std::string>("InstanceName","");
            fVerbosity             = p.get<int>("Verbosity", 0);
            fGlobalTimeCorrection  = p.get<bool>("GlobalTimeCorrection", false);
            fClusterTimeCorrection = p.get<bool>("ClusterTimeCorrection", false);
            fCheatDirection        = p.get<bool>("CheatDirection", false);

            fMaxAngle          = p.get<float>("MaxAngle", 1.0);
            fChisqCut          = p.get<float>("ChisqCut", 1e3);

            fGeo     = gar::providerFrom<geo::GeometryGAr>();
        
            auto detProp   = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fDriftVelocity = detProp->DriftVelocity(detProp->Efield(), detProp->Temperature());

            produces< art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd > >(fInstanceName);
        }

        void TPCECALAssociation2::beginJob() {

            ItsInTulsa[0] = fGeo->TPCXCent();
            ItsInTulsa[1] = fGeo->TPCYCent();
            ItsInTulsa[2] = fGeo->TPCZCent();

            if (fGlobalTimeCorrection && fClusterTimeCorrection) {
                throw cet::exception("TPCECALAssociation2") << " Job started with conflicting configuration!" << std::endl;
            }

        }

        void TPCECALAssociation2::produce(art::Event & e) {
            gErrorIgnoreLevel = kInfo;

            if (fCheatDirection) {
                cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
                BackTrack = const_cast<cheat::BackTrackerCore*>(const_bt);
            }
        
            // Get tracks and clusters.  If either is missing, just skip this event
            // processing.  That's not an exception

            auto TrackHandle = e.getHandle< std::vector<gar::rec::Track> >(fTrackLabel);
            if (!TrackHandle) return;

            art::InputTag itag(fClusterLabel, fInstanceName);
            auto ClusterHandle = e.getHandle< std::vector<gar::rec::Cluster> >(itag);
            // Return in case there's no clusters
            if (!ClusterHandle) return;

            // Here are the associations
            std::unique_ptr<art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd>> ClusterTrackAssns(new art::Assns<gar::rec::Cluster,gar::rec::Track, gar::rec::TrackEnd>);
            auto const clusterPtrMaker = art::PtrMaker<rec::Cluster>(e, ClusterHandle.id());
            auto const   trackPtrMaker = art::PtrMaker<rec::Track>  (e, TrackHandle.id());

            std::unordered_map<int, Association> ClusterToAssociation;

            float ecal_t0 = 0.0;
            if (fGlobalTimeCorrection && (ClusterHandle->size() > 0)) {
                ecal_t0 = get_global_ecal_t0(ClusterHandle);
            }

            for (size_t iTrack=0; iTrack<TrackHandle->size(); ++iTrack) {
                gar::rec::Track track = (*TrackHandle)[iTrack];
                TrackAssn trackassn(iTrack);

                for (gar::rec::TrackEnd iEnd = TrackEndBeg; iEnd >= TrackEndEnd; --iEnd) {
                    float trackPar[5];
                    float trackEnd[3];
                    if (iEnd==TrackEndBeg) {
                        for (int i=0; i<5; ++i) trackPar[i] = track.TrackParBeg()[i];
                        for (int i=0; i<3; ++i) trackEnd[i] = track.Vertex()[i];
                    } else {
                        for (int i=0; i<5; ++i) trackPar[i] = track.TrackParEnd()[i];
                        for (int i=0; i<3; ++i) trackEnd[i] = track.End()[i];
                    }

                    float phi_max = -1.;
                    if(trackPar[2] > 0) {
                        phi_max = trackPar[3]+fMaxAngle*TMath::Pi();
                    } else {
                        phi_max = trackPar[3]-fMaxAngle*TMath::Pi();
                    }

                    for (size_t iCluster=0; iCluster<ClusterHandle->size(); ++iCluster) {
                        gar::rec::Cluster cluster = (*ClusterHandle)[iCluster];

                        TVector3 clusterCenter(cluster.Position());
                        float yClus = clusterCenter[1];
                        float zClus = clusterCenter[2];
                        float rClus = std::hypot(zClus-ItsInTulsa[2],yClus-ItsInTulsa[1]);
                        float tClus = cluster.Time();

                        auto helix_circle_intersections_to_wrap = [this, trackPar, rClus](float phi) {
                            const float y0   = trackPar[0];
                            const float z0   = trackPar[1];
                            const float R    = 1/trackPar[2];
                            const float phi0 = trackPar[3];
                            const float r    = rClus;
                            return this->helix_circle_intersections(phi, y0, z0, R, phi0, r);
                        };

                        ROOT::Math::WrappedFunction<std::function<float(float)>> helix_circle_intersections_wrapped(helix_circle_intersections_to_wrap);

                        ROOT::Math::RootFinder rootFinder(ROOT::Math::RootFinder::kBRENT);
                        rootFinder.SetFunction(helix_circle_intersections_wrapped, trackPar[3], phi_max);
                        float chisq = 0.;
                        try {
                            rootFinder.Solve();
                            float root = rootFinder.Root();

                            float projected[3];
                            if (fGlobalTimeCorrection) {
                                helix(root, trackPar, trackEnd, projected, ecal_t0);
                            } else if (fClusterTimeCorrection) {
                                helix(root, trackPar, trackEnd, projected, tClus);
                            } else {
                                helix(root, trackPar, trackEnd, projected);
                            }

                            for(size_t i=0; i<3; ++i){
                                chisq += TMath::Power(projected[i]-clusterCenter[i], 2);
                            }
                            chisq = chisq/3.;

                        } catch (const std::exception &excpt) {
                            chisq = -1.;

                        }

                        trackassn.set_chisq_data(iEnd, iCluster, chisq, fChisqCut); // this will fill (or not, depending on the cut value) the appropriate chisq data vector (either FWD or BAK)

                    } // end loop clusters

                } // end loop over 2 ends of track

                if (fCheatDirection) {
                    // Use the backtracker to get the corresponding MCParticle
                    // and use the reco end that is closer to the true start
                    std::vector<std::pair<simb::MCParticle*,float>> trakt;
                    trakt = BackTrack->TrackToMCParticles( const_cast<rec::Track*>(&track) );
                    if (trakt.size() > 0) {
                        simb::MCParticle mcp = *trakt[0].first;
                        const TVector3& true_start_pos = mcp.Position(0).Vect();

                        TVector3 end_begin_pos(track.Vertex()[0], track.Vertex()[1], track.Vertex()[2]);
                        TVector3 end_end_pos(track.End()[0], track.End()[1], track.End()[2]);

                        float end_begin_dist = (true_start_pos-end_begin_pos).Mag();
                        float end_end_dist   = (true_start_pos-end_end_pos).Mag();

                        if (end_begin_dist >= end_end_dist) {
                            trackassn.set_direction(TrackEndBeg);
                        } else {
                            trackassn.set_direction(TrackEndEnd);
                        }

                    } else {
                        // If the backtracker fails use the default solution
                        trackassn.check_direction();
                    }
                } else {
                    // Use the track end with more clusters passing the cut
                    trackassn.check_direction();
                }

                int n_candidate = trackassn.get_nchisq_candidate();
                gar::rec::TrackEnd end_candidate = trackassn.get_end_candidate();
                std::vector<std::pair<int, float>> chisq_candidate = trackassn.get_chisq_candidate();

        
                for (size_t i=0; i<static_cast<size_t>(n_candidate); ++i) {
                    Association association;
                    association.TrkID = iTrack;
                    association.TrkEnd = end_candidate;
                    association.ClusID = chisq_candidate[i].first;
                    association.Chisq = chisq_candidate[i].second;

                    if (ClusterToAssociation.find(association.ClusID) == ClusterToAssociation.end()) {
                        // If there is no previous candidate fill the map
                        ClusterToAssociation[association.ClusID] = association;
                    } else {
                        // If there is a previous candidate check chisq
                        float previous_chisq = ClusterToAssociation[association.ClusID].Chisq;
                        if (association.Chisq < previous_chisq) {
                            // New is smaller, fill this
                            ClusterToAssociation[association.ClusID] = association;
                        }
                    }
                    
                } // end loop over candidate clusters

            } // end loop over tracks

            for (auto& [key, value]: ClusterToAssociation) {
                art::Ptr<gar::rec::Track> const trackPtr = trackPtrMaker(value.TrkID);
                art::Ptr<gar::rec::Cluster> const clusterPtr = clusterPtrMaker(value.ClusID);
                ClusterTrackAssns->addSingle(clusterPtr,trackPtr,value.TrkEnd);
            }

            e.put(std::move(ClusterTrackAssns), fInstanceName);
            return;
        }

        float TPCECALAssociation2::get_global_ecal_t0(art::Handle<std::vector<gar::rec::Cluster>> ClusterHandle) {
            float ecal_t0;

            std::vector<float> ClusterTimeVector;
            for (auto const& cluster: *ClusterHandle) {
                ClusterTimeVector.push_back(cluster.Time());
            }

            ecal_t0 = *std::min_element(ClusterTimeVector.begin(), ClusterTimeVector.end());

            return ecal_t0;
        }

        float TPCECALAssociation2::helix_circle_intersections(float phi, float y0, float z0, float R, float phi0, float r) {
                return TMath::Power(y0-R*(TMath::Cos(phi)-TMath::Cos(phi0))-ItsInTulsa[1], 2)+TMath::Power(z0+R*(TMath::Sin(phi)-TMath::Sin(phi0))-ItsInTulsa[2], 2)-TMath::Power(r, 2);
            }

        void TPCECALAssociation2::helix(float phi, float *trackPar, float *trackEnd, float *projected, float t0) {
            
            projected[0] = trackEnd[0] + ((trackEnd[0] > 0.) ? 1. : ((trackEnd[0] < 0.) ? -1. : 0.))*fDriftVelocity*t0 + (1/trackPar[2])*TMath::Tan(trackPar[4])*(phi-trackPar[3]);
            projected[1] = trackPar[0] - (1/trackPar[2])*(TMath::Cos(phi) - TMath::Cos(trackPar[3]));
            projected[2] = trackPar[1] + (1/trackPar[2])*(TMath::Sin(phi) - TMath::Sin(trackPar[3]));
        }

        DEFINE_ART_MODULE(TPCECALAssociation2)

    } // namespace rec
} // namespace gar
