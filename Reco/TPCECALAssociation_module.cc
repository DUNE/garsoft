////////////////////////////////////////////////////////////////////////
// Class:       TPCECALAssociation
// Plugin Type: producer (art v2_11_02)
// File:        TPCECALAssociation_module.cc
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "Utilities/TrackPropagator.h"

#include "art/Framework/Services/Optional/TFileService.h"
#include "TH2D.h"



namespace gar {
    namespace rec {



        class TPCECALAssociation : public art::EDProducer {
        public:
            explicit TPCECALAssociation(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            TPCECALAssociation(TPCECALAssociation const &) = delete;
            TPCECALAssociation(TPCECALAssociation &&) = delete;
            TPCECALAssociation & operator = (TPCECALAssociation const &) = delete;
            TPCECALAssociation & operator = (TPCECALAssociation &&) = delete;

            // Required functions.
            void beginJob();
            void produce(art::Event & e) override;

        private:

            // Declare member data here.
            std::string fTrackLabel;  ///< label to find the reco tracks
            std::string fClusterLabel;  ///< label to find the right reco caloclusters
            int fVerbosity;

            const detinfo::DetectorProperties*  fDetProp;    ///< detector properties
            const geo::GeometryCore*            fGeo;        ///< pointer to the geometry

            TH2D* distOffvsExtrap;                            ///< Histo distance off helix vs extrap distance
        };



        TPCECALAssociation::TPCECALAssociation(fhicl::ParameterSet const & p) {

            fTrackLabel   = p.get<std::string>("TrackLabel", "track");
            fClusterLabel = p.get<std::string>("ClusterLabel","calocluster");
            fVerbosity    = p.get<int>("Verbosity", 0);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            
            produces< art::Assns<gar::rec::Cluster, gar::rec::Track> >();
        }



        void TPCECALAssociation::beginJob() {

            art::ServiceHandle< art::TFileService > tfs;
            distOffvsExtrap = tfs->make<TH2D>("distOffvsExtrap", 
                "Helix-cluster dist vs extrapolation distance", 100,0.0,100.0, 100, 0.0,20.0);
        
        }



        void TPCECALAssociation::produce(art::Event & e) {
        
            // Get tracks and clusters.  If either is missing, just skip this event
            // processing.  That's not an exception
            art::Handle< std::vector<gar::rec::Cluster> > theClusterHandle;
            e.getByLabel(fClusterLabel, theClusterHandle);
            if (!theClusterHandle.isValid()) return;
            art::Handle< std::vector<gar::rec::Track> >   theTrackHandle;
            e.getByLabel(fTrackLabel, theTrackHandle);
            if (!theTrackHandle.isValid()) return;



            std::unique_ptr<art::Assns<gar::rec::Cluster,gar::rec::Track>> ClusterTrackAssns(new art::Assns<gar::rec::Cluster,gar::rec::Track>);
            auto const clusterPtrMaker = art::PtrMaker<rec::Cluster>(e, theClusterHandle.id());
            auto const trackPtrMaker   = art::PtrMaker<rec::Track>  (e, theTrackHandle.id());

            for (size_t iCluster=0; iCluster<theClusterHandle->size(); ++iCluster) {
                gar::rec::Cluster cluster = (*theClusterHandle)[iCluster];
                for (size_t iTrack=0; iTrack<theTrackHandle->size(); ++iTrack) {
                    gar::rec::Track track = (*theTrackHandle)[iTrack];

                    // Which trackend to use for the extrapolation?
                    float begDir[3];    float endDir[3];
                    for (int i=0; i<3; ++i) {
                        begDir[i] = +track.VtxDir()[i];
                        endDir[i] = -track.EndDir()[i];
                    }

                    // Dot product of direction vector with distance from trackend to cluster
                    TVector3 begDelPoint = TVector3(cluster.Position())
                                          -TVector3(track.Vertex());
                    TVector3 endDelPoint = TVector3(cluster.Position())
                                          -TVector3(track.End());
                    float begDot = begDelPoint.Dot( TVector3(begDir) );
                    float endDot = endDelPoint.Dot( TVector3(endDir) );

                    // Want the smallest positive number
                    TrackEnd which;
                    if ( begDot<0 && endDot<0 ) continue;
                    if ( begDot<0 ) {
                        which = TrackEndEnd;
                    } else if ( endDot>0 ) {
                        which = TrackEndBeg;
                    } else if ( begDot<endDot ) {
                        which = TrackEndBeg;
                    } else {
                        which = TrackEndEnd;
                    }

                    // Find the distance off the track, and the throw arm of the 
                    // extrapolation
                    float offTrackDist, throwArmExtrapol;
                    if ( which==TrackEndBeg ) {
                        util::TrackPropagator::DistXYZ(track.TrackParBeg(),track.Vertex(), cluster.Position(), &offTrackDist);
                        throwArmExtrapol = begDot;
                    } else {
                        util::TrackPropagator::DistXYZ(track.TrackParEnd(),track.End(),    cluster.Position(), &offTrackDist);
                        throwArmExtrapol = endDot;
                    }

                    

                    // For now just histogram the thang
                    distOffvsExtrap->Fill(throwArmExtrapol, offTrackDist);

                    
                    // Add the cluster-track association
                    auto const clusterPtr = clusterPtrMaker(iCluster);
                    auto const   trackPtr = trackPtrMaker(iTrack);
                    ClusterTrackAssns->addSingle(clusterPtr,trackPtr);
                }
            }

            e.put(std::move(ClusterTrackAssns));
            return;
        }



        DEFINE_ART_MODULE(TPCECALAssociation)

    } // namespace rec
} // namespace gar
