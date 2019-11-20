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
#include "RecoAlg/TrackPropagator.h"

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
            void beginJob() override;
            void produce(art::Event & e) override;

        private:

            // Declare member data here.
            std::string fTrackLabel;    ///< label to find the reco tracks
            std::string fClusterLabel;  ///< label to find the right reco caloclusters
            int fVerbosity;
            float fExtrapAngleCut;      ///< Cut on extrapolation angle
            float fExtrapDistCut;       ///< Cut on extrapolation distance
            float fMatchHelixCut;       ///< Cut on helix-cluster distance

            const detinfo::DetectorProperties*  fDetProp;    ///< detector properties
            const geo::GeometryCore*            fGeo;        ///< pointer to the geometry

            TH2D* withWithoutDistS;                          ///< Extrap angle vs extrap distance
            TH2D* withWithoutDistL;
            TH2D* distOffvsExtrapS;                          ///< Histo distance off helix vs extrap distance
            TH2D* distOffvsExtrapL;
        };



        TPCECALAssociation::TPCECALAssociation(fhicl::ParameterSet const & p) {

            fTrackLabel     = p.get<std::string>("TrackLabel", "track");
            fClusterLabel   = p.get<std::string>("ClusterLabel","calocluster");
            fVerbosity      = p.get<int>("Verbosity", 0);
            fExtrapAngleCut = p.get<float>("ExtrapAngleCut",  0.8);    // acos(angle) cut
            fExtrapDistCut  = p.get<float>("ExtrapDistCut", 120.0);    // in cm
            fMatchHelixCut  = p.get<float>("MatchHelixCut",  14.0);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            
            produces< art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd > >();
        }



        void TPCECALAssociation::beginJob() {

            if (fVerbosity>0) {
                art::ServiceHandle< art::TFileService > tfs;
                withWithoutDistS = tfs->make<TH2D>("withWithoutDistS",
                    "Extrapolation angle vs extrapolation distance", 100,0.0,200.0, 100, 0.0, 1.0);
                withWithoutDistL = tfs->make<TH2D>("withWithoutDistL",
                    "Extrapolation angle vs extrapolation distance", 100,0.0,500.0, 100, 0.0, 1.0);
                distOffvsExtrapS = tfs->make<TH2D>("distOffvsExtrapS", 
                    "Helix-cluster dist vs extrapolation distance", 100,0.0,200.0, 80, 0.0,80.0);
                distOffvsExtrapL = tfs->make<TH2D>("distOffvsExtrapL", 
                    "Helix-cluster dist vs extrapolation distance", 100,0.0,500.0, 100, 0.0,200.0);
            }
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



            std::unique_ptr<art::Assns<gar::rec::Cluster, gar::rec::Track, gar::rec::TrackEnd>>
                            ClusterTrackAssns(new art::Assns<gar::rec::Cluster,gar::rec::Track, gar::rec::TrackEnd>);
            auto const clusterPtrMaker = art::PtrMaker<rec::Cluster>(e, theClusterHandle.id());
            auto const trackPtrMaker   = art::PtrMaker<rec::Track>  (e, theTrackHandle.id());

            for (size_t iCluster=0; iCluster<theClusterHandle->size(); ++iCluster) {
                gar::rec::Cluster cluster = (*theClusterHandle)[iCluster];
                for (size_t iTrack=0; iTrack<theTrackHandle->size(); ++iTrack) {
                    gar::rec::Track track = (*theTrackHandle)[iTrack];

                    // Which trackend to use for the extrapolation?
                    float begDir[3];    float endDir[3];
                    for (int i=0; i<3; ++i) {
                        begDir[i] = -track.VtxDir()[i];        // Unit vectors continuing the track
                        endDir[i] = -track.EndDir()[i];
                    }

                    // Dot product of direction vector with distance from trackend to cluster
                    TVector3 begDelPoint = TVector3(cluster.Position())
                                          -TVector3(track.Vertex());
                    TVector3 begDelPointHat = begDelPoint.Unit();
                    TVector3 endDelPoint = TVector3(cluster.Position())
                                          -TVector3(track.End());
                    TVector3 endDelPointHat = endDelPoint.Unit();
                    float begDot    = begDelPoint.Dot( TVector3(begDir) );
                    float endDot    = endDelPoint.Dot( TVector3(endDir) );
                    float begDotHat = begDelPointHat.Dot( TVector3(begDir) );
                    float endDotHat = endDelPointHat.Dot( TVector3(endDir) );



                    // Plot extrapolation angle vs distance before cut on anything
                    if (fVerbosity>0) {
                        withWithoutDistS->Fill(begDot,begDotHat);
                        withWithoutDistS->Fill(endDot,endDotHat);
                        withWithoutDistL->Fill(begDot,begDotHat);
                        withWithoutDistL->Fill(endDot,endDotHat);
                    }

                    // Consider only extrapolations that pass the angle cut; neither does, we done.
                    bool considerBeg = begDotHat>fExtrapAngleCut;
                    bool considerEnd = endDotHat>fExtrapAngleCut;
                    if (!considerBeg &&!considerEnd ) continue;

                    TrackEnd whichEnd;
                    if        ( considerBeg &&!considerEnd ) {
                        whichEnd = TrackEndBeg;
                    } else if (!considerBeg && considerEnd ) {
                        whichEnd = TrackEndEnd;
                    } else if ( begDot<endDot ) {
                        // If both ends inside track angle cut, pick the shortest extrapolation distance
                        whichEnd = TrackEndBeg;
                    } else {
                        whichEnd = TrackEndEnd;
                    }



                    // Find the distance off the track, and the throw arm of the extrapolation
                    float matchHelixDist, extrapDist;
                    if ( whichEnd==TrackEndBeg ) {
                        util::TrackPropagator::DistXYZ(track.TrackParBeg(),track.Vertex(), cluster.Position(), matchHelixDist);
                        extrapDist = begDot;
                    } else {
                        util::TrackPropagator::DistXYZ(track.TrackParEnd(),track.End(),    cluster.Position(), matchHelixDist);
                        extrapDist = endDot;
                    }

                    // Histogram distance off track vs extrapolation distance and cut on it
                    if (fVerbosity>0) {
                       distOffvsExtrapS->Fill(extrapDist, matchHelixDist);
                       distOffvsExtrapL->Fill(extrapDist, matchHelixDist);
                    }
                    
                    if ( extrapDist<fExtrapDistCut && matchHelixDist<fMatchHelixCut ) {
                        // Make the cluster-track association
                        art::Ptr<gar::rec::Cluster> const clusterPtr = clusterPtrMaker(iCluster);
                        art::Ptr<gar::rec::Track>   const   trackPtr = trackPtrMaker(iTrack);
                        ClusterTrackAssns->addSingle(clusterPtr,trackPtr,whichEnd);
                    }
                }
            }

            e.put(std::move(ClusterTrackAssns));
            return;
        }



        DEFINE_ART_MODULE(TPCECALAssociation)

    } // namespace rec
} // namespace gar
