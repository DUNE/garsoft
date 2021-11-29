////////////////////////////////////////////////////////////////////////////////
// Class:        MomentumPerformance
// Plugin Type: analyzer (art v2_11_02)
// File:        MomentumPerformance_module.cc
//
// Generated 20 Apr 2020 by Leo Bellantoni
//
// For each e, mu, pi+, p or antiparticle thereof from the primary vertex (but
// no cut on location of primary vertex!) use the backtracker to get a list of
// possible reconstructed tracks and pick the best matching trackend from that
// list.  Provide a little info about the matching and then provide info about
// the reconstructed track particularly but not only the momentum.
//
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <string>
#include <vector>
#include <unordered_map>



namespace gar {

    class MomentumPerformance : public art::EDAnalyzer {
    public:
        explicit MomentumPerformance(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        MomentumPerformance(MomentumPerformance const &) = delete;
        MomentumPerformance(MomentumPerformance &&) = delete;
        MomentumPerformance & operator = (MomentumPerformance const &) = delete;
        MomentumPerformance & operator = (MomentumPerformance &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;



    private:
        void ClearVectors();
        void FillVectors(art::Event const & e);

        void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
                                   float& forwardIonVal, float& backwardIonVal);
        float processOneDirection(std::vector<std::pair<float,float>> SigData,
                                  float ionizeTruncate);
        // Helper method for processOneDirection
        static bool lessThan_byE(std::pair<float,float> a, std::pair<float,float> b)
            {return a.first < b.first;}




        // Position of TPC from geometry service; 1 S Boston Ave.
        double ItsInTulsa[3];

        // Input data labels
        std::string fGeantLabel;
        std::string fTPCClusterLabel;
        std::string fTrackLabel;

        // Ionization truncation parameter
        float fIonizTruncate;

        // the analysis tree
        TTree *fTree;

        // Backtracker service
        cheat::BackTrackerCore* fBack;



        // The analysis tree. Use Rtypes.h here, as these data get used by root
        // global event info
        Int_t    fRun;
        Int_t    fSubRun;
        Int_t    fEvent;

        // MCTruth data; this analysis assumes that there is only 1 neutrino
        Int_t    fNeutrinoType;
        Int_t    fCCNC;
        Int_t    fMode;
        Int_t    fInteractionType;
        Float_t  fMCVertexX;
        Float_t  fMCVertexY;
        Float_t  fMCVertexZ;

        // MCParticle data
        std::vector<Int_t>              fMCP_PDG;
        std::vector<Float_t>            fMCP_X;
        std::vector<Float_t>            fMCP_Y;
        std::vector<Float_t>            fMCP_Z;
        std::vector<Float_t>            fMCP_PX;
        std::vector<Float_t>            fMCP_PY;
        std::vector<Float_t>            fMCP_PZ;
        std::vector<Float_t>            fMCP_Time;

        // track data
        std::vector<Int_t>              fBackMatchNo;
        std::vector<ULong64_t>          fTrackIDNumber;
        std::vector<Float_t>            fTrackX;
        std::vector<Float_t>            fTrackY;
        std::vector<Float_t>            fTrackZ;
        std::vector<Float_t>            fTrackPX;
        std::vector<Float_t>            fTrackPY;
        std::vector<Float_t>            fTrackPZ;
        std::vector<Int_t>              fTrackQ;
        std::vector<Float_t>            fTrackLen;
        std::vector<Int_t>              fNTPCClustersOnTrack;
        std::vector<Float_t>            fTrackAvgIon;

        // Some MC->Track matching info
        std::vector<Int_t>              fNhitPicked;
        std::vector<Int_t>              fNhitMatched;
        std::vector<Float_t>            fTrackMCcosT;
        std::vector<Float_t>            fTrackMCdelX;
    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::MomentumPerformance::MomentumPerformance(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    fGeantLabel         = p.get<std::string>("GEANTLabel","geant");
    fTPCClusterLabel = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel         = p.get<std::string>("TrackLabel","track");

    fIonizTruncate     = p.get<float>("IonizTruncate",    0.70);



    consumesMany<std::vector<simb::MCTruth> >();
    consumes<std::vector<simb::MCParticle> >(fGeantLabel);

    consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);
    consumes<std::vector<rec::TrackIoniz>>(fTrackLabel);
    consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MomentumPerformance::beginJob() {

    art::ServiceHandle<geo::GeometryGAr> euclid;
    ItsInTulsa[0] = euclid->GetOriginX();
    ItsInTulsa[1] = euclid->GetOriginY();
    ItsInTulsa[2] = euclid->GetOriginZ();

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

    cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
    fBack = const_cast<cheat::BackTrackerCore*>(const_bt);



    fTree->Branch("Run",                    &fRun,              "Run/I");
    fTree->Branch("SubRun",                 &fSubRun,  	        "SubRun/I");
    fTree->Branch("Event",                  &fEvent,			"Event/I");

    fTree->Branch("NType",                  &fNeutrinoType,		"Event/I");
    fTree->Branch("CCNC",                   &fCCNC, 	        "Event/I");
    fTree->Branch("Mode",                   &fMode, 		    "Event/I");
    fTree->Branch("InterT",                 &fInteractionType,	"Event/I");
    fTree->Branch("MCVertX",                &fMCVertexX,        "Event/F");
    fTree->Branch("MCVertY",                &fMCVertexY,        "Event/F");
    fTree->Branch("MCVertZ",                &fMCVertexZ,        "Event/F");

    fTree->Branch("MCP_PDG",                &fMCP_PDG);
    fTree->Branch("MCP_X",                  &fMCP_X);
    fTree->Branch("MCP_Y",                  &fMCP_Y);
    fTree->Branch("MCP_Z",                  &fMCP_Z);
    fTree->Branch("MCP_PX",                 &fMCP_PX);
    fTree->Branch("MCP_PY",                 &fMCP_PY);
    fTree->Branch("MCP_PZ",                 &fMCP_PZ);
    fTree->Branch("MCP_Time",               &fMCP_Time);

    fTree->Branch("BackMatchNo",            &fBackMatchNo);
    fTree->Branch("TrackIDNumber",          &fTrackIDNumber);
    fTree->Branch("TrackX",                 &fTrackX);
    fTree->Branch("TrackY",                 &fTrackY);
    fTree->Branch("TrackZ",                 &fTrackZ);
    fTree->Branch("TrackPX",                &fTrackPX);
    fTree->Branch("TrackPY",                &fTrackPY);
    fTree->Branch("TrackPZ",                &fTrackPZ);
    fTree->Branch("TrackQ",                 &fTrackQ);
    fTree->Branch("TrackLen",               &fTrackLen);
    fTree->Branch("NTPCClustersOnTrack",    &fNTPCClustersOnTrack);
    fTree->Branch("TrackAvgIon",            &fTrackAvgIon);

    fTree->Branch("NhitPicked",             &fNhitPicked);
    fTree->Branch("NhitMatched",            &fNhitMatched);
    fTree->Branch("TrackMCcosT",            &fTrackMCcosT);
    fTree->Branch("TrackMCdelX",            &fTrackMCdelX);

    return;
}  // End of :MomentumPerformance::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MomentumPerformance::analyze(art::Event const & event) {

    ClearVectors();
    FillVectors(event);
    fTree->Fill();
    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MomentumPerformance::ClearVectors() {
    // clear out all our vectors
    fMCP_PDG.clear();
    fMCP_X.clear();
    fMCP_Y.clear();
    fMCP_Z.clear();
    fMCP_PX.clear();
    fMCP_PY.clear();
    fMCP_PZ.clear();
    fMCP_Time.clear();

    fBackMatchNo.clear();
    fTrackIDNumber.clear();
    fTrackX.clear();
    fTrackY.clear();
    fTrackZ.clear();
    fTrackPX.clear();
    fTrackPY.clear();
    fTrackPZ.clear();
    fTrackQ.clear();
    fTrackLen.clear();
    fTrackAvgIon.clear();
    fNTPCClustersOnTrack.clear();

    fTrackMCcosT.clear();
    fTrackMCdelX.clear();

    return;
} // end :MomentumPerformance::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::MomentumPerformance::FillVectors(art::Event const& event) {



    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory.  Want to get all
    // MCTruths, regardless of generator label.

    auto mctruthHandles = event.getMany<std::vector<simb::MCTruth> >();

    if (mctruthHandles.size()!=1) {
        throw cet::exception("MomentumPerformance") << " Need just 1 simb::MCTruth"
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    auto MCPHandle = event.getHandle< std::vector<simb::MCParticle> >(fGeantLabel);
    if (!MCPHandle) {
        throw cet::exception("MomentumPerformance") << " No simb::MCParticle"
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // Get handles for Tracks and their ionizations; also Assn's to TPCClusters, TrackIoniz
    art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;

    auto TrackHandle = event.getHandle< std::vector<rec::Track> >(fTrackLabel);
    if (!TrackHandle) {
        throw cet::exception("MomentumPerformance") << " No rec::Track"
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    auto TrackIonHandle = event.getHandle< std::vector<rec::TrackIoniz> >(fTrackLabel);
    if (!TrackIonHandle) {
        throw cet::exception("MomentumPerformance") << " No rec::TrackIoniz"
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    findIonization = new art::FindOneP<rec::TrackIoniz>(TrackHandle,event,fTrackLabel);





    // =============  Pull art handles =========================================
    fRun    = event.run();
    fSubRun = event.subRun();
    fEvent  = event.id().event();



    // save MCTruth info.  This analysis is designed for single-interaction MC
    for ( simb::MCTruth const& mct : *(mctruthHandles.at(0)) ) {
        if (!mct.NeutrinoSet()) return;
        simb::MCNeutrino nuw = mct.GetNeutrino();
        fNeutrinoType     = nuw.Nu().PdgCode();
        fCCNC             = nuw.CCNC();
        fMode             = nuw.Mode();
        fInteractionType = nuw.InteractionType();
        fMCVertexX         = nuw.Nu().EndX();
        fMCVertexY         = nuw.Nu().EndY();
        fMCVertexZ         = nuw.Nu().EndZ();
    }  // end MC info from MCTruth

    typedef int TrkId;
    std::unordered_map<TrkId, Int_t> TrackIdToIndex;
    Int_t index = 0;
    for ( auto const& mcp : *MCPHandle ) {
        int TrackId = mcp.TrackId();
        TrackIdToIndex[TrackId] = index++;
    }



    for ( simb::MCParticle mcp : *MCPHandle ) {
        // If mcp.Mother() == 0, particle is from initial vertex;
        TrkId momTrkId = mcp.Mother();
        int momPDG = 0;
        if (momTrkId>0) {
            if(TrackIdToIndex.find(momTrkId) != TrackIdToIndex.end()){
                Int_t momIndex = TrackIdToIndex[momTrkId];
                momPDG   = (*MCPHandle).at(momIndex).PdgCode();
            } else {
                throw cet::exception("MomentumPerformance") << " mcp trackid "
                    << mcp.TrackId() << " pdg code " << mcp.PdgCode()
                    << " could not find mother trk id " << momTrkId
                    << " in the TrackIdToIndex map; creating process is "
                    << mcp.Process() << "\nLine " << __LINE__ << " in file "
                    << __FILE__ << std::endl;
            }
        }
        if (momPDG != 0) continue;



        // Pick up muons, electrons, charged pions, protons
        int thisPDG = mcp.PdgCode();
        bool isTrackable =
            abs(thisPDG)==13 || abs(thisPDG)==11 || abs(thisPDG)==211 || abs(thisPDG)==2212;
        if (!isTrackable) continue;

        const TLorentzVector& positionMCP = mcp.Position(0);
        const TLorentzVector& momentumMCP = mcp.Momentum(0);



        // OK get the matching tracks from the backtracker.  First collect yer tracks
        std::vector<art::Ptr<rec::Track>> allTracks;
        auto const trackPtrMaker = art::PtrMaker<rec::Track>(event,TrackHandle.id());
        for (size_t iTrack=0; iTrack<TrackHandle->size(); ++iTrack ) {
            allTracks.push_back( trackPtrMaker(iTrack) );
        }

        std::vector<art::Ptr<rec::Track>> matchedTracks;
        matchedTracks = fBack->MCParticleToTracks(&mcp,allTracks);



        // Which track end you want?
        float minDist = 1e6;
        int pickedTrack = -1;
        rec::TrackEnd kate = rec::TrackEndBeg;
        size_t iTrack;
        for (iTrack=0; iTrack<matchedTracks.size(); ++iTrack) {
            rec::Track track = *(matchedTracks[iTrack]);
            float distStart = std::hypot(track.Vertex()[0] -positionMCP[0],
                                         track.Vertex()[1] -positionMCP[1],
                                         track.Vertex()[2] -positionMCP[2]);
            float distEnd   = std::hypot(track.End()[0]       -positionMCP[0],
                                         track.End()[1]       -positionMCP[1],
                                         track.End()[2]       -positionMCP[2]);
            if (distStart<distEnd) {
                if (distStart<minDist) {
                    pickedTrack = iTrack;
                    minDist        = distStart;
                    kate        = rec::TrackEndBeg;
                }
            } else {
                if (distEnd<minDist) {
                    pickedTrack = iTrack;
                    minDist        = distEnd;
                    kate        = rec::TrackEndEnd;
                }
            }
        }
        if (pickedTrack == -1) continue;



        // Save that tracks info; but first you need the corresponding MC info
        fMCP_PDG.push_back(thisPDG);
        fMCP_X.push_back(positionMCP.X());
        fMCP_Y.push_back(positionMCP.Y());
        fMCP_Z.push_back(positionMCP.Z());
        fMCP_PX.push_back(momentumMCP.Px());
        fMCP_PY.push_back(momentumMCP.Py());
        fMCP_PZ.push_back(momentumMCP.Pz());
        fMCP_Time.push_back(mcp.T());

        fBackMatchNo.push_back(iTrack);
        rec::Track theTrack = *(matchedTracks[pickedTrack]);
        fTrackIDNumber.push_back(theTrack.getIDNumber());

        TVector3* trackInPhase;    TVector3* trackInSpace;

        if (kate==rec::TrackEndBeg) {
            fTrackX.push_back  (theTrack.Vertex()[0]);
            fTrackY.push_back  (theTrack.Vertex()[1]);
            fTrackZ.push_back  (theTrack.Vertex()[2]);
            fTrackPX.push_back (theTrack.Momentum_beg()*theTrack.VtxDir()[0]);
            fTrackPY.push_back (theTrack.Momentum_beg()*theTrack.VtxDir()[1]);
            fTrackPZ.push_back (theTrack.Momentum_beg()*theTrack.VtxDir()[2]);
            fTrackQ.push_back  (theTrack.ChargeBeg());
            fTrackLen.push_back(theTrack.LengthBackward());
            trackInPhase = new TVector3( theTrack.VtxDir() );
            trackInSpace = new TVector3( theTrack.Vertex() );
        } else {
            fTrackX.push_back  (theTrack.End()[0]);
            fTrackY.push_back  (theTrack.End()[1]);
            fTrackZ.push_back  (theTrack.End()[2]);
            fTrackPX.push_back (theTrack.Momentum_end()*theTrack.EndDir()[0]);
            fTrackPY.push_back (theTrack.Momentum_end()*theTrack.EndDir()[1]);
            fTrackPZ.push_back (theTrack.Momentum_end()*theTrack.EndDir()[2]);
            fTrackQ.push_back  (theTrack.ChargeEnd());
            fTrackLen.push_back(theTrack.LengthForward());
            trackInPhase = new TVector3( theTrack.EndDir() );
            trackInSpace = new TVector3( theTrack.End() );
        }
        fNTPCClustersOnTrack.push_back(theTrack.NHits());


        // What fraction of the hits in all tracks matched to this MCParticle
        // are in the pickedTrack?
        int num, den;        num = den = 0;
        for (int iTrack_local=0; iTrack_local<(int)matchedTracks.size(); ++iTrack_local) {
            rec::Track thisTrack = *(matchedTracks[iTrack_local]);
            std::vector<art::Ptr<rec::Hit>> umlazi = fBack->TrackToHits(&thisTrack);
            den += umlazi.size();
            if (iTrack_local == pickedTrack) num += umlazi.size();
        }
        fNhitPicked.push_back(num);
        fNhitMatched.push_back(den);

        // And spatial track-matching info
        float cosT = (*trackInPhase)[0]*momentumMCP.Px()
                    +(*trackInPhase)[1]*momentumMCP.Py()
                    +(*trackInPhase)[2]*momentumMCP.Pz();
        cosT /= momentumMCP.P();
        fTrackMCcosT.push_back(cosT);

        TVector3 vecX( (*trackInSpace)[0] -positionMCP.X(),
                       (*trackInSpace)[1] -positionMCP.Y(),
                       (*trackInSpace)[2] -positionMCP.Z() );
        float delX = vecX.Cross(*trackInPhase).Mag();
        fTrackMCdelX.push_back(delX);



        // Need to determine which track in TrackHandle is the one
        // picked from matchedTracks
        int iPickedTrack = -1;
        for (size_t iTrack_local=0; iTrack_local<TrackHandle->size(); ++iTrack_local ) {
            if ( allTracks[iTrack_local]->getIDNumber() == theTrack.getIDNumber() ) {
                iPickedTrack = iTrack_local;
                break;
            }
        }
        if (iPickedTrack>=0 && findIonization->isValid()) {
            // No calibration for now.  Someday this should all be in reco
            rec::TrackIoniz ionization = *(findIonization->at(iPickedTrack));
            float avgIonF, avgIonB;
            processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
            if (kate==rec::TrackEndBeg) {
                fTrackAvgIon.push_back( avgIonF );
            } else {
                fTrackAvgIon.push_back( avgIonB );
            }
        } else {
            // must push_back something so that fTrackAvgIon is of correct size.
            fTrackAvgIon.push_back( 0.0 );
        }

    } // end loop on MCParticles

} // end MomentumPerformance::FillVectors



//==============================================================================
//==============================================================================
//==============================================================================
// Process ionization.  Copy of anatree_module code as of Apr 2020.
void gar::MomentumPerformance::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
float& forwardIonVal, float& backwardIonVal) {

    // NO CALIBRATION SERVICE FOR NOW

    std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
    forwardIonVal = processOneDirection(SigData, ionizeTruncate);

    SigData = ion.getBAK_dSigdXs();
    backwardIonVal = processOneDirection(SigData, ionizeTruncate);

    return;
}



float gar::MomentumPerformance::processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

    std::vector<std::pair<float,float>> dEvsX;    // Ionization vs distance along track

    // The first hit on the track never had its ionization info stored.  Not a problem
    // really.  Each pair is a hit and the step along the track that ends at the hit
    // For the last hit, just take the step from the n-1 hit; don't guess some distance
    // to (nonexistant!) n+1 hit.  Using pointer arithmetic because you are a real K&R
    // C nerd!  Except that C++ doesn't know you are such a nerd and if
    //  SigData.size()==0, then SigData.end()-1 is 0xFFFFFFFFFFFFFFF8.
    if (SigData.size()==0) return 0.0;
    float distAlongTrack = 0;
    std::vector<std::pair<float,float>>::iterator littlebit = SigData.begin();
    for (; littlebit<(SigData.end()-1); ++littlebit) {
        float dE =   std::get<0>(*littlebit);
        // tpctrackfit2_module.cc fills the TrackIoniz data product so that
        // this quantity is really dL > 0 not dX, a coordinate on the drift axis
        float dX  = std::get<1>(*littlebit);
        distAlongTrack += dX;    // But count full step to get hit position on track
        // Take dX to be 1/2 the previous + last segment
        dX += std::get<1>(*(littlebit+1));
        float dEdX = dE/(0.5*dX);

        std::pair pushme = std::make_pair(dEdX,distAlongTrack);
        dEvsX.push_back( pushme );
    }

    // Get the truncated mean; first sort then take mean
    std::sort(dEvsX.begin(),dEvsX.end(), lessThan_byE);

    // Get the dEdX vs length data, truncated.
    int goUpTo = ionizeTruncate * dEvsX.size() +0.5;
    if (goUpTo > (int)dEvsX.size()) goUpTo = dEvsX.size();
    int i = 1;        float returnvalue = 0;
    littlebit = dEvsX.begin();
    for (; littlebit<dEvsX.end(); ++littlebit) {
        returnvalue += std::get<0>(*littlebit);
        ++i;
        if (i>goUpTo) break;
    }
    returnvalue /= goUpTo;
    return returnvalue;
}





DEFINE_ART_MODULE(gar::MomentumPerformance)
