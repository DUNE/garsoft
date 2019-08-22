////////////////////////////////////////////////////////////////////////////////
// Class:       anatree
// Plugin Type: analyzer (art v2_11_02)
// File:        anatree_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Additions from Leo Bellantoni, 2019
////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "RawDataProducts/CaloRawDigit.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <string>
#include <vector>



namespace gar {

    class anatree : public art::EDAnalyzer {
    public:
        explicit anatree(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        anatree(anatree const &) = delete;
        anatree(anatree &&) = delete;
        anatree & operator = (anatree const &) = delete;
        anatree & operator = (anatree &&) = delete;

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

        // Compute T for coherent pion analysis
        float computeT( simb::MCTruth theMCTruth );



        // Input data labels
        std::vector<std::string> fGeneratorLabels;
        std::vector<std::string> fGENIEGeneratorLabels;
        std::string fGeantLabel;
        std::string fHitLabel;
        std::string fTPCClusterLabel;
        std::string fTrackLabel;
        std::string fVertexLabel;
        std::string fRawCaloHitLabel;
        std::string fCaloHitLabel;
        std::string fClusterLabel;
        std::string fECALAssnLabel;

       // Optionally keep/drop parts of the analysis tree
        bool fWriteMCinfo;        ///< Info from MCTruth, GTruth     Default=true
        bool fWriteMCPTrajectory; ///< MCP Trajectory                Default=true
        bool fWriteMCCaloInfo;    ///< Write MC info for calorimeter Default=true

        bool fWriteHits;          ///< Write info about TPC Hits     Default=false
        bool fWriteTPCClusters;   ///< Write TPCClusters info        Default=true
        bool fWriteTracks;        ///< Start/end X, P for tracks     Default=true
        bool fWriteVertices;      ///< Reco vertexes & their tracks  Default=true

        bool fWriteCaloDigits;    ///< Raw digits for calorimetry.   Default=false
        bool fWriteCaloHits;      ///< Write ECAL hits.              Default=true
        bool fWriteCaloClusters;  ///< Write ECAL clusters.          Default=true
        bool fWriteMatchedTracks; ///< Write ECAL-track Assns        Default=true

        // Truncation parameter for dE/dx (average this fraction lowest readings)
        float fIonizTruncate;     ///<                               Default=1.00;
        bool fWriteCohInfo;       ///< MC level t for coherent pi+.  Default=false

        // the analysis tree
        TTree *fTree;



        // global event info
        Int_t fEvent;        ///< number of the event being processed
        Int_t fRun;          ///< number of the run being processed
        Int_t fSubRun;       ///< number of the sub-run being processed

        // MCTruth data.
        // GENIE kinematics computed ignoring Fermi momentum and the off-shellness
        // of the bound nucleon.  Well, that's what the documentation says.  Might
        // not be true!  But to get the right kinematics here, t has to be computed.
        // Mode and InteractionType given in simb::MCNeutrino.h of the nusimdata
        // product.
        // Use Rtypes.h here, as these data get used by root

        std::vector<Int_t>       fNeutrinoType;
        std::vector<Int_t>       fCCNC;
        std::vector<Int_t>       fMode;
        std::vector<Int_t>       fInteractionType;
        std::vector<Float_t>     fQ2;
        std::vector<Float_t>     fW;
        std::vector<Float_t>     fX;
        std::vector<Float_t>     fY;
        std::vector<Float_t>     fTheta;
        std::vector<Float_t>     fT;
        std::vector<Float_t>     fMCVertexX;
        std::vector<Float_t>     fMCVertexY;
        std::vector<Float_t>     fMCVertexZ;
        std::vector<Float_t>     fMCnuPx;
        std::vector<Float_t>     fMCnuPy;
        std::vector<Float_t>     fMCnuPz;

        // GTruth data
        std::vector<Int_t>       fGint;
        std::vector<Int_t>       fTgtPDG;
        std::vector<Float_t>     fWeight;
        std::vector<Float_t>     fgT;

        // MCParticle data
        std::vector<Int_t>       fMCPTrkID;
        std::vector<Int_t>       fMCPDG;
        std::vector<Int_t>       fMCMother;
        std::vector<Int_t>       fMCPDGMother;
        std::vector<Float_t>     fMCPStartX;
        std::vector<Float_t>     fMCPStartY;
        std::vector<Float_t>     fMCPStartZ;
        std::vector<Float_t>     fMCPTime;
        std::vector<Float_t>     fMCPStartPX;
        std::vector<Float_t>     fMCPStartPY;
        std::vector<Float_t>     fMCPStartPZ;
        std::vector<Float_t>     fMCPEndX;
        std::vector<Float_t>     fMCPEndY;
        std::vector<Float_t>     fMCPEndZ;
        std::vector<Float_t>     fMCPEndPX;
        std::vector<Float_t>     fMCPEndPY;
        std::vector<Float_t>     fMCPEndPZ;
        std::vector<std::string> fMCPProc;
        std::vector<std::string> fMCPEndProc;

        // track trajectory of MCP
        std::vector<Float_t>     fTrajMCPX;
        std::vector<Float_t>     fTrajMCPY;
        std::vector<Float_t>     fTrajMCPZ;
        std::vector<Float_t>     fTrajMCPT;
        std::vector<Float_t>     fTrajMCPE;
        std::vector<Int_t>       fTrajMCPTrajIndex;

        // sim calo hit data
        UInt_t                   fSimnHits;
        std::vector<Float_t>     fSimHitX;
        std::vector<Float_t>     fSimHitY;
        std::vector<Float_t>     fSimHitZ;
        std::vector<Float_t>     fSimHitTime;
        std::vector<Float_t>     fSimHitEnergy;
        std::vector<Int_t>       fSimHitTrackID;
        std::vector<ULong64_t>   fSimHitCellID;
        Float_t                  fSimEnergySum;

        // Hit data
        std::vector<Float_t>     fHitX;
        std::vector<Float_t>     fHitY;
        std::vector<Float_t>     fHitZ;
        std::vector<Float_t>     fHitSig;
        std::vector<Float_t>     fHitRMS;

        // TPCCluster data
        std::vector<Float_t>     fTPCClusterX;
        std::vector<Float_t>     fTPCClusterY;
        std::vector<Float_t>     fTPCClusterZ;
        std::vector<Float_t>     fTPCClusterSig;
        std::vector<Float_t>     fTPCClusterRMS;
        std::vector<Int_t>       fTPCClusterTrkIDNumber;

        // track data
        std::vector<ULong64_t>   fTrackIDNumber;
        std::vector<Float_t>     fTrackStartX;
        std::vector<Float_t>     fTrackStartY;
        std::vector<Float_t>     fTrackStartZ;
        std::vector<Float_t>     fTrackStartPX;
        std::vector<Float_t>     fTrackStartPY;
        std::vector<Float_t>     fTrackStartPZ;
        std::vector<Int_t>       fTrackStartQ;

        std::vector<Float_t>     fTrackEndX;
        std::vector<Float_t>     fTrackEndY;
        std::vector<Float_t>     fTrackEndZ;
        std::vector<Float_t>     fTrackEndPX;
        std::vector<Float_t>     fTrackEndPY;
        std::vector<Float_t>     fTrackEndPZ;
        std::vector<Int_t>       fTrackEndQ;

        std::vector<Float_t>     fTrackLenF;
        std::vector<Float_t>     fTrackLenB;
        std::vector<Int_t>       fNTPCClustersOnTrack;
        std::vector<Float_t>     fTrackAvgIonF;
        std::vector<Float_t>     fTrackAvgIonB;

        // vertex branches
        std::vector<ULong64_t>   fVertexIDNumber;
        std::vector<Float_t>     fVertexX;
        std::vector<Float_t>     fVertexY;
        std::vector<Float_t>     fVertexZ;
        std::vector<ULong64_t>   fVertexT;
        std::vector<Int_t>       fVertexN;
        std::vector<Int_t>       fVertexQ;

        std::vector<ULong64_t>   fVTAssn_VertIDNumber;     // Being the vertex which this Assn belongs to
        std::vector<ULong64_t>   fVTAssn_TrackIDNumber;
        std::vector<int>         fVTAssn_TrackEnd;

        // raw calo digits data
        UInt_t                   fDiginHits;
        std::vector<Float_t>     fDigiHitX;
        std::vector<Float_t>     fDigiHitY;
        std::vector<Float_t>     fDigiHitZ;
        std::vector<Float_t>     fDigiHitTime;
        std::vector<UInt_t>      fDigiHitADC;        // UInt_t is unsigned 32 bit integer
        std::vector<ULong64_t>   fDigiHitCellID;     // ULong64_t is size_t on 64 bit machines

        // reco calo hit data
        std::vector<ULong64_t>   fReconHitIDNumber;
        std::vector<Float_t>     fRecoHitX;
        std::vector<Float_t>     fRecoHitY;
        std::vector<Float_t>     fRecoHitZ;
        std::vector<Float_t>     fRecoHitTime;
        std::vector<Float_t>     fRecoHitEnergy;
        std::vector<ULong64_t>   fRecoHitCellID;
        Float_t                  fRecoEnergySum;

        // calo cluster data
        std::vector<ULong64_t>   fClusterIDNumber;
        std::vector<UInt_t>      fClusterNhits;
        std::vector<Float_t>     fClusterEnergy;
        std::vector<Float_t>     fClusterX;
        std::vector<Float_t>     fClusterY;
        std::vector<Float_t>     fClusterZ;
        std::vector<Float_t>     fClusterTheta;
        std::vector<Float_t>     fClusterPhi;
        std::vector<Float_t>     fClusterPID;
        // std::vector<Float_t> fClusterShape;
        std::vector<Float_t>     fClusterMainAxisX;
        std::vector<Float_t>     fClusterMainAxisY;
        std::vector<Float_t>     fClusterMainAxisZ;

        // ECAL cluster to track association info
        std::vector<ULong64_t>   fCALAssn_ClusIDNumber;   // Being the cluster which this Assn belongs to
        std::vector<ULong64_t>   fCALAssn_TrackIDNumber;  // The rec::TrackEnd (see Track.h) that extrapolated to cluster
        std::vector<int>         fCALAssn_TrackEnd;
    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatree::anatree(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    bool usegenlabels =
        p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
    if (!usegenlabels) fGeneratorLabels.clear();

    bool usegeniegenlabels  =
        p.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
    if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

    fGeantLabel       = p.get<std::string>("GEANTLabel","geant");
    fHitLabel         = p.get<std::string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<std::string>("TrackLabel","track");
    fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
    fRawCaloHitLabel  = p.get<std::string>("RawCaloHitLabel","daqecal");
    fCaloHitLabel     = p.get<std::string>("CaloHitLabel","calohit");
    fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");
    fECALAssnLabel    = p.get<std::string>("ECALAssnLabel","trkecalassn");



    // What to write
    fWriteMCinfo              = p.get<bool>("WriteMCinfo",       true);
    fWriteMCPTrajectory       = p.get<bool>("WriteMCPTrajectory",true);
    fWriteMCCaloInfo          = p.get<bool>("WriteMCCaloInfo",   true);

    fWriteHits                = p.get<bool>("WriteHits",         false);
    fWriteTPCClusters         = p.get<bool>("WriteTPCClusters",  true);
    fWriteTracks              = p.get<bool>("WriteTracks",       true);
    fWriteVertices            = p.get<bool>("WriteVertices",     true);

    fWriteCaloDigits          = p.get<bool>("WriteCaloDigits",   false);
    fWriteCaloHits            = p.get<bool>("WriteCaloHits",     true);
    fWriteCaloClusters        = p.get<bool>("WriteCaloClusters", true);
    fWriteMatchedTracks       = p.get<bool>("WriteMatchedTracks",true);

    fIonizTruncate            = p.get<float>("IonizTruncate",    0.70);
    fWriteCohInfo             = p.get<bool> ("WriteCohInfo",     false);

    consumes<std::vector<simb::MCParticle> >(fGeantLabel);

    if (usegenlabels) {
        for (size_t i=0; i<fGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::MCTruth> >();
    }

    if (usegeniegenlabels) {
        for (size_t i=0; i<fGENIEGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::GTruth> >(fGENIEGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::GTruth> >();
    }

    //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);

    consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Hit> >(fHitLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);
    consumes<std::vector<rec::Vertex> >(fVertexLabel);
    consumes<art::Assns<rec::Track, rec::Vertex> >(fVertexLabel);
    consumes<std::vector<rec::TrackIoniz>>(fTrackLabel);
    consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

    consumes<std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
    consumes<std::vector<raw::CaloRawDigit> >(fRawCaloHitLabel);
    consumes<std::vector<rec::CaloHit> >(fCaloHitLabel);
    consumes<std::vector<rec::Cluster> >(fClusterLabel);
    consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

    fTree->Branch("Run",           &fRun,         "Run/I");
    fTree->Branch("SubRun",        &fSubRun,      "SubRun/I");
    fTree->Branch("Event",         &fEvent,       "Event/I");

    if (fWriteMCinfo) {
        fTree->Branch("NType",       &fNeutrinoType);
        fTree->Branch("CCNC",        &fCCNC);
        fTree->Branch("Mode",        &fMode);
        fTree->Branch("InterT",      &fInteractionType);
        fTree->Branch("MC_Q2",       &fQ2);
        fTree->Branch("MC_W",        &fW);
        fTree->Branch("MC_X",        &fX);
        fTree->Branch("MC_Y",        &fY);
        fTree->Branch("MC_Theta",    &fTheta);
        if (fWriteCohInfo) fTree->Branch("MC_T", &fT);
        fTree->Branch("MCVertX",     &fMCVertexX);
        fTree->Branch("MCVertY",     &fMCVertexY);
        fTree->Branch("MCVertZ",     &fMCVertexZ);
        fTree->Branch("MCNuPx",      &fMCnuPx);
        fTree->Branch("MCNuPy",      &fMCnuPy);
        fTree->Branch("MCNuPz",      &fMCnuPz);

        fTree->Branch("Gint",        &fGint);
        fTree->Branch("TgtPDG",      &fTgtPDG);
        fTree->Branch("Weight",      &fWeight);
        fTree->Branch("GT_T",        &fgT);

        fTree->Branch("MCPTrkID",    &fMCPTrkID);
        fTree->Branch("PDG",         &fMCPDG);
        fTree->Branch("Mother",      &fMCMother);
        fTree->Branch("PDGMother",   &fMCPDGMother);
        fTree->Branch("MCPStartX",   &fMCPStartX);
        fTree->Branch("MCPStartY",   &fMCPStartY);
        fTree->Branch("MCPStartZ",   &fMCPStartZ);
        fTree->Branch("MCPTime",     &fMCPTime);
        fTree->Branch("MCPStartPX",  &fMCPStartPX);
        fTree->Branch("MCPStartPY",  &fMCPStartPY);
        fTree->Branch("MCPStartPZ",  &fMCPStartPZ);
        fTree->Branch("MCPEndX",     &fMCPEndX);
        fTree->Branch("MCPEndY",     &fMCPEndY);
        fTree->Branch("MCPEndZ",     &fMCPEndZ);
        fTree->Branch("MCPEndPX",    &fMCPEndPX);
        fTree->Branch("MCPEndPY",    &fMCPEndPY);
        fTree->Branch("MCPEndPZ",    &fMCPEndPZ);
        fTree->Branch("MCPProc",     &fMCPProc);
        fTree->Branch("MCPEndProc",  &fMCPEndProc);
    }

    if (fWriteMCPTrajectory) {
        // Checking for parameter file validity only once to simplify
        // later code in anatree::process etc.
        if (!fWriteMCinfo) {
            throw cet::exception("anatree")
            << " fWriteMCPTrajectory, but !fWriteMCinfo."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        } else {
            //Write MCP Trajectory
            fTree->Branch("TrajMCPX",          &fTrajMCPX);
            fTree->Branch("TrajMCPY",          &fTrajMCPY);
            fTree->Branch("TrajMCPZ",          &fTrajMCPZ);
            fTree->Branch("TrajMCPT",          &fTrajMCPT);
            fTree->Branch("TrajMCPE",          &fTrajMCPE);
            fTree->Branch("TrajMCPTrajIndex",  &fTrajMCPTrajIndex);
        }
    }

    if (fWriteMCCaloInfo) {
        if (!fWriteMCinfo) {
            throw cet::exception("anatree")
            << " fWriteMCCaloInfo, but !fWriteMCinfo."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        } else {
            // Write calorimetry MC information
            fTree->Branch("SimnHits",     &fSimnHits);
            fTree->Branch("SimHitX",      &fSimHitX);
            fTree->Branch("SimHitY",      &fSimHitY);
            fTree->Branch("SimHitZ",      &fSimHitZ);
            fTree->Branch("SimHitTime",   &fSimHitTime);
            fTree->Branch("SimHitEnergy", &fSimHitEnergy);
            fTree->Branch("SimHitTrkID",  &fSimHitTrackID);
            fTree->Branch("SimHitCellID", &fSimHitCellID);
            fTree->Branch("SimEnergySum", &fSimEnergySum);
        }
    }

    if (fWriteHits) {
        fTree->Branch("HitX",        &fHitX);
        fTree->Branch("HitY",        &fHitY);
        fTree->Branch("HitZ",        &fHitZ);
        fTree->Branch("HitSig",      &fHitSig);
        fTree->Branch("HitRMS",      &fHitRMS);
    }

    if (fWriteTPCClusters) {
        fTree->Branch("TPCClusterX",           &fTPCClusterX);
        fTree->Branch("TPCClusterY",           &fTPCClusterY);
        fTree->Branch("TPCClusterZ",           &fTPCClusterZ);
        fTree->Branch("TPCClusterSig",         &fTPCClusterSig);
        fTree->Branch("TPCClusterRMS",         &fTPCClusterRMS);
        fTree->Branch("TPCClusterTrkIDNumber", &fTPCClusterTrkIDNumber);
    }

    if (fWriteTracks) {                         // All position, momentum, etc
        fTree->Branch("TrackIDNumber",      &fTrackIDNumber);

        fTree->Branch("TrackStartX",        &fTrackStartX);
        fTree->Branch("TrackStartY",        &fTrackStartY);
        fTree->Branch("TrackStartZ",        &fTrackStartZ);
        fTree->Branch("TrackStartPX",       &fTrackStartPX);
        fTree->Branch("TrackStartPY",       &fTrackStartPY);
        fTree->Branch("TrackStartPZ",       &fTrackStartPZ);
        fTree->Branch("TrackStartQ",        &fTrackStartQ);

        fTree->Branch("TrackEndX",          &fTrackEndX);
        fTree->Branch("TrackEndY",          &fTrackEndY);
        fTree->Branch("TrackEndZ",          &fTrackEndZ);
        fTree->Branch("TrackEndPX",         &fTrackEndPX);
        fTree->Branch("TrackEndPY",         &fTrackEndPY);
        fTree->Branch("TrackEndPZ",         &fTrackEndPZ);
        fTree->Branch("TrackEndQ",          &fTrackEndQ);

        fTree->Branch("TrackLenF",          &fTrackLenF);
        fTree->Branch("TrackLenB",          &fTrackLenB);
        fTree->Branch("NTPCClustersOnTrack",&fNTPCClustersOnTrack);
        fTree->Branch("TrackAvgIonF",       &fTrackAvgIonF);
        fTree->Branch("TrackAvgIonB",       &fTrackAvgIonB);

    }

    if (fWriteVertices) {                     // Reco'd verts & their track-ends
        if (!fWriteTracks) {
            throw cet::exception("anatree")
            << " fWriteVertices, but !fWriteTracks."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        } else {
            fTree->Branch("VertIDNumber",    &fVertexIDNumber);
            fTree->Branch("VertX",           &fVertexX);
            fTree->Branch("VertY",           &fVertexY);
            fTree->Branch("VertZ",           &fVertexZ);
            fTree->Branch("VertT",           &fVertexT);
            fTree->Branch("VertN",           &fVertexN);
            fTree->Branch("VertQ",           &fVertexQ);

            fTree->Branch("VT_VertIDNumber", &fVTAssn_VertIDNumber);
            fTree->Branch("VT_TrackIDNumber",&fVTAssn_TrackIDNumber);
            fTree->Branch("VT_TrackEnd",     &fVTAssn_TrackEnd);
        }
    }

    if (fWriteCaloDigits) {                       // Write calorimetry digits
        fTree->Branch("DiginHits",        &fDiginHits);
        fTree->Branch("DigiHitX",         &fDigiHitX);
        fTree->Branch("DigiHitY",         &fDigiHitY);
        fTree->Branch("DigiHitZ",         &fDigiHitZ);
        fTree->Branch("DigiHitTime",      &fDigiHitTime);
        fTree->Branch("DigiHitADC",       &fDigiHitADC);
        fTree->Branch("DigiHitCellID",    &fDigiHitCellID);
    }

    if (fWriteCaloHits) {                         // Write calorimetry hits
        fTree->Branch("ReconHitIDNumber", &fReconHitIDNumber);
        fTree->Branch("RecoHitX",         &fRecoHitX);
        fTree->Branch("RecoHitY",         &fRecoHitY);
        fTree->Branch("RecoHitZ",         &fRecoHitZ);
        fTree->Branch("RecoHitTime",      &fRecoHitTime);
        fTree->Branch("RecoHitEnergy",    &fRecoHitEnergy);
        fTree->Branch("RecoHitCellID",    &fRecoHitCellID);
        fTree->Branch("RecoEnergySum",    &fRecoEnergySum);
    }

    if (fWriteCaloClusters) {                     // Write calorimetry clusters
        fTree->Branch("ClusterIDHits",    &fClusterIDNumber);
        fTree->Branch("ClusterNhits",     &fClusterNhits);
        fTree->Branch("ClusterEnergy",    &fClusterEnergy);
        fTree->Branch("ClusterX",         &fClusterX);
        fTree->Branch("ClusterY",         &fClusterY);
        fTree->Branch("ClusterZ",         &fClusterZ);
        fTree->Branch("ClusterTheta",     &fClusterTheta);
        fTree->Branch("ClusterPhi",       &fClusterPhi);
        fTree->Branch("ClusterPID",       &fClusterPID);
        // fTree->Branch("ClusterShape",  &fClusterShape);
        fTree->Branch("ClusterMainAxisX", &fClusterMainAxisX);
        fTree->Branch("ClusterMainAxisY", &fClusterMainAxisY);
        fTree->Branch("ClusterMainAxisZ", &fClusterMainAxisZ);
    }

    if (fWriteMatchedTracks) {
        if (!fWriteTracks || !fWriteCaloClusters) {
            throw cet::exception("anatree")
            << " fWriteMatchedTracks, but (!fWriteTracks || !fWriteCaloClusters)."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        } else {
            fTree->Branch("ECALAssn_ClusIDNumber",   &fCALAssn_ClusIDNumber);
            fTree->Branch("ECALAssn_TrackIDNumber",  &fCALAssn_TrackIDNumber);
            fTree->Branch("ECALAssn_TrackEnd",       &fCALAssn_TrackEnd);
        }
    }
    return;
}  // End of :anatree::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::analyze(art::Event const & e) {

    ClearVectors();
    FillVectors(e);
    fTree->Fill();
    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::ClearVectors() {
    // clear out all our vectors
    if (fWriteMCinfo) {
        fNeutrinoType.clear();
        fCCNC.clear();
        fMode.clear();
        fInteractionType.clear();
        fQ2.clear();
        fW.clear();
        fX.clear();
        fY.clear();
        fTheta.clear();
        if (fWriteCohInfo) fT.clear();
        fMCVertexX.clear();
        fMCVertexY.clear();
        fMCVertexZ.clear();
        fMCnuPx.clear();
        fMCnuPy.clear();
        fMCnuPz.clear();

        fGint.clear();
        fTgtPDG.clear();
        fWeight.clear();
        fgT.clear();

        fMCPTrkID.clear();
        fMCPDG.clear();
        fMCMother.clear();
        fMCPDGMother.clear();
        fMCPStartX.clear();
        fMCPStartY.clear();
        fMCPStartZ.clear();
        fMCPTime.clear();
        fMCPStartPX.clear();
        fMCPStartPY.clear();
        fMCPStartPZ.clear();
        fMCPEndX.clear();
        fMCPEndY.clear();
        fMCPEndZ.clear();
        fMCPEndPX.clear();
        fMCPEndPY.clear();
        fMCPEndPZ.clear();
        fMCPProc.clear();
        fMCPEndProc.clear();
    }

    if (fWriteMCPTrajectory) {
        fTrajMCPX.clear();
        fTrajMCPY.clear();
        fTrajMCPZ.clear();
        fTrajMCPT.clear();
        fTrajMCPE.clear();
        fTrajMCPTrajIndex.clear();
    }

    if (fWriteMCCaloInfo) {
        fSimnHits = 0;
        fSimHitX.clear();
        fSimHitY.clear();
        fSimHitZ.clear();
        fSimHitTime.clear();
        fSimHitEnergy.clear();
        fSimHitTrackID.clear();
        fSimHitCellID.clear();
        fSimEnergySum = 0.;
    }

    if (fWriteHits) {
        fHitX.clear();
        fHitY.clear();
        fHitZ.clear();
        fHitSig.clear();
        fHitRMS.clear();
    }

    if (fWriteTPCClusters) {
        fTPCClusterX.clear();
        fTPCClusterY.clear();
        fTPCClusterZ.clear();
        fTPCClusterSig.clear();
        fTPCClusterRMS.clear();
        fTPCClusterTrkIDNumber.clear();
    }

    if (fWriteTracks) {
        fTrackIDNumber.clear();
        fTrackStartX.clear();
        fTrackStartY.clear();
        fTrackStartZ.clear();
        fTrackStartPX.clear();
        fTrackStartPY.clear();
        fTrackStartPZ.clear();
        fTrackStartQ.clear();
        fTrackEndX.clear();
        fTrackEndY.clear();
        fTrackEndZ.clear();
        fTrackEndPX.clear();
        fTrackEndPY.clear();
        fTrackEndPZ.clear();
        fTrackEndQ.clear();
        fTrackLenF.clear();
        fTrackLenB.clear();
        fTrackAvgIonF.clear();
        fTrackAvgIonB.clear();
        fNTPCClustersOnTrack.clear();
    }

    if (fWriteVertices) {
        fVertexIDNumber.clear();
        fVertexX.clear();
        fVertexY.clear();
        fVertexZ.clear();
        fVertexT.clear();
        fVertexN.clear();
        fVertexQ.clear();

        fVTAssn_VertIDNumber.clear();
        fVTAssn_TrackIDNumber.clear();
        fVTAssn_TrackEnd.clear();
    }

    if (fWriteCaloDigits) {
        fDiginHits = 0;
        fDigiHitX.clear();
        fDigiHitY.clear();
        fDigiHitZ.clear();
        fDigiHitTime.clear();
        fDigiHitADC.clear();
        fDigiHitCellID.clear();
    }

    if (fWriteCaloHits) {
        fReconHitIDNumber.clear();
        fRecoHitX.clear();
        fRecoHitY.clear();
        fRecoHitZ.clear();
        fRecoHitTime.clear();
        fRecoHitEnergy.clear();
        fRecoHitCellID.clear();
        fRecoEnergySum = 0.;
    }

    if (fWriteCaloClusters) {
        fClusterIDNumber.clear();
        fClusterNhits.clear();
        fClusterEnergy.clear();
        fClusterX.clear();
        fClusterY.clear();
        fClusterZ.clear();
        fClusterTheta.clear();
        fClusterPhi.clear();
        fClusterPID.clear();
        // fClusterShape.clear();
        fClusterMainAxisX.clear();
        fClusterMainAxisY.clear();
        fClusterMainAxisZ.clear();
    }

    if (fWriteMatchedTracks) {
        fCALAssn_ClusIDNumber.clear();
        fCALAssn_TrackIDNumber.clear();
        fCALAssn_TrackIDNumber.clear();
    }
} // end :anatree::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillVectors(art::Event const & e) {



    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory
    std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
    std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;
    art::Handle< std::vector<simb::MCParticle> > MCPHandle;
    art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;

    if (fWriteMCinfo) {
        if (fGeneratorLabels.size()<1) {
            e.getManyByType(mcthandlelist);    // get them all (even if there are none)
        } else {
            mcthandlelist.resize(fGeneratorLabels.size());
            for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
                // complain if we wanted a specific one but didn't find it
                if (!e.getByLabel(fGeneratorLabels.at(i),mcthandlelist.at(i)))
                    throw cet::exception("anatree") << " No simb::MCTruth branch." 
                         << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }

        if (fGENIEGeneratorLabels.size()<1) {
            e.getManyByType(gthandlelist);  // get them all (even if there are none)
        } else {
            gthandlelist.resize(fGENIEGeneratorLabels.size());
            for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i) {
                // complain if we wanted a specific one but didn't find it
                if (!e.getByLabel(fGENIEGeneratorLabels.at(i),gthandlelist.at(i)))
                    throw cet::exception("anatree") << " No simb::GTruth branch."
                        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }

        if (!e.getByLabel(fGeantLabel, MCPHandle)) {
            throw cet::exception("anatree") << " No simb::MCParticle branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handles for MCCaloInfo
    if (fWriteMCCaloInfo) {
        if (!e.getByLabel(fGeantLabel, SimHitHandle)) {
            throw cet::exception("anatree") << " No gar::sdp::CaloDeposit branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handle for TPC hit data
    art::Handle< std::vector<rec::Hit> > HitHandle;
    if (fWriteHits) {
        if (!e.getByLabel(fHitLabel, HitHandle)) {
            throw cet::exception("anatree") << " No rec::Hit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handle for TPCClusters
    art::Handle< std::vector<rec::TPCCluster> > TPCClusterHandle;
    if (fWriteTPCClusters) {
        if (!e.getByLabel(fTPCClusterLabel, TPCClusterHandle)) {
            throw cet::exception("anatree") << " No rec::TPCCluster branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handles for Tracks and their ionizations; also Assn's to TPCClusters, TrackIoniz
    art::Handle< std::vector<rec::Track> > TrackHandle;
    art::Handle< std::vector<rec::TrackIoniz> > TrackIonHandle;
    art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
    art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
    if(fWriteTracks) {
        if (!e.getByLabel(fTrackLabel, TrackHandle)) {
            throw cet::exception("anatree") << " No rec::Track branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (!e.getByLabel(fTrackLabel, TrackIonHandle)) {
            throw cet::exception("anatree") << " No rec::TrackIoniz branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);
        findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,e,fTrackLabel);
    }

    // Get handle for Vertices; also Assn's to Tracks
    art::Handle< std::vector<rec::Vertex> > VertexHandle;
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
    if (fWriteVertices) {
        if (!e.getByLabel(fVertexLabel, VertexHandle)) {
            throw cet::exception("anatree") << " No rec::Vertex branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);
    }
    
    // Get handle for CaloDigits
    art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
    if (fWriteCaloDigits) {
        if (!e.getByLabel(fRawCaloHitLabel, RawHitHandle)) {
            throw cet::exception("anatree") << " No :raw::CaloRawDigit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handle for CaloHits
    art::Handle< std::vector<rec::CaloHit> > RecoHitHandle;
    if (fWriteCaloHits) {
        if (!e.getByLabel(fCaloHitLabel, RecoHitHandle)) {
            throw cet::exception("anatree") << " No rec::CaloHit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    // Get handle for CaloClusters; also Assn for matching tracks
    art::Handle< std::vector<rec::Cluster> > RecoClusterHandle;
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
    if (fWriteCaloClusters) {
        if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {
            throw cet::exception("anatree") << " No rec::Cluster branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fWriteTracks) {
            findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>
                                           (RecoClusterHandle,e,fECALAssnLabel);
        }
    }



    // =============  Use art handles ==========================================
    // Pull on the handles now!
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();



    if (fWriteMCinfo) {

        // save MCTruth info    
        for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
            for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
                if (mct.NeutrinoSet()) {
                    simb::MCNeutrino nuw = mct.GetNeutrino();
                    fNeutrinoType.push_back(nuw.Nu().PdgCode());
                    fCCNC.push_back(nuw.CCNC());
                    fMode.push_back(nuw.Mode());
                    fInteractionType.push_back(nuw.InteractionType());
                    fQ2.push_back(nuw.QSqr());
                    fW.push_back(nuw.W());
                    fX.push_back(nuw.X());
                    fY.push_back(nuw.Y());
                    fTheta.push_back(nuw.Theta());
                    if (fWriteCohInfo) {
                        double getT = computeT(mct);
                        fT.push_back( static_cast<Float_t>(getT) );
                    }
                    fMCVertexX.push_back(nuw.Nu().EndX());
                    fMCVertexY.push_back(nuw.Nu().EndY());
                    fMCVertexZ.push_back(nuw.Nu().EndZ());
                    fMCnuPx.push_back(nuw.Nu().Px());
                    fMCnuPy.push_back(nuw.Nu().Py());
                    fMCnuPz.push_back(nuw.Nu().Pz());
                }  // end MC info from MCTruth
            }
        }

        // save GTruth info
        for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl) {
            for ( auto const& gt : (*gthandlelist.at(igthl)) ) {
                fGint.push_back(gt.fGint);
                fTgtPDG.push_back(gt.ftgtPDG);
                fWeight.push_back(gt.fweight);
                fgT.push_back(gt.fgT);
            }
        }

        // save MCParticle info (post-GEANT; for pre-GEANT, maybe try MCTruth?)
        for ( auto const& mcp : (*MCPHandle) ) {
            fMCPTrkID.push_back(mcp.TrackId());
            fMCPDG.push_back(mcp.PdgCode());
            // if mcp.Mother() == 0, particle is from initial vertex; MCParticles in
            // the event per se are indexed from 1 via TrackID.
            // (*MCPHandle) is a vector<simb::MCParticle> and so is indexed from 0.
            int momNumber = mcp.Mother();    int momPDG = 0;
            // Short loop to find the mother.  MCParticle.fmother appears to be the TrackID of
            // of the mother particle, minus 1.  However, not all TrackID values appear
            // sequentially in the event record; if they did you could look at the MCParticle
            // (*MCPHandle)[momNumber-1].  Too bad!  BTW, StatusCode==1 at this point, no point
            // checking that.
            if (momNumber>0) {
                for ( auto const& mcp_short : (*MCPHandle) ) {
                    if (mcp_short.TrackId() == momNumber) {
                        momPDG = mcp_short.PdgCode();
                        break;
                    }
                }
            }
            fMCMother.push_back(momNumber);
            fMCPDGMother.push_back(momPDG);

            const TLorentzVector& position = mcp.Position(0);
            const TLorentzVector& momentum = mcp.Momentum(0);
            fMCPStartX.push_back(position.X());
            fMCPStartY.push_back(position.Y());
            fMCPStartZ.push_back(position.Z());
            fMCPTime.push_back(mcp.T());
            fMCPStartPX.push_back(momentum.Px());
            fMCPStartPY.push_back(momentum.Py());
            fMCPStartPZ.push_back(momentum.Pz());

            const TLorentzVector& positionEnd = mcp.EndPosition();
            const TLorentzVector& momentumEnd = mcp.EndMomentum();
            fMCPEndX.push_back(positionEnd.X());
            fMCPEndY.push_back(positionEnd.Y());
            fMCPEndZ.push_back(positionEnd.Z());
            fMCPEndPX.push_back(momentumEnd.Px());
            fMCPEndPY.push_back(momentumEnd.Py());
            fMCPEndPZ.push_back(momentumEnd.Pz());
            fMCPProc.push_back(mcp.Process());
            fMCPEndProc.push_back(mcp.EndProcess());
        }
    }



    if(fWriteMCPTrajectory) {
        // It's in the MCParticle table
        Int_t mcpIndex = 0;
        for ( auto const& mcp : (*MCPHandle) ) {
            const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
            const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
            //No charge don't store the trajectory
            if (definition==nullptr || definition->Charge() == 0) continue;
    
            for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++) {
                fTrajMCPX.push_back(mcp.Trajectory().X(iTraj));
                fTrajMCPY.push_back(mcp.Trajectory().Y(iTraj));
                fTrajMCPZ.push_back(mcp.Trajectory().Z(iTraj));
                fTrajMCPT.push_back(mcp.Trajectory().T(iTraj));
                fTrajMCPE.push_back(mcp.Trajectory().E(iTraj));
                fTrajMCPTrajIndex.push_back(mcpIndex);
            }
            mcpIndex++;
        }
    }
    


    if (fWriteMCCaloInfo) {
        // Save simulation hit info
        for ( auto const& SimHit : (*SimHitHandle) ) {
            fSimnHits++;
            fSimHitX.push_back(SimHit.X());
            fSimHitY.push_back(SimHit.Y());
            fSimHitZ.push_back(SimHit.Z());
            fSimHitTime.push_back(SimHit.Time());
            fSimHitEnergy.push_back(SimHit.Energy());
            fSimHitTrackID.push_back(SimHit.TrackID());
            fSimHitCellID.push_back(SimHit.CellID());
            fSimEnergySum += SimHit.Energy();
        }
    }



    // save hits in the TPC
    if (fWriteHits) {
        for ( auto const& Hit : (*HitHandle) ) {
            fHitX.push_back(Hit.Position()[0]);
            fHitY.push_back(Hit.Position()[1]);
            fHitZ.push_back(Hit.Position()[2]);
            fHitSig.push_back(Hit.Signal());
            fHitRMS.push_back(Hit.RMS());
        }
    }



    // save clusters in the TPC. For some reason, can't get FindOneP<rec::Track> or
    // FindManyP<rec::Track> to work; seems the underlying Assn isn't found.  Have
    // to FindManyP<TPCCluster> instead and  iterate if (fWriteTracks).  :(
    if (fWriteTPCClusters) {
        for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
            fTPCClusterX.push_back(TPCCluster.Position()[0]);
            fTPCClusterY.push_back(TPCCluster.Position()[1]);
            fTPCClusterZ.push_back(TPCCluster.Position()[2]);
            fTPCClusterSig.push_back(TPCCluster.Signal());
            fTPCClusterRMS.push_back(TPCCluster.RMS());

            Int_t trackForThisTPCluster = -1;
            if (fWriteTracks) {
                size_t iTrack = 0;
                for ( auto const& track : (*TrackHandle) ) {
                    for (size_t iCluster=0; iCluster<track.NHits(); iCluster++) {
                        auto const& trackedCluster = 
                                  *(findManyTPCClusters->at(iTrack).at(iCluster));
                        if (TPCCluster==trackedCluster) {
                            trackForThisTPCluster = track.getIDNumber();
                            // No cluster is in 2 tracks (don't mess up dE/dx!)
                            goto pushit;   // break 2 loops
                        }
                    }
                    iTrack++;
                }
            }
            pushit:
            fTPCClusterTrkIDNumber.push_back(trackForThisTPCluster);
        }
    }



    // save per-track info
    if (fWriteTracks) {
        size_t iTrack = 0;
        for ( auto const& track : (*TrackHandle) ) {
            // track is a rec::Track, not a rec::TrackPar
            fTrackIDNumber.push_back(track.getIDNumber());

            fTrackStartX.push_back(track.Vertex()[0]);
            fTrackStartY.push_back(track.Vertex()[1]);
            fTrackStartZ.push_back(track.Vertex()[2]);
            fTrackStartPX.push_back(track.Momentum_beg()*track.VtxDir()[0]);
            fTrackStartPY.push_back(track.Momentum_beg()*track.VtxDir()[1]);
            fTrackStartPZ.push_back(track.Momentum_beg()*track.VtxDir()[2]);
            fTrackStartQ.push_back(track.ChargeBeg());

            fTrackEndX.push_back(track.End()[0]);
            fTrackEndY.push_back(track.End()[1]);
            fTrackEndZ.push_back(track.End()[2]);
            fTrackEndPX.push_back(track.Momentum_end()*track.EndDir()[0]);
            fTrackEndPY.push_back(track.Momentum_end()*track.EndDir()[1]);
            fTrackEndPZ.push_back(track.Momentum_end()*track.EndDir()[2]);
            fTrackEndQ.push_back(track.ChargeEnd());

            fTrackLenF.push_back(track.LengthForward());
            fTrackLenB.push_back(track.LengthBackward());
            fNTPCClustersOnTrack.push_back(track.NHits());

            if (findIonization->isValid()) {
                // No calibration for now.  Someday this should all be in reco
                rec::TrackIoniz ionization = *(findIonization->at(iTrack));
                float avgIonF, avgIonB;
                processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
                fTrackAvgIonF.push_back( avgIonF );
                fTrackAvgIonB.push_back( avgIonB );
            } else {
                // must push_back something so that fTrackAvgIonF,B are of correct size.
                fTrackAvgIonF.push_back( 0.0 );
                fTrackAvgIonB.push_back( 0.0 );
            }

            iTrack++;
        } // end loop over TrackHandle
    }



    // save Vertex and Track-Vertex association info
    if (fWriteVertices) {
        size_t iVertex = 0;
        for ( auto const& vertex : (*VertexHandle) ) {
            fVertexIDNumber.push_back(vertex.getIDNumber());
            fVertexX.push_back(vertex.Position()[0]);
            fVertexY.push_back(vertex.Position()[1]);
            fVertexZ.push_back(vertex.Position()[2]);
            fVertexT.push_back(vertex.Time());

            int nVertexedTracks = 0;
            if ( findManyTrackEnd->isValid() ) {
                nVertexedTracks = findManyTrackEnd->at(iVertex).size();
            }
            fVertexN.push_back(nVertexedTracks);

            int vertexCharge = 0;
            for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack) {
                fVTAssn_VertIDNumber.push_back(vertex.getIDNumber());

                // Get this vertexed track.
                rec::Track track = *(findManyTrackEnd->at(iVertex).at(iVertexedTrack));
                fVTAssn_TrackIDNumber.push_back(track.getIDNumber());

                // Get the end of the track in the vertex.  It isn't that odd for the end
                // of the track not used in the vertex to be closer to the vertex than the
                // one actually used; you might have a very short stub track in a 3 track
                // vertex with small opening angles and the other 2 tracks might pull the
                // vertex some distance towards the far end of the stub track
                rec::TrackEnd fee = *(findManyTrackEnd->data(iVertex).at(iVertexedTrack));
                // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
                fVTAssn_TrackEnd.push_back(fee);

                if (fee==rec::TrackEndBeg) {
                    vertexCharge += track.ChargeBeg();
                } else {
                    vertexCharge += track.ChargeEnd();
                }
            }
            fVertexQ.push_back(vertexCharge);
            ++iVertex;
        } // end loop over VertexHandle
    }



    // save calorimetry raw digits info
    if (fWriteCaloDigits) {
        //Save Digit Hit info
        for ( auto const& DigiHit : (*RawHitHandle) ) {
            fDiginHits++;
            fDigiHitX.push_back(DigiHit.X());
            fDigiHitY.push_back(DigiHit.Y());
            fDigiHitZ.push_back(DigiHit.Z());
            fDigiHitTime.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
            fDigiHitADC.push_back(DigiHit.ADC());
            fDigiHitCellID.push_back(DigiHit.CellID());
        }
    }



    // save reco'd Calorimetry hits
    if (fWriteCaloHits) {
        for ( auto const& Hit : (*RecoHitHandle) ) {
            fReconHitIDNumber.push_back(Hit.getIDNumber());
            fRecoHitX.push_back(Hit.Position()[0]);
            fRecoHitY.push_back(Hit.Position()[1]);
            fRecoHitZ.push_back(Hit.Position()[2]);
            fRecoHitTime.push_back(Hit.Time());
            fRecoHitEnergy.push_back(Hit.Energy());
            fRecoHitCellID.push_back(Hit.CellID());
            fRecoEnergySum += Hit.Energy();
        }
    }



    // save Cluster info
    if (fWriteCaloClusters) {
        for ( auto const& cluster : (*RecoClusterHandle) ) {
            fClusterIDNumber.push_back(cluster.getIDNumber());
            fClusterNhits.push_back(cluster.CalorimeterHits().size());
            fClusterEnergy.push_back(cluster.Energy());
            fClusterX.push_back(cluster.Position()[0]);
            fClusterY.push_back(cluster.Position()[1]);
            fClusterZ.push_back(cluster.Position()[2]);
            fClusterTheta.push_back(cluster.ITheta());
            fClusterPhi.push_back(cluster.IPhi());
            fClusterPID.push_back(cluster.ParticleID());
            // fClusterShape.push_back(cluster.Shape());
            fClusterMainAxisX.push_back(cluster.EigenVectors()[0]);
            fClusterMainAxisY.push_back(cluster.EigenVectors()[1]);
            fClusterMainAxisZ.push_back(cluster.EigenVectors()[2]);
        }
    }



    // Write info for ECAL-matched tracks
    if (fWriteMatchedTracks) {
        size_t iCluster = 0;
        for ( auto const& cluster : (*RecoClusterHandle) ) {
            int nCALedTracks(0);
            if ( findManyCALTrackEnd->isValid() ) {
                nCALedTracks = findManyCALTrackEnd->at(iCluster).size();
            }
            for (int iCALedTrack=0; iCALedTrack<nCALedTracks; ++iCALedTrack) {
                fCALAssn_ClusIDNumber.push_back(cluster.getIDNumber());
                rec::Track track  = *(findManyCALTrackEnd->at(iCluster).at(iCALedTrack));
                fCALAssn_TrackIDNumber.push_back( track.getIDNumber() );

                rec::TrackEnd fee = *(findManyCALTrackEnd->data(iCluster).at(iCALedTrack));
                fCALAssn_TrackEnd.push_back(fee);    // The rec::TrackEnd (see Track.h) that extrapolated to cluster
            }
            iCluster++;
        }



    } // end branch on fWriteCaloInfo
} // end :anatree::FillVectors



//==============================================================================
//==============================================================================
//==============================================================================
// Process ionization.  Eventually this moves into the reco code.
void gar::anatree::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
                           float& forwardIonVal, float& backwardIonVal) {

    // NO CALIBRATION SERVICE FOR NOW

    std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
    forwardIonVal = processOneDirection(SigData, ionizeTruncate);

    SigData = ion.getBAK_dSigdXs();
    backwardIonVal = processOneDirection(SigData, ionizeTruncate);

    return;
}



float gar::anatree::processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

    std::vector<std::pair<float,float>> dEvsX;    // Will be the ionization vs distance along tra

    // The first hit on the track never had its ionization info stored.  Not a problem
    // really.  Each pair is a hit and the step along the track that ends at the hit
    // For the last hit, just take the step from the n-1 hit; don't guess some distance to
    // (nonexistant!) n+1 hit.  Using pointer arithmetic because you are a real K&R C nerd!
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



//==============================================================================
//==============================================================================
//==============================================================================
// Coherent pion analysis specific code
float gar::anatree::computeT( simb::MCTruth theMCTruth ) {
    // Warning.  You probably want the absolute value of t, not t.
    int nPart = theMCTruth.NParticles();
    enum { nu, mu, pi};
    float E[3], Px[3], Py[3], Pz[3];
    E[nu] = E[mu] = E[pi] = -1e42;

    for (int i=0; i<3;++i) {
        Px[i] = 0; 
        Py[i] = 0;
        Pz[i] = 0;
        E[i]  = 0;
    }
    // Find t from the MCParticles via the
    for (int iPart=0; iPart<nPart; iPart++) {
        simb::MCParticle Part = theMCTruth.GetParticle(iPart);
        int code = Part.PdgCode();
        int mom  = Part.Mother();

        // get the neutrino
        if ( abs(code) == 12 || abs(code) == 14 || abs(code) == 16 ) {
            if (mom == -1) {
                E[nu] = Part.E();   Px[nu] = Part.Px();   Py[nu] = Part.Py();   Pz[nu] = Part.Pz();
            }
        }

        // get the lepton
        if ( abs(code) == 11 || abs(code) == 13 || abs(code) == 15 ) {
            if (mom == 0) {
                E[mu] = Part.E();   Px[mu] = Part.Px();   Py[mu] = Part.Py();   Pz[mu] = Part.Pz();
            }
        }

        // get the pion
        if ( code==111 || abs(code)==211 ) {
            if (mom == 1) {
                E[pi] = Part.E();   Px[pi] = Part.Px();   Py[pi] = Part.Py();   Pz[pi] = Part.Pz();
            }
        }

        // get outa here
        if ( E[nu]!=0 && E[mu]!=0 && E[pi]!=0) break;

    }

    // Compute t; reuse nu 4-vector to get first q, then t.
    E[nu] -= E[mu];   Px[nu] -= Px[mu];   Py[nu] -= Py[mu];   Pz[nu] -= Pz[mu];
    E[nu] -= E[pi];   Px[nu] -= Px[pi];   Py[nu] -= Py[pi];   Pz[nu] -= Pz[pi];
    float t = E[nu]*E[nu] -Px[nu]*Px[nu] -Py[nu]*Py[nu] -Pz[nu]*Pz[nu];
    return t;
}


DEFINE_ART_MODULE(gar::anatree)
