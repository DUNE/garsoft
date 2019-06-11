////////////////////////////////////////////////////////////////////////
// Class:       anatree
// Plugin Type: analyzer (art v2_11_02)
// File:        anatree_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Additions from Leo Bellantoni, 2019
////////////////////////////////////////////////////////////////////////

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

#include "Ana/anautil.h"

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

        // Truncation parameter for dE/dx (average this fraction lowest readings)
        float fIonizTruncate;

        // Optionally keep/drop parts of the analysis tree
        bool fWriteMCinfo;               ///< Info from MCTruth, GTruth into tree.  Default=true
        bool fWriteMCPTrajectory;        ///< MCP Trajectory                        Default=false

        bool fWriteHits;                 ///< Write info about Hits into tree.      Default=false
        bool fWriteTPCClusters;          ///< Write TPCClusters info into tree.     Default=true
        bool fWriteAllTracks;            ///< Start/end X, P for all tracks         Default=true
        bool fWriteTPCClustersInTracks;  ///< TPCClusters in tracks with Assns      Default=false
        bool fWriteVertsNtracks;         ///< Reco vertexes & associated tracks     Default=true
        bool fWriteCaloDigits;           ///< Low-level variables for calorimetry.  Default=false
        bool fWriteCaloInfo;             ///< Write ECAL hits & clusters.           Default=true
        bool fWriteCohInfo;              ///< MC level variables for coherent pi+.  Default=false

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

        // track data
        std::vector<Float_t>     fTrackStartX;
        std::vector<Float_t>     fTrackStartY;
        std::vector<Float_t>     fTrackStartZ;
        std::vector<Float_t>     fTrackStartPX;
        std::vector<Float_t>     fTrackStartPY;
        std::vector<Float_t>     fTrackStartPZ;

        std::vector<Float_t>     fTrackEndX;
        std::vector<Float_t>     fTrackEndY;
        std::vector<Float_t>     fTrackEndZ;
        std::vector<Float_t>     fTrackEndPX;
        std::vector<Float_t>     fTrackEndPY;
        std::vector<Float_t>     fTrackEndPZ;

        std::vector<Float_t>     fTrackLenF;
        std::vector<Float_t>     fTrackLenB;
        std::vector<Float_t>     fTrackAvgIon;
        std::vector<Int_t>       fNTPCClustersOnTrack;

        // TPCClusters belonging to tracks data
        std::vector<Float_t>     fTrkTPCClusterX;
        std::vector<Float_t>     fTrkTPCClusterY;
        std::vector<Float_t>     fTrkTPCClusterZ;
        std::vector<Float_t>     fTrkTPCClusterSig;
        std::vector<Float_t>     fTrkTPCClusterRMS;
        std::vector<Int_t>       fTrkTPCClusterTrkIndex;

        // vertex branches
        std::vector<Float_t>     fVertexX;
        std::vector<Float_t>     fVertexY;
        std::vector<Float_t>     fVertexZ;
        std::vector<ULong64_t>   fVertexT;
        std::vector<Int_t>       fVertexN;
        std::vector<Int_t>       fVertexQ;

        // vertex-track Assn branches
        std::vector<Int_t>       fVTAssn_Vertex;
        std::vector<Float_t>     fVTAssn_TrkEndPx;
        std::vector<Float_t>     fVTAssn_TrkEndPy;
        std::vector<Float_t>     fVTAssn_TrkEndPz;
        std::vector<Float_t>     fVTAssn_TrkEndLen;
        std::vector<Float_t>     fVTAssn_TrkEndChiSq;
        std::vector<Float_t>     fVTAssn_TrkEndIon;
        std::vector<Int_t>       fVTAssn_TrkEndNclus;

        // // raw calo hit data
        UInt_t                   fDiginHits;
        std::vector<Float_t>     fDigiHitX;
        std::vector<Float_t>     fDigiHitY;
        std::vector<Float_t>     fDigiHitZ;
        std::vector<Float_t>     fDigiHitTime;
        std::vector<UInt_t>      fDigiHitADC;
        std::vector<ULong64_t>   fDigiHitCellID;

        // reco calo hit data
        UInt_t                   fReconHits;
        std::vector<Float_t>     fRecoHitX;
        std::vector<Float_t>     fRecoHitY;
        std::vector<Float_t>     fRecoHitZ;
        std::vector<Float_t>     fRecoHitTime;
        std::vector<Float_t>     fRecoHitEnergy;
        std::vector<ULong64_t>   fRecoHitCellID;
        Float_t                  fRecoEnergySum;

        // calo cluster data
        UInt_t                   fnCluster;
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
    };
}



//==============================================================================
// constructor
gar::anatree::anatree(fhicl::ParameterSet const & p) : EDAnalyzer(p) {
    bool usegenlabels  = p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
    if (!usegenlabels) fGeneratorLabels.clear();
    bool usegeniegenlabels  = p.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
    if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

    fGeantLabel       = p.get<std::string>("GEANTLabel","geant");
    fHitLabel         = p.get<std::string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<std::string>("TrackLabel","track");
    fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
    fRawCaloHitLabel  = p.get<std::string>("RawCaloHitLabel","daqecal");
    fCaloHitLabel     = p.get<std::string>("CaloHitLabel","calohit");
    fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");

    //MC Truth
    fWriteMCinfo              = p.get<bool>("WriteMCinfo",true);
    fWriteMCPTrajectory       = p.get<bool>("WriteMCPTrajectory",false);

    fWriteHits                = p.get<bool>("WriteHits",false);
    fWriteTPCClusters         = p.get<bool>("WriteTPCClusters",true);
    fWriteAllTracks           = p.get<bool>("WriteAllTracks",true);
    fWriteTPCClustersInTracks = p.get<bool>("WriteTPCClustersInTracks",true);
    fWriteVertsNtracks        = p.get<bool>("WriteVertsNtracks",true);
    fIonizTruncate            = p.get<float>("IonizTruncate",0.70);
    fWriteCaloDigits          = p.get<bool>("WriteCaloInfo",false);
    fWriteCaloInfo            = p.get<bool>("WriteCaloInfo",true);
    fWriteCohInfo             = p.get<bool>("WriteCohInfo",false);

    consumes<std::vector<simb::MCParticle> >(fGeantLabel);

    if (usegenlabels)
    {
        for (size_t i=0; i<fGeneratorLabels.size(); ++i)
        {
            consumes<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
        }
    }
    else
    {
        consumesMany<std::vector<simb::MCTruth> >();
    }

    if (usegeniegenlabels)
    {
        for (size_t i=0; i<fGENIEGeneratorLabels.size(); ++i)
        {
            consumes<std::vector<simb::GTruth> >(fGENIEGeneratorLabels.at(i));
        }
    }
    else
    {
        consumesMany<std::vector<simb::GTruth> >();
    }

    //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);

    consumes<std::vector<gar::rec::TPCCluster> >(fTPCClusterLabel);
    consumes<art::Assns<gar::rec::Track, gar::rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<gar::rec::Hit> >(fHitLabel);
    consumes<std::vector<gar::rec::Track> >(fTrackLabel);
    consumes<std::vector<gar::rec::Vertex> >(fVertexLabel);
    consumes<art::Assns<gar::rec::Track, gar::rec::Vertex> >(fVertexLabel);
    consumes<std::vector<gar::rec::TrackIoniz>>(fTrackLabel);
    consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

    consumes<std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
    consumes<std::vector<gar::raw::CaloRawDigit> >(fRawCaloHitLabel);
    consumes<std::vector<gar::rec::CaloHit> >(fCaloHitLabel);
    consumes<std::vector<gar::rec::Cluster> >(fClusterLabel);
} // end constructor



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

        if (fWriteMCPTrajectory) {
            //Write MCP Trajectory
            fTree->Branch("TrajMCPX",          &fTrajMCPX);
            fTree->Branch("TrajMCPY",          &fTrajMCPY);
            fTree->Branch("TrajMCPZ",          &fTrajMCPZ);
            fTree->Branch("TrajMCPT",          &fTrajMCPT);
            fTree->Branch("TrajMCPE",          &fTrajMCPE);
            fTree->Branch("TrajMCPTrajIndex",  &fTrajMCPTrajIndex);
        }

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

    if (fWriteHits) {
        fTree->Branch("HitX",        &fHitX);
        fTree->Branch("HitY",        &fHitY);
        fTree->Branch("HitZ",        &fHitZ);
        fTree->Branch("HitSig",      &fHitSig);
        fTree->Branch("HitRMS",      &fHitRMS);
    }

    if (fWriteTPCClusters) {
        fTree->Branch("TPCClusterX",        &fTPCClusterX);
        fTree->Branch("TPCClusterY",        &fTPCClusterY);
        fTree->Branch("TPCClusterZ",        &fTPCClusterZ);
        fTree->Branch("TPCClusterSig",      &fTPCClusterSig);
        fTree->Branch("TPCClusterRMS",      &fTPCClusterRMS);
    }

    if (fWriteAllTracks) {                         // All position, momentum, etc
        fTree->Branch("TrackStartX",     &fTrackStartX);
        fTree->Branch("TrackStartY",     &fTrackStartY);
        fTree->Branch("TrackStartZ",     &fTrackStartZ);
        fTree->Branch("TrackStartPX",    &fTrackStartPX);
        fTree->Branch("TrackStartPY",    &fTrackStartPY);
        fTree->Branch("TrackStartPZ",    &fTrackStartPZ);

        fTree->Branch("TrackEndX",       &fTrackEndX);
        fTree->Branch("TrackEndY",       &fTrackEndY);
        fTree->Branch("TrackEndZ",       &fTrackEndZ);
        fTree->Branch("TrackEndPX",      &fTrackEndPX);
        fTree->Branch("TrackEndPY",      &fTrackEndPY);
        fTree->Branch("TrackEndPZ",      &fTrackEndPZ);

        fTree->Branch("TrackLenF",       &fTrackLenF);
        fTree->Branch("TrackLenB",       &fTrackLenB);
        fTree->Branch("TrackAvgIon",     &fTrackAvgIon);
        fTree->Branch("NTPCClustersOnTrack",    &fNTPCClustersOnTrack);
    }

    if (fWriteTPCClustersInTracks) {               // Write TPCClusters in reco tracks
        fTree->Branch("TrkTPCClusterX",        &fTrkTPCClusterX);
        fTree->Branch("TrkTPCClusterY",        &fTrkTPCClusterY);
        fTree->Branch("TrkTPCClusterZ",        &fTrkTPCClusterZ);
        fTree->Branch("TrkTPCClusterSig",      &fTrkTPCClusterSig);
        fTree->Branch("TrkTPCClusterRMS",      &fTrkTPCClusterRMS);
        fTree->Branch("TrkTPCClusterTrkIndex", &fTrkTPCClusterTrkIndex);
    }

    if (fWriteVertsNtracks) {                     // Reco'd verts & their track-ends
        fTree->Branch("VertX",           &fVertexX);
        fTree->Branch("VertY",           &fVertexY);
        fTree->Branch("VertZ",           &fVertexZ);
        fTree->Branch("VertT",           &fVertexT);
        fTree->Branch("VertN",           &fVertexN);
        fTree->Branch("VertQ",           &fVertexQ);

        fTree->Branch("VT_Vertex",       &fVTAssn_Vertex);
        fTree->Branch("VT_TrkEndPx",     &fVTAssn_TrkEndPx);
        fTree->Branch("VT_TrkEndPy",     &fVTAssn_TrkEndPy);
        fTree->Branch("VT_TrkEndPz",     &fVTAssn_TrkEndPz);
        fTree->Branch("VT_TrkEndLen",    &fVTAssn_TrkEndLen);
        fTree->Branch("VT_TrkEndChiSq",  &fVTAssn_TrkEndChiSq);
        fTree->Branch("VT_TrkEndIon",    &fVTAssn_TrkEndIon);
        fTree->Branch("VT_TrkEndNclus",  &fVTAssn_TrkEndNclus);
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

    if (fWriteCaloInfo) {                         // Write calorimetry hits & clusters
        fTree->Branch("ReconHits",        &fReconHits);
        fTree->Branch("RecoHitX",         &fRecoHitX);
        fTree->Branch("RecoHitY",         &fRecoHitY);
        fTree->Branch("RecoHitZ",         &fRecoHitZ);
        fTree->Branch("RecoHitTime",      &fRecoHitTime);
        fTree->Branch("RecoHitEnergy",    &fRecoHitEnergy);
        fTree->Branch("RecoHitCellID",    &fRecoHitCellID);
        fTree->Branch("RecoEnergySum",    &fRecoEnergySum);

        fTree->Branch("nCluster",         &fnCluster);
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
}  // End of gar::anatree::beginJob



//==============================================================================
void gar::anatree::analyze(art::Event const & e)
{
    // clear out all our vectors
    gar::anatree::ClearVectors();

    // Fill the vectors
    gar::anatree::FillVectors(e);

    fTree->Fill();
}




//==========================================================================
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

        if (fWriteMCPTrajectory) {
            fTrajMCPX.clear();
            fTrajMCPY.clear();
            fTrajMCPZ.clear();
            fTrajMCPT.clear();
            fTrajMCPE.clear();
            fTrajMCPTrajIndex.clear();
        }

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
    }

    if (fWriteAllTracks) {
        fTrackStartX.clear();
        fTrackStartY.clear();
        fTrackStartZ.clear();
        fTrackStartPX.clear();
        fTrackStartPY.clear();
        fTrackStartPZ.clear();
        fTrackEndX.clear();
        fTrackEndY.clear();
        fTrackEndZ.clear();
        fTrackEndPX.clear();
        fTrackEndPY.clear();
        fTrackEndPZ.clear();
        fTrackLenF.clear();
        fTrackLenB.clear();
        fTrackAvgIon.clear();
        fNTPCClustersOnTrack.clear();
    }

    if (fWriteTPCClustersInTracks) {
        fTrkTPCClusterX.clear();
        fTrkTPCClusterY.clear();
        fTrkTPCClusterZ.clear();
        fTrkTPCClusterSig.clear();
        fTrkTPCClusterRMS.clear();
        fTrkTPCClusterTrkIndex.clear();
    }

    if (fWriteVertsNtracks) {
        fVertexX.clear();
        fVertexY.clear();
        fVertexZ.clear();
        fVertexT.clear();
        fVertexN.clear();
        fVertexQ.clear();

        fVTAssn_Vertex.clear();
        fVTAssn_TrkEndPx.clear();
        fVTAssn_TrkEndPy.clear();
        fVTAssn_TrkEndPz.clear();
        fVTAssn_TrkEndLen.clear();
        fVTAssn_TrkEndChiSq.clear();
        fVTAssn_TrkEndIon.clear();
        fVTAssn_TrkEndNclus.clear();
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

    if (fWriteCaloInfo) {
        fReconHits = 0;
        fRecoHitX.clear();
        fRecoHitY.clear();
        fRecoHitZ.clear();
        fRecoHitTime.clear();
        fRecoHitEnergy.clear();
        fRecoHitCellID.clear();
        fRecoEnergySum = 0.;

        fnCluster = 0;
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
} // end gar::anatree::ClearVectors



//==========================================================================
void gar::anatree::FillVectors(art::Event const & e) {

    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();

    // Get a good grip on that data!

    std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
    std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;

    art::Handle< std::vector<simb::MCParticle> > MCPHandle;
    art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;

    if (fWriteMCinfo) {
    // Handles for MC info might be there, might not
    if (fGeneratorLabels.size()<1)
    {
        // get them all (even if there are none)
        e.getManyByType(mcthandlelist);
    }
    else
    {
        mcthandlelist.resize(fGeneratorLabels.size());
        for (size_t i=0; i< fGeneratorLabels.size(); ++i)
        {
            // complain if we wanted a specific one but didn't find it
            if (!e.getByLabel(fGeneratorLabels.at(i),mcthandlelist.at(i)))
            throw cet::exception("anatree") << " No simb::MCTruth branch - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }

    if (fGENIEGeneratorLabels.size()<1)
    {
        e.getManyByType(gthandlelist);  // get them all (even if there are none)
    }
    else
    {
        gthandlelist.resize(fGENIEGeneratorLabels.size());
        for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i)
        {
            // complain if we wanted a specific one but didn't find it
            if (!e.getByLabel(fGENIEGeneratorLabels.at(i),gthandlelist.at(i)))
            throw cet::exception("anatree")
            << " No simb::GTruth branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;

        }
    }

    if (!e.getByLabel(fGeantLabel, MCPHandle)) {
        throw cet::exception("anatree")
        << " No simb::MCParticle branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (!e.getByLabel(fGeantLabel, SimHitHandle)) {
        throw cet::exception("anatree")
        << " No gar::sdp::CaloDeposit branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}

art::Handle< std::vector<gar::rec::Hit> > HitHandle;
art::Handle< std::vector<gar::rec::TPCCluster> > TPCClusterHandle;
art::Handle< std::vector<gar::rec::Track> > TrackHandle;
art::Handle< std::vector<gar::rec::TrackIoniz> > TrackIonHandle;
art::Handle< std::vector<gar::rec::Vertex> > VertexHandle;
if (fWriteHits) {
    if (!e.getByLabel(fHitLabel, HitHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::Hit branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}

if (fWriteTPCClusters) {
    if (!e.getByLabel(fTPCClusterLabel, TPCClusterHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::TPCCluster branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}

if(fWriteAllTracks)
{
    if (!e.getByLabel(fTrackLabel, TrackHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::Track branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (!e.getByLabel(fTrackLabel, TrackIonHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::TrackIoniz branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}

if(fWriteVertsNtracks)
{

    if (!e.getByLabel(fVertexLabel, VertexHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::Vertex branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}

art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
art::Handle< std::vector<gar::rec::CaloHit> >      RecoHitHandle;
art::Handle< std::vector<gar::rec::Cluster> >      RecoClusterHandle;
if (fWriteCaloInfo) {
    if (!e.getByLabel(fRawCaloHitLabel, RawHitHandle)) {
        throw cet::exception("anatree")
        << " No gar::raw::CaloRawDigit branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (!e.getByLabel(fCaloHitLabel, RecoHitHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::CaloHit branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if (!e.getByLabel(fClusterLabel, RecoClusterHandle)) {
        throw cet::exception("anatree")
        << " No gar::rec::Cluster branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
}



// Pull on the handles now!
if (fWriteMCinfo) {
    // save MCTruth info

    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl)
    {
        for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
            if (mct.NeutrinoSet())
            {

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
                    double getT = gar::computeT(mct);
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
    for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl)
    {
        for ( auto const& gt : (*gthandlelist.at(igthl)) ) {
            fGint.push_back(gt.fGint);
            fTgtPDG.push_back(gt.ftgtPDG);
            fWeight.push_back(gt.fweight);
            fgT.push_back(gt.fgT);
        }
    }

    // save MCParticle info (post-GEANT; for pre-GEANT, maybe try MCTruth?)
    Int_t mcpIndex = 0;
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

        if(fWriteMCPTrajectory) {
            const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
            const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
            //No charge don't store the trajectory
            if(nullptr == definition || definition->Charge() == 0) continue;

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
    }  // end MC info from MCParticle

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
} // end if MCWriteInfo

// save TPC Hit info
if (fWriteHits) {
    for ( auto const& Hit : (*HitHandle) ) {
        fHitX.push_back(Hit.Position()[0]);
        fHitY.push_back(Hit.Position()[1]);
        fHitZ.push_back(Hit.Position()[2]);
        fHitSig.push_back(Hit.Signal());
        fHitRMS.push_back(Hit.RMS());
    }
}

// save TPCCluster info
if (fWriteTPCClusters) {
    for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
        fTPCClusterX.push_back(TPCCluster.Position()[0]);
        fTPCClusterY.push_back(TPCCluster.Position()[1]);
        fTPCClusterZ.push_back(TPCCluster.Position()[2]);
        fTPCClusterSig.push_back(TPCCluster.Signal());
        fTPCClusterRMS.push_back(TPCCluster.RMS());
    }
}

// save Track info.  Need track-TPCCluster-Ionization Assns outside if () for later use
if (fWriteAllTracks) {
    const art::FindManyP<gar::rec::TPCCluster> findManyTPCClusters(TrackHandle,e,fTrackLabel);
    const art::FindOneP<gar::rec::TrackIoniz> findIonization(TrackHandle,e,fTrackLabel);

    size_t iTrack = 0;
    for ( auto const& track : (*TrackHandle) ) {
        // track is a gar::rec::Track, not a gar::rec::TrackPar
        fTrackStartX.push_back(track.Vertex()[0]);
        fTrackStartY.push_back(track.Vertex()[1]);
        fTrackStartZ.push_back(track.Vertex()[2]);
        fTrackStartPX.push_back(track.Momentum_beg()*track.VtxDir()[0]);
        fTrackStartPY.push_back(track.Momentum_beg()*track.VtxDir()[1]);
        fTrackStartPZ.push_back(track.Momentum_beg()*track.VtxDir()[2]);

        fTrackEndX.push_back(track.End()[0]);
        fTrackEndY.push_back(track.End()[1]);
        fTrackEndZ.push_back(track.End()[2]);
        fTrackEndPX.push_back(track.Momentum_end()*track.EndDir()[0]);
        fTrackEndPY.push_back(track.Momentum_end()*track.EndDir()[1]);
        fTrackEndPZ.push_back(track.Momentum_end()*track.EndDir()[2]);

        fTrackLenF.push_back(track.LengthForward());
        fTrackLenB.push_back(track.LengthBackward());

        if (findIonization.isValid()) {
            // No calibration for now
            rec::TrackIoniz ionization = *(findIonization.at(iTrack));
            fTrackAvgIon.push_back( gar::processIonizationInfo(ionization, fIonizTruncate) );
        } else {
            // must push_back something so that fTrackAvgIon is of correct size.
            fTrackAvgIon.push_back( 0.0 );
        }

        int nTrackedTPCClusters = 0;
        if (findManyTPCClusters.isValid()) {
            // TPCClusters is a vector of gar::rec::TPCCluster
            auto const& TPCClusters = findManyTPCClusters.at(iTrack);
            nTrackedTPCClusters = TPCClusters.size();
        }
        fNTPCClustersOnTrack.push_back(nTrackedTPCClusters);

        if (fWriteTPCClustersInTracks) {
            int iTrackedTPCCluster=0;
            for (; iTrackedTPCCluster<nTrackedTPCClusters; iTrackedTPCCluster++) {
                auto const& TPCCluster = *(findManyTPCClusters.at(iTrack).at(iTrackedTPCCluster));
                fTrkTPCClusterX.push_back(TPCCluster.Position()[0]);
                fTrkTPCClusterY.push_back(TPCCluster.Position()[1]);
                fTrkTPCClusterZ.push_back(TPCCluster.Position()[2]);
                fTrkTPCClusterSig.push_back(TPCCluster.Signal());
                fTrkTPCClusterRMS.push_back(TPCCluster.RMS());
                fTrkTPCClusterTrkIndex.push_back((Int_t)iTrack);
            }
        }
        iTrack++;
    } // end loop over TrackHandle
}

// save Vertex and Track-Vertex association info
if (fWriteVertsNtracks) {
    const art::FindManyP<gar::rec::TPCCluster> findManyTPCClusters(TrackHandle,e,fTrackLabel);
    const art::FindOneP<gar::rec::TrackIoniz> findIonization(TrackHandle,e,fTrackLabel);
    const art::FindMany<gar::rec::Track, gar::rec::TrackEnd> findManyTrack(VertexHandle,e,fVertexLabel);
    size_t iVertex = 0;
    for ( auto const& vertex : (*VertexHandle) ) {
        fVertexX.push_back(vertex.Position()[0]);
        fVertexY.push_back(vertex.Position()[1]);
        fVertexZ.push_back(vertex.Position()[2]);
        fVertexT.push_back(vertex.Time());

        int nVertexedTracks = 0;
        if ( findManyTrack.isValid() ) {
            // tracks is a vector of gar::rec::Track
            auto const& tracks = findManyTrack.at(iVertex);
            nVertexedTracks = tracks.size();
        }
        fVertexN.push_back(nVertexedTracks);

        int vertexCharge = 0;
        for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack) {
            fVTAssn_Vertex.push_back(iVertex);  // 1 push of iVertex for each track in vertex

            // Get this vertexed track.  findManyTrack.at(iVertex) is a
            // vector<const gar::rec::Track*, std::allocator<const gar::rec::Track*> >
            gar::rec::Track track = *(findManyTrack.at(iVertex).at(iVertexedTrack));

            // Get the end of the track in the vertex.  It isn't that odd for the end
            // of the track not used in the vertex to be closer to the vertex than the
            // one actually used; you might have a very short stub track in a 3 track
            // vertex with small opening angles and the other 2 tracks might pull the
            // vertex some distance towards the far end of the stub track
            gar::rec::TrackEnd fee = *(findManyTrack.data(iVertex).at(iVertexedTrack));

            if (fee==gar::rec::TrackEndBeg) {
                fVTAssn_TrkEndPx.push_back( track.Momentum_beg()*track.VtxDir()[0] );
                fVTAssn_TrkEndPy.push_back( track.Momentum_beg()*track.VtxDir()[1] );
                fVTAssn_TrkEndPz.push_back( track.Momentum_beg()*track.VtxDir()[2] );
                fVTAssn_TrkEndLen.push_back( track.LengthForward() );
                fVTAssn_TrkEndChiSq.push_back( track.ChisqForward() );
                vertexCharge += track.ChargeBeg();
            } else {
                fVTAssn_TrkEndPx.push_back( track.Momentum_end()*track.EndDir()[0] );
                fVTAssn_TrkEndPy.push_back( track.Momentum_end()*track.EndDir()[1] );
                fVTAssn_TrkEndPz.push_back( track.Momentum_end()*track.EndDir()[2] );
                fVTAssn_TrkEndLen.push_back( track.LengthBackward() );
                fVTAssn_TrkEndChiSq.push_back( track.ChisqBackward() );
                vertexCharge += track.ChargeEnd();
            }

            // Stupid kludge to find index of this track, after all.  *yeesh*
            // We have assns from Vert to Track and Track to Ionize (and TPCClusters).  At this
            // point in the code we have a Track from a Vert of interest and want the Ionize for
            // it.  But FindOne, FindMany, FindOneP and FindManyP can't be indexed by the values
            // they store; the at() and get() and for that matter, the data() methods only have
            // size_type indexing.  A better design would have the TPCCluster & Ionize data in the
            // track rather than associated to it; or provide an equality operator in Track, possibly
            // by incorporating an identifier number in the Track class.
            int iTrack = 0;   bool found = false;
            float X =track.Vertex()[0];            // A single floating point == should do it
            for ( auto const& soughtTrack : (*TrackHandle) ) {
                if ( X==soughtTrack.Vertex()[0] ) {
                    found = true;
                    break;
                }
                ++iTrack;
            }
            if (!found) {
                throw cet::exception("anatree")
                << " Fail to find vertexted track in TrackHandle "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }

            if (found && findIonization.isValid()) {
                // No calibration for now
                rec::TrackIoniz ionization = *(findIonization.at(iTrack));
                fVTAssn_TrkEndIon.push_back( gar::processIonizationInfo(ionization, fIonizTruncate) );
            } else {
                fVTAssn_TrkEndIon.push_back( 0.0 );
            }

            if (found && findManyTPCClusters.isValid()) {
                // TPCClusters is a vector of gar::rec::TPCCluster
                auto const& TPCClusters = findManyTPCClusters.at(iTrack);
                fVTAssn_TrkEndNclus.push_back( TPCClusters.size() );
            } else {
                fVTAssn_TrkEndNclus.push_back( 0 );
            }
        }

        fVertexQ.push_back(vertexCharge);
        ++iVertex;
    } // end loop over VertexHandle
}

// save calorimetry info
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

if (fWriteCaloInfo) {
    //Save Reco Hit info
    for ( auto const& Hit : (*RecoHitHandle) ) {
        fReconHits++;
        fRecoHitX.push_back(Hit.Position()[0]);
        fRecoHitY.push_back(Hit.Position()[1]);
        fRecoHitZ.push_back(Hit.Position()[2]);
        fRecoHitTime.push_back(Hit.Time());
        fRecoHitEnergy.push_back(Hit.Energy());
        fRecoHitCellID.push_back(Hit.CellID());
        fRecoEnergySum += Hit.Energy();
    }

    // save Cluster info
    for ( auto const& Cluster : (*RecoClusterHandle) ) {
        fnCluster++;
        fClusterNhits.push_back(Cluster.CalorimeterHits().size());
        fClusterEnergy.push_back(Cluster.Energy());
        fClusterX.push_back(Cluster.Position()[0]);
        fClusterY.push_back(Cluster.Position()[1]);
        fClusterZ.push_back(Cluster.Position()[2]);
        fClusterTheta.push_back(Cluster.ITheta());
        fClusterPhi.push_back(Cluster.IPhi());
        fClusterPID.push_back(Cluster.ParticleID());
        // fClusterShape.push_back(Cluster.Shape());
        fClusterMainAxisX.push_back(Cluster.EigenVectors()[0]);
        fClusterMainAxisY.push_back(Cluster.EigenVectors()[1]);
        fClusterMainAxisZ.push_back(Cluster.EigenVectors()[2]);
    }
} // end branch on fWriteCaloInfo
} // end gar::anatree::FillVectors



DEFINE_ART_MODULE(gar::anatree)
