////////////////////////////////////////////////////////////////////////////////
// Class:       
// Plugin Type: analyzer (art v2_11_02)
// File:        StructuredTree_module.cc
//
// Copied from anatree_module.cc on 2020 Dec 18 by C. Hilgenberg (chilgenb@fnal.gov)
// This module extracts info useful for general analysis or reco R&D based
// on the fcl config. Several class based trees are writen to file.
//
////////////////////////////////////////////////////////////////////////////////

// framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//simb::
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

//garsoft includes
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/Vee.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "SummaryDataProducts/RunData.h"
#include "SummaryDataProducts/POTSummary.h"

#include "AnaUtils.cxx"

#include "MCCheater/BackTracker.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH2.h"
#include "TFile.h"

#include "CLHEP/Random/RandFlat.h"

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>

//save ourselves some typing
using std::vector;
using std::string;
using std::pair;
using simb::MCParticle;
using art::InputTag;
using art::Assns;
using art::Event;
using art::Handle;
using art::ValidHandle;

namespace gar {

    typedef pair<float, string> P;

    class StructuredTree : public art::EDAnalyzer {
    public:
        explicit StructuredTree(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        StructuredTree(StructuredTree const &) = delete;
        StructuredTree(StructuredTree &&) = delete;
        StructuredTree & operator = (StructuredTree const &) = delete;
        StructuredTree & operator = (StructuredTree &&) = delete;

        virtual void beginJob() override;
	virtual void endRun(art::Run const& run) override;
        virtual void endJob() override;

        // Required functions.
        void analyze( Event const & e ) override;

    private:
        void ClearVectors();///< clear all vectors and reset counters

        //Working Horses
        void FillGenTree( Event const & e ); ///< extracts generator(truth) level sim products
        void FillG4Tree( Event const & e );  ///< extracts G4(truth) level sim products
        void FillDetTree( Event const & e ); ///< extracts readout sim products
        void FillRecoInfo( Event const & e ); ///< extracts base level reco products
        void FillHighLevelRecoInfo( Event const & e ); ///< extracts high level reco products

        //Helpers
        std::string PointToRegion(const TVector3& point);
        void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
        float& forwardIonVal, float& backwardIonVal);
        float processOneDirection(vector<pair<float,float>> SigData,
        float ionizeTruncate);
        // Helper method for processOneDirection
        static bool lessThan_byE(std::pair<float,float> a, std::pair<float,float> b)
        {return a.first < b.first;}

        // Calculate the track PID based on reco momentum and Tom's parametrization
        vector< pair<int, float> > processPIDInfo( float p );

        float xTPC;
        float rTPC;

        // Input data labels
        string fPOTtag;
        string fGeantLabel; ///< module label for geant4 simulated hits
        string fGeantInstanceCalo; ///< Instance name ECAL for sdp::CaloDeposit
        string fGeantInstanceMuID; ///< Instance name MuID for sdp::CaloDeposit

        string fHitLabel; ///< module label for reco TPC hits rec::Hit
        string fTPCClusterLabel; ///< module label for TPC Clusters rec::TPCCluster
        string fTrackLabel; ///< module label for TPC Tracks rec:Track
        string fVertexLabel; ///< module label for vertexes rec:Vertex
        string fVeeLabel; ///< module label for conversion/decay vertexes rec:Vee

        string fTPCDigitLabel; ///< module label for digitized TPC hits raw::RawDigit
        string fCalDigitLabel; ///< module label for digitized calo/muID hits raw::CaloRawDigit
        string fCalDigitInstance; ///< Instance name ECAL for raw::CaloRawDigit
        string fMuDigitInstance; ///< Instance name MuID for raw::CaloRawDigit

        string fCaloHitLabel; ///< module label for reco calo hits rec::CaloHit
        string fCaloHitInstanceCalo; ///< Instance name ECAL for rec::CaloHit
        string fMuIDHitLabel;
        string fCaloHitInstanceMuID; ///< Instance name MuID for rec::CaloHit

        string fClusterLabel; ///< module label for calo clusters rec::Cluster
        string fClusterInstanceCalo; ///< Instance name ECAL for rec::Cluster
        string fPFLabel; ///< module label for reco particles rec::PFParticle
        string fECALAssnLabel; ///< module label for track-clusters associations

        //config
        string fAnaMode; ///< analysis configs: (gen)eral, (reco) R&D or (readout) sim R&D, default="gen"
        bool   fWriteDisplay; ///< include sim/reco traj points for event display, default=false


        // Truncation parameter for dE/dx (average this fraction lowest readings)
        float fIonizTruncate;      ///<                               Default=1.00;

        // the analysis output trees
        TTree* fHeaderTree;
        TTree* fGenTree;
        TTree* fG4Tree;
        TTree* fDetTree;
        TTree* fRecoTree;
        TTree* fDisplayTree;

        //Geometry
        const geo::GeometryCore* fGeo; ///< pointer to the geometry service

        //backtracker for determining inter product associations
        cheat::BackTrackerCore* fBt;

        // global event info
        Int_t fEvent;        ///< number of the event being processed
        Int_t fRun;          ///< number of the run being processed
        Int_t fSubRun;       ///< number of the sub-run being processed
	Int_t fPOT;
        Int_t fNSpills;
	TLorentzVector fTpcCenter;
        std::string fGeometry;      
 
        // headerTree
        string   fTreeType; ///<"structured" (in our case yes) or "flat" (no)

        // genTree
        vector<Int_t>                      fGIndex;       ///< index of GTruth object associated with MCTruth in this row (-1 if not associated, e.g. CRY event)
        vector<garana::GTruth>             fGTruth;       ///< GENIE interction-level info
        vector<vector<garana::FSParticle>> fFSParticles;  ///< FS particles/4-vecs

        // g4Tree
        vector<garana::G4Particle>         fG4Particles;  ///< 'condensed' MCParticles from G4
        vector<UInt_t>                     fG4TruthIndex; ///< index of GTruth objects matched to this particle
        vector<UInt_t>                     fG4FSIndex;    ///< index of FSParticle objects matched to this particle
        vector<vector<sdp::EnergyDeposit>> fSimTPCHits;   ///< several "hits" per MCParticle
        vector<vector<sdp::CaloDeposit>>   fSimCalHits;  
        vector<vector<sdp::CaloDeposit>>   fSimMuHits;
        //vector<UInt_t>                     fSimTPCHitsG4P;
        //vector<UInt_t>                     fSimCalHitsG4P;
        //vector<UInt_t>                     fSimMuHitsG4P;

        // detTree
        vector<raw::RawDigit>              fTPCDigits;
        vector<raw::CaloRawDigit>          fCaloDigits;
        vector<raw::CaloRawDigit>          fMuDigits;
        //TODO: add Assns digit -> simDeposit (need to add to BackTrackerCore)

        // recoTree
        //   Reco Hit data
        vector<rec::Hit>           fTPCHits;
        vector<rec::CaloHit>       fCalHits;
        vector<rec::CaloHit>       fMuHits;
        //TODO: add Assns hit -> digit (need to add to BackTrackerCore)

        //   Cluster data
        vector<rec::TPCCluster>    fTPCClusters;  ///< reco TPC clusters
        vector<garana::CaloCluster>       fCalClusters;  ///< reco ECal clusters
        vector<garana::CaloCluster>       fMuClusters;   ///< reco MuID clusters
        //TODO: add Assns Cluster->Hit

        //   Track data
        vector<garana::Track>    fTracks;           ///< reco TPC tracks
        vector<vector<UInt_t>>   fTrackG4PIndices;  ///< index of fG4Particles associated with fTracks

        // Vertex data
        vector<garana::Vertex>    fVertices;          ///< reco TPC vertices
        vector<vector<UInt_t>>    fVertTrackIndices;  ///< index of fTracks associated with fVertices
        vector<vector<Int_t>>     fVertTrackEnds;      ///< which fTrack end belongs to this vertex

        // Vee data
        vector<garana::Vee>       fVees;             ///< reco Vees (reco 4-mom, 4-pos from decay vertex)
        vector<vector<UInt_t>>    fVeeTrackIndices;  ///< index of fTracks associated with fVees
        vector<vector<Int_t>>     fVeeTrackEnds;      ///< track end associated with fVees

        // ECal cluster data
        vector<vector<UInt_t>>    fCalTrackIndices; ///< index of fTracks associated with fCalClusters
        vector<vector<Int_t>>     fCalTrackEnds;   ///< Track end associated with the cluster
        vector<vector<UInt_t>>    fCalG4PIndices; ///< index of fG4Particles associated with fCalClusters

        // displayTree
        vector<vector<TVector3>>  fMCTraj;       ///< sparsified MCTrajectory
        vector<vector<TVector3>>  fTrackTraj;    ///< reconstructed track trajectory points (using only FWD fit for now)
	vector<UInt_t>            fDispG4Index;
        vector<UInt_t>            fDispTrkIndex;

        //map of PID TH2 per momentum value
        CLHEP::HepRandomEngine &fEngine;  ///< random engine
        std::unordered_map<int, TH2F*> m_pidinterp;

        //data members used for associations bookkeeping between fill methods
        vector<const simb::MCTruth*> fMCTruths; //needed for G4Particle -> MCTruth (FS particles) matching
        vector<const simb::MCParticle*> fMCParticles; //needed for reco object -> MCParticle matching

	std::map<std::string,int> fRegionNameToID;

    };//class
}//namespace gar

//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::StructuredTree::StructuredTree(fhicl::ParameterSet const & p)
: EDAnalyzer(p),
  fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed"))
{
    fGeo     = providerFrom<geo::GeometryGAr>();

    //Sim Hits
    fGeantLabel        = p.get<string>("GEANTLabel","geant");
    fGeantInstanceCalo = p.get<string>("GEANTInstanceCalo","ECAL");
    fGeantInstanceMuID = p.get<string>("GEANTInstanceMuID","MuID");

    fHitLabel         = p.get<string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<string>("TrackLabel","track");
    fVertexLabel      = p.get<string>("VertexLabel","vertex");
    fVeeLabel         = p.get<string>("VeeLabel","veefinder1");

    //readout simulation labels
    fTPCDigitLabel    = p.get<string>("TPCDigitLabel","daq");
    fCalDigitLabel    = p.get<string>("CalDigitLabel","daqsipm"); //same label for ECal/MuID
    fCalDigitInstance = p.get<string>("CalDigitInstance","ECAL");
    fMuDigitInstance  = p.get<string>("MuDigitInstance","MuID");

    fCaloHitLabel        = p.get<string>("CaloHitLabel","sipmhit");
    fCaloHitInstanceCalo = p.get<string>("CaloHitInstanceCalo","ECAL");
    fMuIDHitLabel        = p.get<string>("MuIDHitLabel","sipmhit");
    fCaloHitInstanceMuID = p.get<string>("CaloHitInstanceMuID","MuID");

    fClusterLabel     = p.get<string>("ClusterLabel","calocluster");
    fClusterInstanceCalo = p.get<string>("ClusterInstanceCalo","ECAL");
    fPFLabel          = p.get<string>("PFLabel","pandora");
    fECALAssnLabel    = p.get<string>("ECALAssnLabel","trkecalassn");

    // analysis configuration
    fAnaMode                  = p.get<string>("AnalysisMode",     "gen");
    fWriteDisplay             = p.get<bool>("WriteDisplay",       false);

    fIonizTruncate            = p.get<float>("IonizTruncate",    0.70);

    //required input products
    // gen
    consumesMany<vector<simb::MCTruth> >();
    consumesMany<vector<simb::GTruth> >();

    // g4
    consumes<Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
    consumes<vector<simb::MCParticle> >(fGeantLabel);
    consumes<vector<sdp::EnergyDeposit> >(fGeantLabel);
    InputTag ecalgeanttag(fGeantLabel, fGeantInstanceCalo);
    consumes<vector<sdp::CaloDeposit> >(ecalgeanttag); 

    // reco
    consumes<vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<vector<rec::Hit> >(fHitLabel);
    consumes<vector<rec::Track> >(fTrackLabel);
    consumes<vector<rec::Vertex> >(fVertexLabel);
    consumes<Assns<rec::Track, rec::Vertex> >(fVertexLabel);
    consumes<vector<rec::Vee> >(fVeeLabel);
    consumes<Assns<rec::Track, rec::Vee> >(fVeeLabel);
    consumes<vector<rec::TrackIoniz>>(fTrackLabel);
    consumes<Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

    // det
    if(fAnaMode == "readout"){

        //TPC
        InputTag tpcrawtag(fTPCDigitLabel);
        consumes<vector<raw::RawDigit> >(tpcrawtag);

        //ECal
        InputTag ecalrawtag(fCalDigitLabel, fCalDigitInstance);
        consumes<vector<raw::CaloRawDigit> >(ecalrawtag);

        //MuonID system 
        if (fGeo->HasMuonDetector()) {
            InputTag murawtag(fCalDigitLabel, fMuDigitInstance);
            consumes<vector<raw::CaloRawDigit> >(murawtag);
        }
    }

    InputTag ecalhittag(fCaloHitLabel, fCaloHitInstanceCalo);
    consumes<vector<rec::CaloHit> >(ecalhittag);

    //Muon system related
    if (fGeo->HasMuonDetector()) {
        InputTag muidgeanttag(fGeantLabel, fGeantInstanceMuID);
        consumes<vector<sdp::CaloDeposit> >(muidgeanttag);

        InputTag muidhittag(fCaloHitLabel, fCaloHitInstanceMuID);
        consumes<vector<rec::CaloHit> >(muidhittag);
    }

    InputTag clustertag(fClusterLabel, fClusterInstanceCalo);
    consumes<vector<rec::Cluster> >(clustertag); //fClusterLabel);
    consumes<Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

    return;
} // end constructor

//==============================================================================
//==============================================================================
//==============================================================================
void gar::StructuredTree::beginJob() {

    fTreeType = "structured";
    fTpcCenter.SetXYZT(fGeo->TPCXCent(),fGeo->TPCYCent(),fGeo->TPCZCent(),0.);
    xTPC = fGeo->TPCLength() / 2.;
    rTPC = fGeo->TPCRadius();
    fGeometry = fGeo->GDMLFile();

    art::ServiceHandle<art::TFileService> tfs;

    // headerTree
    fHeaderTree  = tfs->make<TTree>("headerTree", "sample provenance information");
    fHeaderTree->Branch("Run",           &fRun,            "Run/I");
    fHeaderTree->Branch("SubRun",        &fSubRun,         "SubRun/I");
    fHeaderTree->Branch("TreeType",      "std::string",    &fTreeType);
    fHeaderTree->Branch("POT",           &fPOT,            "POT/I");
    fHeaderTree->Branch("NSpills",       &fNSpills,        "NSpills/I");
    fHeaderTree->Branch("Geometry",      "std::string",    &fGeometry);
    fHeaderTree->Branch("TpcCenter",     "TLorentzVector", &fTpcCenter); 

    // genTree
    fGenTree     = tfs->make<TTree>("genTree",    "generator level info");
    fGenTree->Branch("Event",        &fEvent,       "Event/I");
    fGenTree->Branch("GIndex",       "std::vector<Int_t>",  &fGIndex);
    fGenTree->Branch("GTruth",       "std::vector<garana::GTruth>",  &fGTruth);
    fGenTree->Branch("FSParticles",  "std::vector<std::vector<garana::FSParticle>>", &fFSParticles);

    // g4Tree
    fG4Tree      = tfs->make<TTree>("g4Tree",     "GEANT4 level info");
    fG4Tree->Branch("Event",        &fEvent,       "Event/I");
    fG4Tree->Branch("TruthIndex",   "std::vector<UInt_t>", &fG4TruthIndex);
    fG4Tree->Branch("FSIndex",      "std::vector<UInt_t>", &fG4FSIndex);
    fG4Tree->Branch("G4Particles",  "std::vector<garana::G4Particle>",  &fG4Particles);
    if(fAnaMode == "reco" || fAnaMode == "readout"){
        fG4Tree->Branch("CalHits.", "std::vector<sdp::CaloDeposit>", &fSimCalHits);
        fG4Tree->Branch("TPCHits", "std::vector<sdp::EnergyDeposit>", &fSimTPCHits);
        //fG4Tree->Branch("CalHitsG4Indices.", "vector<UInt_t>", &fSimCalHitsG4P);
        //fG4Tree->Branch("TPCHitsG4Indices", "vector<UInt_t>", &fSimTPCHitsG4P);
        if(fGeo->HasMuonDetector()){
            fG4Tree->Branch("MuIDHits.", "std::vector<sdp::CaloDeposit>", &fSimMuHits);
            //fG4Tree->Branch("MuIDHitsG4Indices.", "vector<sdp::CaloDeposit>", &fSimMuHitsG4P);
        }
    }

    if (fWriteDisplay) {
        fDisplayTree = tfs->make<TTree>("displayTree","truth/reco level traj points for event display");
        fDisplayTree->Branch("Event",         &fEvent,                             "Event/I");
        fDisplayTree->Branch("MCTraj",        "std::vector<std::vector<TVector3>>", &fMCTraj);
        fDisplayTree->Branch("TrackTraj",     "std::vector<std::vector<TVector3>>", &fTrackTraj);
        fDisplayTree->Branch("G4Index",       "std::vector<UInt_t>",                &fDispG4Index); 
        fDisplayTree->Branch("TrackIndex",    "std::vector<UInt_t>",                &fDispTrkIndex);

        /*if(fAnaMode!="readout"){
            fDisplayTree->Branch("RecoTracks", "std::vector<garana::Track>",   &fRecoDisplay);
        }*/

    }
    //Digit. or Digit (no '.')??
    if(fAnaMode=="readout"){
        fDetTree     = tfs->make<TTree>("detTree",    "detector readout level info");
        fDetTree->Branch("Event",         &fEvent,       "Event/I");
        fDetTree->Branch("TPCDigits",   "std::vector<gar::raw::RawDigit>",     &fTPCDigits);
        fDetTree->Branch("CaloDigits.", "std::vector<gar::raw::CaloDigit>",    &fCaloDigits);
        if(fGeo->HasMuonDetector()){
            fDetTree->Branch("MuDigits.", "std::vector<gar::raw::CaloDigit>",    &fMuDigits);
        }

        return; //probably don't want reco info for readout simulation analysis unless doing hit reco ana
    }

    //recoTree
    fRecoTree    = tfs->make<TTree>("recoTree",   "reconstruction level info");
    fRecoTree->Branch("Event",            &fEvent,                      "Event/I");
    fRecoTree->Branch("Tracks",           "std::vector<garana::Track>",      &fTracks);
    fRecoTree->Branch("Vees",             "std::vector<garana::Vee>",      &fVees);
    fRecoTree->Branch("Vertices",         "std::vector<garana::Vertex>",   &fVertices);
    fRecoTree->Branch("CalClusters",      "std::vector<garana::CaloCluster>",  &fCalClusters);
    fRecoTree->Branch("TrackG4Indices",   "std::vector<std::vector<UInt_t>>",     &fTrackG4PIndices);
    fRecoTree->Branch("VertTrackIndices", "std::vector<std::vector<UInt_t>>",     &fVertTrackIndices);
    fRecoTree->Branch("VertTrackEnds",    "std::vector<std::vector<Int_t>>",      &fVertTrackEnds);
    fRecoTree->Branch("VeeTrackIndices",  "std::vector<std::vector<UInt_t>>",     &fVeeTrackIndices);
    fRecoTree->Branch("VeeTrackEnds",     "std::vector<std::vector<Int_t>>",      &fVeeTrackEnds);
    fRecoTree->Branch("CalTrackIndices",  "std::vector<std::vector<UInt_t>>",     &fCalTrackIndices);
    //fRecoTree->Branch("CalTrackEnds",   "std::vector<std::vector<Int_t>>", &fCalTrackEnds);
    fRecoTree->Branch("CalG4Indices",     "std::vector<std::vector<UInt_t>>",     &fCalG4PIndices);
    if(fGeo->HasMuonDetector()){
            fRecoTree->Branch("MuIDClusters", "std::vector<garana::CaloCluster>",   &fMuClusters);
    }

    //TODO: add Assns for reco r&d products
    if(fAnaMode == "reco"){
        fRecoTree->Branch("TPCHits",       "std::vector<gar::rec::Hit>",        &fTPCHits);
        fRecoTree->Branch("TPCClusters",   "std::vector<gar::rec::TPCCluster>", &fTPCClusters);
        fRecoTree->Branch("CalHits",       "std::vector<gar::rec::CaloHit>",    &fCalHits);
        if(fGeo->HasMuonDetector()){
            fRecoTree->Branch("MuIDHits",       "std::vector<gar::rec::CaloHit>", &fMuHits);
        }
    }


    //auxilliary info for PID calculation (not yet in reco)
    string filename = "${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
    TFile infile(filename.c_str(), "READ");

    m_pidinterp.clear();
    char str[11];
    for (int q = 0; q < 501; ++q)
    {
        sprintf(str, "%d", q);
        std::string s = "pidmatrix";
        s.append(str);
        // read the 500 histograms one by one; each histogram is a
        // 6 by 6 matrix of probabilities for a given momentum value
        m_pidinterp.insert( std::make_pair(q, (TH2F*) infile.Get(s.c_str())->Clone("pidinterp")) );
    }

    //map region names to id code for easy access from the trees
    // TODO this is for the current (as of 3/25/2021) default geometry
    //      -> add handling for arbitrary geometry
    fRegionNameToID = { 
                        {"Unknown",        -1 },
                        {"TPC1",           0  },
                        {"TPC2",           1  },
                        {"GArInactive",    2  },
                        {"BarrelECal",     3  },
                        {"EndcapECal",     4  },
                        {"Cryostat",       5  },
                        {"CryostatEndcap", 6  },
                        {"YokeBarrel",     7  }, 
                        {"YokeEndcap",     8  },
                        {"MPD",            9  }
                      };

    return;

}  // End of StructuredTree::beginJob

////==============================================================================
////==============================================================================
////==============================================================================
void gar::StructuredTree::endRun(art::Run const& run) {

   auto const& ID = run.id();

    art::Handle<sumdata::POTSummary> summaryHandle;
    if (!run.getByLabel(fPOTtag, summaryHandle)) {
        MF_LOG_DEBUG("anatree") << " No sumdata::POTSummary branch for run " << ID
        <<" Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        return;
    }

    sumdata::POTSummary const& RunPOT = *summaryHandle;
    fPOT = RunPOT.TotalPOT();
    fNSpills = RunPOT.TotalSpills();

    //fHeaderTree->Fill();

    MF_LOG_INFO("anatree") << "POT for this file is " << fPOT
    << " The number of spills is " << fNSpills;

}

void gar::StructuredTree::endJob() {

    fHeaderTree->Fill();
}

//==============================================================================
//==============================================================================
//==============================================================================
void gar::StructuredTree::analyze(art::Event const & e) {

    // Is the backtracker good to go?  Need to fill fBack on per-event basis
    // not in beginJob()
    cheat::BackTrackerCore const* const_bt = providerFrom<cheat::BackTracker>();
    fBt = const_cast<cheat::BackTrackerCore*>(const_bt);
    /*if ( !fBt->HasTracks() || !fBack->HasClusters() ) {
        throw cet::exception("MatchingPerformance") << " Not enough backtracker info"
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }*/

    ClearVectors();

    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();

    //Fill generator and MC Information
    mf::LogDebug("StructuredTree") << "fill truth level info";
    FillGenTree(e);
    FillG4Tree(e);
    fGenTree->Fill();
    fG4Tree->Fill();

    if(fWriteDisplay){
        mf::LogDebug("StructuredTree") << "fill display tree";
        fDisplayTree->Fill();
    }

    //Fill Readout Information
    if(fAnaMode=="readout") {
        mf::LogDebug("StructuredTree") << "fill det info";
        FillDetTree(e);
        fDetTree->Fill();
        mf::LogDebug("StructuredTree") << "done.";
        return; //no need for reco info (possible exception could be hit reco ana)
    }

    //Fill Reco Information
    if(fAnaMode=="reco"){
        mf::LogDebug("StructuredTree") << "fill base level reco info";
        FillRecoInfo(e);
    }

    //Fill HL Reco Information
    mf::LogDebug("StructuredTree") << "fill high level reco info";
    FillHighLevelRecoInfo(e);
    fRecoTree->Fill();

    mf::LogDebug("StructuredTree") << "done.";

    return;
}

//==============================================================================
//==============================================================================
//==============================================================================
void gar::StructuredTree::ClearVectors() {

    //Assns bookkeeing
    fMCTruths.clear();
    fMCParticles.clear();

    // genTree
    fGIndex.clear();
    fGTruth.clear();
    fFSParticles.clear();

    // g4Tree
    fG4TruthIndex.clear();
    fG4FSIndex.clear();
    fG4Particles.clear();

    // displayTree
    if (fWriteDisplay) {
        fMCTraj.clear();     
        fTrackTraj.clear();  
        fDispG4Index.clear();
        fDispTrkIndex.clear();
    }

    // detTree
    if(fAnaMode == "readout") {
        fTPCDigits.clear();
        fCaloDigits.clear();
        fMuDigits.clear();
    }

    // recoTree
    if(fAnaMode == "reco") {

        //g4 tree
        fSimTPCHits.clear();
        fSimCalHits.clear();
        //fSimTPCHitsG4P.clear(); //may not be needed. use vector of
        //fSimCalHitsG4P.clear(); //hits for each G4Particle (assuming each G4Part has one...)

        // Hit data
        fTPCHits.clear();
        fCalHits.clear();
        if(fGeo->HasMuonDetector()){
            fSimMuHits.clear(); //g4
            //fSimMuHitsG4P.clear();
            fMuHits.clear(); //reco
        }

        // TPC Cluster data
        fTPCClusters.clear();
    }

    //   Track data
    fTracks.clear();
    fTrackG4PIndices.clear();
    //TODO: for reco R&D, add Assns to TPC clusters

    // Vertex data
    fVertices.clear();
    fVertTrackIndices.clear();
    fVertTrackEnds.clear();

    // Vee data
    fVees.clear();
    fVeeTrackIndices.clear();
    fVeeTrackEnds.clear();

    //ECal/MuID clusters
    fCalClusters.clear();
    if(fGeo->HasMuonDetector()){
        fMuClusters.clear();
     }

    // ECAL cluster to track association info
    fCalTrackIndices.clear(); // index of fTracks associated with fCalClusters
    fCalTrackEnds.clear();   // Track end associated with the cluster
    fCalG4PIndices.clear();

    return;

} // end :StructuredTree::ClearVectors

//==============================================================================
void gar::StructuredTree::FillGenTree(art::Event const & e) {

    //Need to keep enough info for GENIE reweighting: GTruth + MCTruth

    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory

    auto mcthandlelist = e.getMany<vector<simb::MCTruth> >();
    auto gthandlelist  = e.getMany<vector<simb::GTruth> >();

    vector<string> genInstances; //instances for different generator products

    //all generators produce MCTruths, some produce associated GTruth
    for (size_t imcthl = 0; imcthl < mcthandlelist.size(); ++imcthl) {

        auto mcthandle = mcthandlelist.at(imcthl);
        string instance = mcthandle.provenance()->productInstanceName();
        if( std::find(genInstances.begin(),genInstances.end(),instance) == genInstances.end()) {
            mf::LogDebug("StructuredTree::FillGenTree")
                << "found new generator instance, '" << instance << "'";
            genInstances.push_back(instance);
        }
        else
            mf::LogWarning("StructuredTree::FillGenTree")
                << "found repeated generator instance, '" << instance << "'";

        mf::LogDebug("StructuredTree::FillGenTree")
            << "processing " << mcthandle->size() << " MCTruths";

        //loop over MCTruths
        for (size_t imct = 0; imct < mcthandle->size(); imct++) {

            auto const& mct = mcthandle->at(imct); //MCTruth
            fMCTruths.push_back(&mct);
            //if(fEvent==994)
		//std::cout << "neutrino in event 994 has CCNC = " << mct.GetNeutrino().CCNC() << std::endl;

            //check if MCTruth produced by GENIE
            //TO DO: we might not want to make all nu samples reweightable
            //       TPC - yes of course; 
            //       ECal - probably
            //       Dirt - probably not
            if(mct.GeneratorInfo().generator == simb::_ev_generator::kGENIE) {
                // for GTruth handles...
                //for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl) {
                for (auto const& gth : gthandlelist) {

                    auto const& gt = (*gth).at(0);
			
                    // only expect one GTruth per MCTruth
                    //if((*gthandlelist.at(igthl)).size()!=1) 
                    if((*gth).size()!=1)
                        throw cet::exception("StructuredTree") 
                                << "more than 1 GTruth in MCTruth! ("
                                //<< (*gthandlelist.at(igthl)).size() << ")";
                                << (*gth).size() << ")";

                    //adp::GarGTruth gt((*gthandlelist.at(igthl)).at(0));
                    // need to convert vertex units from meters to cm for GENIE data
                    garana::GTruth gtout = MakeAnaGTruth( gt, fRegionNameToID.at( PointToRegion(100*gt.fVertex.Vect()) ) );
                    gtout.fVertex -= fTpcCenter;
                    fGTruth.push_back(gtout); //MakeAnaGTruth( gt, fRegionNameToID.at( PointToRegion(100*gt.fVertex.Vect()) ) ) );
                    fGIndex.push_back(fGIndex.size());
                }
                mf::LogDebug("StructuredTree::FillGenTree")
                    << "filled vector with " << fGTruth.size() << " GTruths";
 
            }//if MCTruth from GENIE
            else{
                fGIndex.push_back(-1); //non-GENIE gen flag
            }
            
            fFSParticles.push_back({});//only one GTruth per MCTruth

            //loop over gen-level MCParticles, select FS particles only
            for(int imcp = 0; imcp<mct.NParticles(); imcp++) {

                auto const& mcp = mct.GetParticle(imcp);
                if(mcp.StatusCode() == 1){ //GENIE specific???

                    mf::LogDebug("StructuredTree::FillGenTree") 
                        << "FS particle with trackID " << (int)mcp.TrackId() << ", PDG "
                        << mcp.PdgCode() << ", " << "4-position (" << mcp.Position(0).X() << ", "
                        << mcp.Position(0).Y() << ", " << mcp.Position(0).Z() << ", "
                        << mcp.Position(0).T() << "), and "
                        << (int)mcp.NumberTrajectoryPoints() << " trajectory points";

                    fFSParticles.back().push_back(MakeFSParticle(mcp));

                }//if FS particle
            }//for MCParticles

            mf::LogDebug("StructuredTree::FillGenTree")
                << "Processed MCTruth with " << mct.NParticles() << " particles of which "
                << fFSParticles.back().size() << " are in the final state";
            
        }//for MCTruths
    }//for MCTruth handles

}//end FillGenTree()

//===============================================================================
void gar::StructuredTree::FillG4Tree( Event const & e) {

    // handle to MCParticles produced in G4 stage
    InputTag mcptag(fGeantLabel);
    auto MCPHandle = e.getValidHandle<vector<MCParticle>>(mcptag);

    //loop over MCParticles
    for ( auto mcp : *MCPHandle ) {

	//if(abs(mcp.PdgCode())==12 || abs(mcp.PdgCode())==14)
	//    std::cout << "neutrino MCParticle??? (PDG = " << mcp.PdgCode() 
        //              << " | trackID = " << mcp.TrackId() << ")" << std::endl;

        fMCParticles.push_back(&mcp);
	//std::cout << "MCParticle created in " << fGeo->VolumeName(mcp.Position(0).Vect()) << std::endl;

        int parentPdg = INT_MAX;
        int progenitorId = INT_MAX;
        int progenitorPdg = INT_MAX;
        vector<pair<TLorentzVector,TLorentzVector>> positions;
        vector<pair<TLorentzVector,TLorentzVector>> momenta; 
        vector<int> regions;
        vector<size_t> npointsPerRegion;

        if(mcp.Mother()>0) {
            parentPdg = (fBt->FindMother(&mcp))->PdgCode();
            auto const& progenitor = fBt->FindEve(&mcp);
            progenitorId = progenitor->TrackId();
            progenitorPdg = progenitor->PdgCode();
        }
        else {
            //std::cout 
            mf::LogDebug("StructuredTree") << "G4 primary with trackID " << mcp.TrackId() << ", PDG " << mcp.PdgCode() 
                      << ", initial 4-position (" << mcp.Position(0).X()
                      << "," << mcp.Position(0).Y() << ","
                      << mcp.Position(0).Z() << ", " << mcp.Position(0).T() << ")" ;
                      //<< std::endl;
        }
       

        // find regions of interest and get the entry and exit 4-vectors
        //std::cout << "looking for regions..." << std::endl;
        bool enter = false;
        size_t nptsreg=0;
        TLorentzVector *pEnter=0, *pExit=0, *xEnter=0, *xExit=0;
        for(size_t ipt = 0; ipt < mcp.NumberTrajectoryPoints(); ipt++){

            //if particle is not in a region of interest, skip to next traj point
            std::string regionName = PointToRegion(mcp.Position(ipt).Vect());  //fGeo->VolumeName(mcp.Position(ipt).Vect());
            if(fRegionNameToID.at(regionName)==-1)//==fRegionNameToID.end())
                continue;

            //std::cout << "trajectory point in " << regionName << " (" << 
            //             fRegionNameToID.at(regionName) << ")" << std::endl;;
            //is this the last trajectory point?
            bool last = (ipt==mcp.NumberTrajectoryPoints()-1);
            std::string nextName;
            if(!last) nextName = PointToRegion(mcp.Position(ipt+1).Vect());
            //if(nextName!=regionName) std::cout << "   next region is " << nextName << std::endl;

            // if the particle is entering into a region of interest
            if(!enter) {
                enter = true;
                nptsreg=0;
                nptsreg++;
                pEnter = new TLorentzVector(mcp.Momentum(ipt));
                xEnter = new TLorentzVector(mcp.Position(ipt));
                regions.push_back(fRegionNameToID.at(regionName));

                // is this the only point in this region?
                if(regionName != nextName) {

                    enter = false;
                    pExit = pEnter;
                    xExit = xEnter;
                    momenta.push_back({*pEnter,*pExit});
                    positions.push_back({*xEnter,*xExit});
                    npointsPerRegion.push_back(nptsreg);

                    /*if(momenta.size()==0 && !last)//first region and not last traj point
                        std::cout << "found MCParticle with single point in " << regionName << std::endl;
                    else if(momenta.size()==0 && last) //first region and not last traj point
                        std::cout << "found MCParticle with single and only point in " << regionName << std::endl;
                    else if(last) //not first region and last traj point
                        std::cout << "   last point in  " << regionName << " with a single point" << std::endl;
                    else //not first region and not last traj point
                        std::cout << "   now it's in  " << regionName << " with a single point" << std::endl; */
                }
            }//if now entering ROI

            if(enter)
                nptsreg++;

            //if last point or next region is not this one (this is the exit point)
            if(enter &&                         //if it entered, it can exit and ...
               (last || ( !last &&              //it's the last point or if not...
                regionName != nextName ) ) ) {  //the next region isn't this one

                /*if(momenta.size()==0 && !last)     //first region but not last traj point
                    std::cout << "found MCParticle in " << regionName << std::endl;
                else if(momenta.size()==0 && last) //first region and last traj point
                    std::cout << "found MCParticle with first and last point in " << regionName << std::endl;
                else if(last)//not first region but last traj point
                    std::cout << "   last point in " << regionName << std::endl;
                else //not first region and not last traj point
                    std::cout << "   now it's in  " << regionName << std::endl;

                std::cout << "number of traj. points in this region = " << nptsreg << std::endl;*/

                enter = false;
                pExit = new TLorentzVector(mcp.Momentum(ipt));
                xExit = new TLorentzVector(mcp.Position(ipt));

                momenta.push_back({*pEnter,*pExit});
                positions.push_back({*xEnter,*xExit});
                npointsPerRegion.push_back(nptsreg);

            }//if now exiting ROI

        }//for trajectory points
        //std::cout << "done." << std::endl;
        //if(momenta.size()!=0)
        //    std::cout << "  filled momenta/positions with " << momenta.size() << " pairs" << std::endl;

        // sanity check
        if(momenta.size()!=positions.size()) std::cout << "G4Particle: positions-momenta size mismatch" << std::endl;
        if(momenta.size()!=regions.size()) std::cout << "G4Particle: positions/momenta - regions size mismatch" << std::endl;

        // make a G4Particle
        fG4Particles.push_back(MakeG4Particle(mcp,parentPdg,progenitorPdg,progenitorId,positions,momenta,regions,npointsPerRegion));

        /*bool foundfs = false;
        for(UInt_t index=0; index<fFSParticles.size(); index++) {
            if(fFSParticles[index].TrackId() == progenitorId) {
                fG4FSIndex.push_back(index);
                foundfs = true;
                break;
            }
        }
        if(!foundfs) {
           std::cout  <<
           "WARNING in StructuredTree_module::FillG4Tree: FSParticle association not found!" 
           << std::endl;
           fG4FSIndex.push_back(UINT_MAX); //need to push back something or association ordering is broken
        }*/

        //find associated MCTruth
        auto mcmatch = fBt->ParticleToMCTruth(&mcp);
        bool found = false;
        for(size_t i=0; i<fMCTruths.size(); i++){
            if(fMCTruths[i]==mcmatch.get()) { //compare pointers to check MCTruth equality
                found = true;
                fG4TruthIndex.push_back(i);

                bool foundfs = false;
                for(UInt_t index=0; index<fFSParticles.size(); index++) {
                    if(fFSParticles[i][index].TrackId() == progenitorId) {
                        fG4FSIndex.push_back(index);
                        foundfs = true;
                        break;
                    }
                }
                if(!foundfs) {
                   //std::cout  <<
                   //"WARNING in StructuredTree_module::FillG4Tree: FSParticle association not found!"
                   //<< std::endl;
                   fG4FSIndex.push_back(UINT_MAX); //need to push back something or association ordering is broken
                }
        
                break;
            } 
        }
        if(!found){
            mf::LogWarning("StructuredTree::FillG4Tree") 
                << "no MCTruth found for MCParticle";
            fG4TruthIndex.push_back(UINT32_MAX);
        }

        if(fWriteDisplay){
            const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
            const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
            //No charge don't store the trajectory
            if (definition==nullptr || definition->Charge() == 0) continue;

            //std::cout << "number of trajectory points pre sparsify: " << mcp.NumberTrajectoryPoints() << std::endl;
	    mcp.SparsifyTrajectory();
            //std::cout << "number of trajectory points post sparsify: " << mcp.NumberTrajectoryPoints() << std::endl;
            fDispG4Index.push_back(fMCParticles.size()-1);
            fMCTraj.emplace_back(vector<TVector3>(mcp.NumberTrajectoryPoints(),TVector3()));
            for(size_t i=0; i<mcp.NumberTrajectoryPoints(); i++){
                fMCTraj.back()[i] = mcp.Position(i).Vect();
            }
            //std::cout << "done filling truth trajectories" << std::endl;
        }

    }//for MCParticles

    // Get handles for SimHits used for reco R&D
    if(fAnaMode == "reco") {

        InputTag geanttag(fGeantLabel); //TPC
        InputTag calgeanttag(fGeantLabel, fGeantInstanceCalo); //ECal
        auto TPCSimHitHandle = e.getValidHandle< vector<sdp::EnergyDeposit> >(geanttag); //TPC
        auto CalSimHitHandle = e.getValidHandle< vector<sdp::CaloDeposit> >(calgeanttag); //ECal
        // Note: no backtracker for sim hit(s) <-> MCParticle so we do it the hard way...
        //  

        fSimTPCHits.resize(fG4Particles.size()); //want MCParticle to have same index as its hits
        fSimCalHits.resize(fG4Particles.size());
        if(fGeo->HasMuonDetector())
            fSimMuHits.resize(fG4Particles.size());

        //try to find SimHits with matching MCParticle
        // there can (and usually will) be multiple hits/particle
        size_t nnomatch = 0;
        for(size_t ig4p=0; ig4p<fG4Particles.size(); ig4p++){

            bool hasmatched = false;
            fSimTPCHits.push_back({});
            fSimCalHits.push_back({});

            // TPC
            for( auto const& hit : *TPCSimHitHandle ){
                if(hit.TrackID() == fG4Particles[ig4p].TrackID()){
                    fSimTPCHits[ig4p].push_back(hit);
                    hasmatched = true;
                }
            }//for TPC hits

            // Save simulation ecal hit info
            for ( auto const& hit : *CalSimHitHandle ) {
                if(hit.TrackID() == fG4Particles[ig4p].TrackID()){
                    fSimCalHits[ig4p].push_back(hit);
                    hasmatched = true;
                }
            }//for ECal hits

            // Save simulation muon system hit info
            if(fGeo->HasMuonDetector()) {

                InputTag mugeanttag(fGeantLabel, fGeantInstanceMuID); //MuID
                auto MuSimHitHandle = e.getValidHandle< vector<sdp::CaloDeposit> >(mugeanttag); //MuID
                fSimMuHits.push_back({});

                for ( auto const& hit : *MuSimHitHandle ) {
                    if(hit.TrackID() == fG4Particles[ig4p].TrackID()) {
                        fSimMuHits[ig4p].push_back(hit);
                        hasmatched = true;
                    }
                }//for MuID hits

            }//if MuID
            if(!hasmatched)
                nnomatch++;

        }//for G4Particles


        if(nnomatch>0){
            mf::LogWarning("StructuredTree::FillG4Tree")
                << " found " << nnomatch << " G4Particle(s) with no matching SimHit."
                << " SimHit -> G4Particle associations corrupted!";
        }

    }//if mode reco

}//FillG4Tree

//==============================================================================

void gar::StructuredTree::FillDetTree( Event const & e) {

        // Get handles for RawDigits
        InputTag tpcdigittag(fTPCDigitLabel);
        InputTag caldigittag(fCalDigitLabel, fCalDigitInstance);

        // request valid handle -> throws exception if product not found
        auto TPCDigitHandle = e.getValidHandle< vector<raw::RawDigit> >(tpcdigittag);
        auto CalDigitHandle = e.getValidHandle< vector<raw::CaloRawDigit> >(caldigittag);

        //TPC digits
        for( auto const& digit : *TPCDigitHandle)
            fTPCDigits.push_back(digit);

        // save ecal raw digits info
        for ( auto const& digit : *CalDigitHandle ) {
            fCaloDigits.push_back(digit);
        }

        // save muon system raw digits info if MuID present
        if(fGeo->HasMuonDetector()) {

            InputTag mudigittag(fCalDigitLabel, fMuDigitInstance);
            auto MuDigitHandle = e.getValidHandle< vector<raw::CaloRawDigit> >(mudigittag);

            for ( auto const& digit : *MuDigitHandle ) {
                fMuDigits.push_back(digit);
            }
        }
}//FillDetTree

//==============================================================================
// base level reconstruction products (used for reco R&D mode)
void gar::StructuredTree::FillRecoInfo( Event const & e) {


    // reco hits
    // Get handles
    InputTag tpchittag(fHitLabel);
    InputTag calhittag(fCaloHitLabel, fCaloHitInstanceCalo);
    auto TPCHitHandle = e.getValidHandle< vector<rec::Hit> >(tpchittag);
    auto CalHitHandle = e.getValidHandle< vector<rec::CaloHit> >(calhittag);

    // save reco hits from the TPC
    for ( auto const& hit : *TPCHitHandle ) {
        fTPCHits.push_back(hit);
    }


    // save reco hits from the ECal
    for ( auto const& hit : *CalHitHandle ) {
        fCalHits.push_back(hit);
    }
 
    // save reco hits from the MuID (if present)
    if(fGeo->HasMuonDetector()) {

        InputTag muhittag(fMuIDHitLabel, fCaloHitInstanceMuID);
        auto MuHitHandle = e.getValidHandle< vector<rec::CaloHit> >(muhittag);

        for ( auto const& hit : *MuHitHandle ) {
            fMuHits.push_back(hit);
        }
    }

    // reco clusters
    // Get handles
    InputTag tpcclustertag(fTPCClusterLabel);
    InputTag tracktag(fTrackLabel);
    auto TPCClusterHandle = e.getValidHandle< vector<rec::TPCCluster> >(tpcclustertag);
    //auto TrackHandle      = e.getValidHandle< vector<rec::Track> >(tracktag);
    //art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
    //findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);

    //FIX ME: ASSOCIATIONS
    // save clusters in the TPC. For some reason, can't get FindOneP<rec::Track> or
    // FindManyP<rec::Track> to work; seems the underlying Assn isn't found.  Have
    // to FindManyP<TPCCluster> instead and  iterate if (fWriteTracks).  :(
    for ( auto const& TPCCluster : (*TPCClusterHandle) ) {

        fTPCClusters.push_back(TPCCluster);
        //TODO: fix assns
        /*Int_t trackForThisTPCluster = -1;
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
        pushit:
            fTPCClusterTrkIDNumber.push_back(trackForThisTPCluster);*/
    }//for TPCClusters

}//FillReco

//==============================================================================

void gar::StructuredTree::FillHighLevelRecoInfo( Event const & e) {


    // Get handles for Track, TrackIoniz, Vertex, Vee and also Assn's
    InputTag tracktag(fTrackLabel);
    InputTag iontag(fTrackLabel);
    InputTag verttag(fVertexLabel);
    InputTag veetag(fVeeLabel);
    InputTag calclustertag(fClusterLabel,fClusterInstanceCalo); //instance too?
    InputTag caltrackassntag(fECALAssnLabel,fClusterInstanceCalo);
    auto TrackHandle      = e.getValidHandle< vector<rec::Track> >(tracktag);
    //auto TrackIonHandle   = e.getValidHandle< vector<rec::TrackIoniz> >(iontag);
    auto VertexHandle     = e.getValidHandle< vector<rec::Vertex> >(verttag);
    auto VeeHandle        = e.getValidHandle< vector<rec::Vee> >(veetag);
    auto CalClusterHandle = e.getValidHandle< vector<rec::Cluster> >(calclustertag);

    //associations (FIX ME)
    art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
    findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,e,fTrackLabel);

    // Vertices Assn's to Tracks
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
    findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);

    // Vees Assn's to Tracks
    art::FindManyP<rec::Track, rec::TrackEnd>* findManyVeeTrackEnd = NULL;
    findManyVeeTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VeeHandle,e,fVeeLabel);

    // ECal cluster Assn's to Tracks
    // TODO: add track end assns when available (add on to BackTrackerCore)
    //art::FindManyP<rec::Track, rec::TrackEnd>* findManyCalTrackEnd = NULL;
    art::FindManyP<rec::Track, int>* findManyCalTrackEnd = NULL;
    //findManyCalTrackEnd = new art::FindManyP<rec::Track,rec::TrackEnd>(CalClusterHandle,e,calclustertag); //fECALAssnLabel);
    findManyCalTrackEnd = new art::FindManyP<rec::Track,int>(CalClusterHandle,e,caltrackassntag);
    //art::FindManyP<MCParticle>* findManyCalMCP = NULL;
    //findManyCalMCP = new art::FindManyP<MCParticle>(CalClusterHandle,e,fECALAssnLabel);
    if(findManyCalTrackEnd==nullptr) std::cout << "WARNING: ECAL cluster to track assns return nullptr!" << std::endl;
    if(findManyCalTrackEnd->size()==0) std::cout << "WARNING: ECal cluster to track assns returns empty vector!" << std::endl;

    // save per-track info
    size_t iTrack = 0;
    vector<size_t> trackIDs;
    for ( auto /*const&*/ track : *TrackHandle ) { //backtracker takes ptr to non-const obj

        //get track -> G4Particle Assns
        vector<pair<MCParticle*,float>> trackmcps = fBt->TrackToMCParticles(&track); //MCParticle w/fraction energy contributed to track
        //std::cout << "size of track <-> G4p vec from BT: " << trackmcps.size() << std::endl;
        vector<UInt_t> indices;
        for(auto const& mcpe : trackmcps){
            //std::cout << "looking for G4P match to associated match MCP with track ID, " << mcpe.first->TrackId() << std::endl;
            for(size_t i=0; i< fG4Particles.size(); i++) {
                //std::cout << "  check G4Particle " << i << " with track ID, " << fG4Particles[i].TrackID() << std::endl;
                if( mcpe.first->TrackId() == fG4Particles[i].TrackID()) //fMCParticles and fG4Particles have same mapping
                    indices.push_back(i);
            }
        }

        // std::cout << "found " << indices.size() << " track <-> G4P matches" << std::endl;
        fTrackG4PIndices.push_back(indices);

        //Reconstructed momentum forward and backward
        vector< pair<int, float> > pidF = processPIDInfo( track.Momentum_beg() );
        vector< pair<int, float> > pidB = processPIDInfo( track.Momentum_end() );
 
        //average ionization
        float avgIonF = 0., avgIonB=0.;
        if (findIonization->isValid()) {
            // No calibration for now.  Someday this should all be in reco
            rec::TrackIoniz ionization = *(findIonization->at(iTrack));
            processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
        } 

        // track truth-matching info
        TVector3 trkVtx(track.Vertex()[0],track.Vertex()[1],track.Vertex()[2]);
        vector<art::Ptr<rec::Hit>> hits = fBt->TrackToHits(&track); // not time ordered
        std::sort(hits.begin(),hits.end(),
          [trkVtx] ( const art::Ptr<rec::Hit>& hl, const art::Ptr<rec::Hit>& hr ) -> bool {
                 TVector3 rhitl(hl->Position()[0],hl->Position()[1],hl->Position()[2]);
                 TVector3 rhitr(hr->Position()[0],hr->Position()[1],hr->Position()[2]);
                 return (trkVtx-rhitl).Mag() < (trkVtx-rhitr).Mag(); 
          });

	// first and last rec::Hit
        auto hitIdesBeg = fBt->HitToHitIDEs(*(hits.begin()));
        auto hitIdesEnd = fBt->HitToHitIDEs(*(hits.end()-1));
        vector<pair<UInt_t,TLorentzVector>> truePosBeg, trueMomBeg, truePosEnd, trueMomEnd;
        vector<std::pair<int,float>> trueEnergy;

        for(auto const& ide : hitIdesBeg) {
            truePosBeg.push_back(std::make_pair(ide.trackID , ide.position));
            trueMomBeg.push_back(std::make_pair(ide.trackID , ide.momentum));
        }

        for(auto const& ide : hitIdesEnd) {
            truePosEnd.push_back(std::make_pair(ide.trackID , ide.position));
            trueMomEnd.push_back(std::make_pair(ide.trackID , ide.momentum));
        }

        std::vector<std::pair<simb::MCParticle*,float>> trkmcps = fBt->TrackToMCParticles(&track);
	double etot = fBt->TrackToTotalEnergy(&track);
	for(auto const& mcpfrac : trkmcps) {
            trueEnergy.push_back(std::make_pair(mcpfrac.first->TrackId(),etot*mcpfrac.second));
        }

        // make a garana::Track
        trackIDs.push_back(track.getIDNumber());
        fTracks.push_back(MakeAnaTrack(track, pidF, pidB, avgIonF, avgIonB,truePosBeg,truePosEnd,trueMomBeg,
                                       trueMomEnd, trueEnergy) );

        if(fWriteDisplay){
            art::FindOneP<rec::TrackTrajectory>*  findTrajectory = NULL;
            findTrajectory  = new art::FindOneP<rec::TrackTrajectory>(TrackHandle,e,fTrackLabel);

            const vector<TVector3> traj = findTrajectory->at(iTrack)->getFWDTrajectory();
            fTrackTraj.emplace_back(vector<TVector3>(traj.size(),TVector3()));
            fDispTrkIndex.push_back(iTrack);

	    for(size_t itrj = 0; itrj < fTrackTraj.back().size(); itrj++) {
                fTrackTraj.back()[itrj] = traj.at(itrj);
            }
        }

        iTrack++;
    } //for tracks

    // save Vertex and Track-Vertex association info
    size_t iVertex = 0;
    for ( auto const& vertex : *VertexHandle ) {

        fVertices.push_back(MakeAnaVtx(vertex));
        fVertTrackIndices.push_back({});
        fVertTrackEnds.push_back({});


        int nVertexedTracks = 0;
        if ( findManyTrackEnd->isValid() ) {
            nVertexedTracks = findManyTrackEnd->at(iVertex).size();
        }
        //fVertexN.push_back(nVertexedTracks); //number of tracks assced with this vertex

        //int vertexCharge = 0;
        for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack) {

            // Get this vertexed track.
            rec::Track track = *(findManyTrackEnd->at(iVertex).at(iVertexedTrack));
            rec::TrackEnd fee = *(findManyTrackEnd->data(iVertex).at(iVertexedTrack));

            // find its index in the track list
            for(size_t itrk=0; itrk<fTracks.size(); itrk++){
                if(trackIDs[itrk] == track.getIDNumber()){
                    fVertTrackEnds.back().push_back(fee);
                    fVertTrackIndices.back().push_back(itrk);
                    break;
                }
            }

            // add Vertex charge?
            /*if (fee==rec::TrackEndBeg) {
                vertexCharge += track.ChargeBeg();
            } else {
                vertexCharge += track.ChargeEnd();
            }*/
        }

        //fVertexQ.push_back(vertexCharge);
        ++iVertex;
    } // end loop over VertexHandle

    // save Vee and Track-Vee association info
    size_t iVee = 0;
    for ( auto const& vee : *VeeHandle ) {

        fVees.push_back(MakeAnaVee(vee));
        fVeeTrackIndices.push_back({});
        fVeeTrackEnds.push_back({});

        int nVeeTracks = 0;
        if ( findManyVeeTrackEnd->isValid() ) {
            nVeeTracks = findManyVeeTrackEnd->at(iVee).size();
        }

        for (int iVeeTrack=0; iVeeTrack<nVeeTracks; ++iVeeTrack) {

            // Get this vertexed track.
            rec::Track track = *(findManyVeeTrackEnd->at(iVee).at(iVeeTrack));
            rec::TrackEnd fee = *(findManyVeeTrackEnd->data(iVee).at(iVeeTrack));

            // find its index in the track list
            for(size_t itrk=0; itrk<fTracks.size(); itrk++){
                if(trackIDs[itrk] == track.getIDNumber()){
                    fVeeTrackIndices.back().push_back(itrk);
                    fVeeTrackEnds.back().push_back(fee);
                    break;
                }
            }
        }

        ++iVee;
    } // for Vees

    // Get info for ECal clusters and for ECAL-matched tracks
    size_t iCluster = 0;
    size_t n0edep = 0; // number of clusters with no associated true energy deposits
    for ( auto cluster : *CalClusterHandle ) {

	size_t n0ide=0; // number of hits with no associated CalIDEs
        vector<std::pair<int,float>> edeps; // trackID , true deposited energy [GeV]
	std::string regionName = PointToRegion(TVector3(cluster.Position()));  //fGeo->VolumeName(mcp.Position(ipt).Vect());
        int clustreg = fRegionNameToID.at(regionName);
        const vector<art::Ptr<rec::CaloHit>> hits = fBt->ClusterToCaloHits(&cluster);
        //if(hits.size()==0) std::cout << "cluster " << iCluster << ": hits is empty!" << std::endl;
        //else std::cout << "cluster " << iCluster << " with " << hits.size() << " hits" << std::endl;
        for(auto const& hit : hits) {

            auto const& ides = fBt->CaloHitToCalIDEs(hit);
            if(ides.size()==0) n0ide++; //std::cout << "have a hit, but ides is empty!" << std::endl;
            for(auto const& ide : ides) {

                std::pair<int,float> trkdep = std::make_pair(ide.trackID,ide.energyTot);
                if(std::find(edeps.begin(),edeps.end(),trkdep) == edeps.end()) {
                    edeps.push_back(trkdep);
                }
                else {
                    (*(std::find(edeps.begin(),edeps.end(),trkdep))).second+=trkdep.second;
                }

            }// for calo IDEs

        }// for calo hits

	if(hits.size()!=0 && n0ide>0) mf::LogDebug("FillHighLevelReco") 
                                                << "     hits with no IDEs: " << n0ide << " of " << hits.size() 
                                                << " hits affected (" << 100.0*(hits.size()-n0ide)/hits.size() 
                                                << " %)"; // << std::endl;

	if(edeps.size()==0) {
            mf::LogDebug("FillHighLevelReco") << "have a cluster, but edeps is empty!";
            edeps.push_back({-1,-1});
            n0edep++;
        }
	//else std::cout << "  found " << edeps.size() << "  truth edeps" << std::endl;

        fCalClusters.push_back(MakeAnaCalCluster(cluster,clustreg,edeps));
        fCalTrackIndices.push_back({});
        fCalTrackEnds.push_back({});
        fCalG4PIndices.push_back({});

        size_t nMatch = 0;
        if ( findManyCalTrackEnd->isValid() ) {
            nMatch = findManyCalTrackEnd->at(iCluster).size();
        }
        else{
            std::cout << "WARNING: findManyCalTrackEnd is invalid" << std::endl;
        }
        for (size_t itrk=0; itrk<nMatch; itrk++) {
            rec::Track track  = *(findManyCalTrackEnd->at(iCluster).at(itrk));
            //TODO: add track ends when they become available
            //rec::TrackEnd fee = *(findManyCalTrackEnd->data(iCluster).at(itrk));

            for(size_t imatch=0; imatch<fTracks.size(); imatch++){
                if(trackIDs[imatch] == track.getIDNumber()){
                    fCalTrackIndices.back().push_back(imatch);
                    //fCalTrackEnds.back().push_back(fee);
                    break;
                }
            }

        }//for matched tracks
        //std::cout << "filled ECal cluster -> track associations. on to MCParticles..." << std::endl;

        // list of MCParticles associated with this cluster and the fraction of the total energy contributed
	vector<std::pair<MCParticle*,float>> clustToMCPs = fBt->ClusterToMCParticles(&cluster);

	//need to find MCParticle index
	for(size_t imcp=0; imcp<clustToMCPs.size(); imcp++) {

	    vector<int> pdgs;
	    for(size_t index=0; index<fG4Particles.size(); index++){        

                if(fG4Particles.at(index).TrackID() == clustToMCPs.at(imcp).first->TrackId()){
		    pdgs.push_back(fG4Particles.at(index).PDG());
    	            fCalG4PIndices.back().push_back(index);
                    break;
                }// if cluster and MCParticle track IDs match
            }// for G4Particles
        }// for cluster-MCParticle matches
	//std::cout << "filled ECal cluster -> MCParticle associations" << std::endl;

        /*const vector<art::Ptr<rec::CaloHit>> hits = fBt->ClusterToCaloHits(&cluster);
        for(auto const& hit : hits) {

            const vector<sdp::CalIDE> ides = CaloHitToCalIDEs(hit);
            for(auto const& ide : ides) {

                

            }

        }// for calo hits*/


        iCluster++;
    }// for ECal clusters
    //std::cout << "number ECal clusters: " << iCluster << std::endl;
    //std::cout << "cluster truth matching 'efficiency': " << 1.0*(iCluster-n0edep)/iCluster << std::endl;


    return;

} // end FillHighLevelReco()



//==============================================================================
//==============================================================================
// Helper Functions
//==============================================================================

// determine if point could be in a ROI
std::string gar::StructuredTree::PointToRegion(const TVector3& point) {

    std::string region = "Unknown";
    std::string volName = fGeo->VolumeName(point);

    //only care about volumes in MPD
    if( fGeo->PointInMPD(point) ) {

        region = "MPD";

        //GArTPC->TPCChamber->TPCGas->TPCDrift
        if( fGeo->PointInGArTPC(point) ) {
            region = "TPC";
            if(volName.find("1")!=std::string::npos)
                region += "1";
            else if(volName.find("2")!=std::string::npos)
                region += "2";
            else
                region = "GArInactive";
        }
        else if ( fGeo->PointInECALBarrel(point) ) region = "BarrelECal";
        else if ( fGeo->PointInECALEndcap(point) ) region = "EndcapECal";

        //FIXME: probably missing the daughter volumes below
        else if ( volName == "volCryostat" )
            region = "Cryostat";

        else if ( volName == "volEndcapCryostat" )
            region = "EndcapCryostat";

        else if ( volName == "volYokeBarrel" )
            region = "YokeBarrel";

        else if ( volName == "volYokeEndcap" )
            region = "YokeEndcap";
    }

    return region;
}

// Process ionization.  Eventually this moves into the reco code.
void gar::StructuredTree::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
           float& forwardIonVal, float& backwardIonVal) {

    // NO CALIBRATION SERVICE FOR NOW

    std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
    forwardIonVal = processOneDirection(SigData, ionizeTruncate);

    SigData = ion.getBAK_dSigdXs();
    backwardIonVal = processOneDirection(SigData, ionizeTruncate);

    return;
}



float gar::StructuredTree::processOneDirection(std::vector<std::pair<float,float>> SigData, 
               float ionizeTruncate) {

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

//==============================================================================
//==============================================================================
//==============================================================================
vector< pair<int, float> > gar::StructuredTree::processPIDInfo( float p )
{
    vector<string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};
    vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
    vector< pair<int, float> > pid;
    pid.resize(6);

    int qclosest = 0;
    float dist = 100000000.;
    CLHEP::RandFlat FlatRand(fEngine);

    for (int q = 0; q < 501; ++q)
    {
        //Check the title and the reco momentum take only the one that fits
        string fulltitle = m_pidinterp[q]->GetTitle();
        unsigned first = fulltitle.find("=");
        unsigned last = fulltitle.find("GeV");
        string substr = fulltitle.substr(first+1, last - first-1);
        float pidinterp_mom = std::atof(substr.c_str());
        //calculate the distance between the bin and mom, store the q the closest
        float disttemp = std::abs(pidinterp_mom - p);

        if( disttemp < dist ){
            dist = disttemp;
            qclosest = q;
        }
    } // closes the "pidmatrix" loop

    //////Compute all the probabities for each type of true to reco
    //////loop over the columns (true pid)
    ////for (int pidm = 0; pidm < 6; ++pidm)
    ////    //calculate the distance between the bin and mom, store the q the closest
    ////    float disttemp = std::abs(pidinterp_mom - p);

    ////    if( disttemp < dist ){
    ////        dist = disttemp;
    ////        qclosest = q;
    ////    }
    ////} // closes the "pidmatrix" loop

    //Compute all the probabities for each type of true to reco
    //loop over the columns (true pid)
    for (int pidm = 0; pidm < 6; ++pidm)
    {
        //loop over the columns (true pid)
        vector< pair<float, string> > v_prob;

        //loop over the rows (reco pid)
        for (int pidr = 0; pidr < 6; ++pidr)
        {
            string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
            float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1, pidr+1);
            //Need to check random number value and prob value then associate the recopdg to the reco prob
            v_prob.push_back( std::make_pair(prob, recoparticlename) );
        }

        //Compute the pid from it
        if(v_prob.size() > 1)
        {
            //Order the vector of prob
            std::sort(v_prob.begin(), v_prob.end());
            //Throw a random number between 0 and 1
            float random_number = FlatRand.fire();
            //Make cumulative sum to get the range
            std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});

            for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++)
            {
                if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first )
                {
                    pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ), v_prob.at(ivec+1).first) );
                }
            }//for probs
        }//if found
        else
        {
            pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ), v_prob.at(0).first) );
        }
    }//compute probabilities

    //return a vector of pid and prob
    return pid;

}//processPIDInfo()

/*string gar::StructuredTree::GetGenCodeString(const simb::MCTruth& mct){

    string gen = "";
    auto code = mct.GeneratorInfo().generator;

    if (code == simb::_ev_generator::kGENIE) gen = "genie";
    else if (code == simb::_ev_generator::kCRY) gen = "cry";
    else if (code == simb::_ev_generator::kCRY) gen = "single";
    else {
        gen = "other";
        std::cout << "Warning: unrecognized generator code encountered!" << std::endl;
    }

    return gen;
}*/

DEFINE_ART_MODULE(gar::StructuredTree)
