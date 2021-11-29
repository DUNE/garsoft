////////////////////////////////////////////////////////////////////////////////
// Class:       anatest
// Plugin Type: analyzer (art v2_11_02)
// File:        anatest_module.cc
//
// Generated at Tue Mar 3 11:11:11 2020 by Leo Bellantoni
//
// a plugin for testing reco output - derived from anatree_module.cc
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

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"
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
#include <unordered_map>

#include <chrono>
using namespace std::chrono;



namespace gar {

    class anatest : public art::EDAnalyzer {
    public:
        explicit anatest(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        anatest(anatest const &) = delete;
        anatest(anatest &&) = delete;
        anatest & operator = (anatest const &) = delete;
        anatest & operator = (anatest &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;



    private:

        void ClearVectors();
        void FillVectors(art::Event const & e);

        // Position of TPC from geometry service; 1 S Boston Ave.
        double ItsInTulsa[3];

        // Input data labels
        std::vector<std::string> fGeneratorLabels;
        std::string fGeantLabel;
        std::string fHitLabel;
        std::string fTPCClusterLabel;
        std::string fTrackLabel;
        std::string fVertexLabel;
        std::string fRawCaloHitLabel;
        std::string fCaloHitLabel;
        std::string fClusterLabel;
        std::string fECALAssnLabel;

        // the analysis output tree
        TTree *fTree;

        // global event info
        Int_t fEvent;        ///< number of the event being processed
        Int_t fRun;          ///< number of the run being processed
        Int_t fSubRun;       ///< number of the sub-run being processed

        // Some random fields for now
        std::vector<Int_t>              fNeutrinoType;
        std::vector<Int_t>              fInteractionType;
        std::vector<Float_t>            fQ2;
    };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatest::anatest(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

    bool usegenlabels =
        p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
    if (!usegenlabels) fGeneratorLabels.clear();

    fGeantLabel       = p.get<std::string>("GEANTLabel","geant");
    fHitLabel         = p.get<std::string>("HitLabel","hit");
    fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
    fTrackLabel       = p.get<std::string>("TrackLabel","track");
    fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
    fRawCaloHitLabel  = p.get<std::string>("RawCaloHitLabel","daqecal");
    fCaloHitLabel     = p.get<std::string>("CaloHitLabel","calohit");
    fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");
    fECALAssnLabel    = p.get<std::string>("ECALAssnLabel","trkecalassn");

    fTree = nullptr;

    if (usegenlabels) {
        for (size_t i=0; i<fGeneratorLabels.size(); ++i) {
            consumes<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
        }
    } else {
        consumesMany<std::vector<simb::MCTruth> >();
    }



    consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);

    consumes<std::vector<rec::Hit> >(fHitLabel);
    consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Track> >(fTrackLabel);
    consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
    consumes<std::vector<rec::Vertex> >(fVertexLabel);
    consumes<art::Assns<rec::Track, rec::Vertex> >(fVertexLabel);

    consumes<std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
    consumes<std::vector<raw::CaloRawDigit> >(fRawCaloHitLabel);
    consumes<std::vector<rec::CaloHit> >(fCaloHitLabel);
    consumes<std::vector<rec::Cluster> >(fClusterLabel);
    consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

    return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::beginJob() {

    art::ServiceHandle<geo::GeometryGAr> euclid;
    ItsInTulsa[0] = euclid->TPCXCent();
    ItsInTulsa[1] = euclid->TPCYCent();
    ItsInTulsa[2] = euclid->TPCZCent();



    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("MPDtestTree","MPDtestTree");



    fTree->Branch("Run",           &fRun,         "Run/I");
    fTree->Branch("SubRun",        &fSubRun,      "SubRun/I");
    fTree->Branch("Event",         &fEvent,       "Event/I");

    fTree->Branch("NType",         &fNeutrinoType);
    fTree->Branch("InterT",        &fInteractionType);
    fTree->Branch("MC_Q2",         &fQ2);

    return;
}  // End of :anatest::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::analyze(art::Event const & e) {

    ClearVectors();
    FillVectors(e);
    fTree->Fill();
    return;
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::ClearVectors() {

    // clear out all our vectors
    fNeutrinoType.clear();
    fInteractionType.clear();
    fQ2.clear();

    return;
} // end :anatest::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatest::FillVectors(art::Event const & e) {



    // =============  Get art handles ==========================================
    // Get handles for MCinfo, also good for MCPTrajectory
    std::vector< art::Handle<std::vector<simb::MCTruth>> > mcthandlelist;
    if (fGeneratorLabels.size()<1) {
        mcthandlelist = e.getMany<std::vector<simb::MCTruth> >();    // get them all (even if there are none)
    } else {
        mcthandlelist.resize(fGeneratorLabels.size());
        for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
	  mcthandlelist.at(i) = e.getHandle<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
            if (!mcthandlelist.at(i))
                throw cet::exception("anatest") << " No simb::MCTruth branch."
                     << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
    }



    auto MCPHandle = e.getHandle<std::vector<simb::MCParticle> >(fGeantLabel);
    if (!MCPHandle) {
        throw cet::exception("anatest") << " No simb::MCParticle branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handles for MCCaloInfo
    auto SimHitHandle = e.getHandle< std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
    if (!SimHitHandle) {
        throw cet::exception("anatest") << " No gar::sdp::CaloDeposit branch."
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for TPC hit data
    auto HitHandle = e.getHandle< std::vector<rec::Hit> >(fHitLabel);
    if (!HitHandle) {
        throw cet::exception("anatest") << " No rec::Hit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handle for TPCClusters
    auto TPCClusterHandle = e.getHandle< std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    if (!TPCClusterHandle) {
        throw cet::exception("anatest") << " No rec::TPCCluster branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }



    // Get handles for Tracks and also Assn's to TPCClusters, TrackIoniz
    auto TrackHandle = e.getHandle< std::vector<rec::Track> >(fTrackLabel);
    if (!TrackHandle) {
        throw cet::exception("anatest") << " No rec::Track branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
    //findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);



    // Get handle for Vertices; also Assn's to Tracks
    auto VertexHandle = e.getHandle< std::vector<rec::Vertex> >(fVertexLabel);
    if (!VertexHandle) {
        throw cet::exception("anatest") << " No rec::Vertex branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
    //findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);


    // Get handle for CaloDigits
    auto RawHitHandle = e.getHandle< std::vector<gar::raw::CaloRawDigit> >(fRawCaloHitLabel);
    if (!RawHitHandle) {
        throw cet::exception("anatest") << " No :raw::CaloRawDigit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }


    // Get handle for CaloHits
    auto RecoHitHandle = e.getHandle< std::vector<rec::CaloHit> >(fCaloHitLabel);
    if (!RecoHitHandle) {
        throw cet::exception("anatest") << " No rec::CaloHit branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }


    // Get handle for CaloClusters; also Assn for matching tracks
    auto RecoClusterHandle = e.getHandle< std::vector<rec::Cluster> >(fClusterLabel);
    if (!RecoClusterHandle) {
        throw cet::exception("anatest") << " No rec::Cluster branch."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    //art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
    //findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>
    //                                   (RecoClusterHandle,e,fECALAssnLabel);



    // =============  Pull art handles =========================================
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();

    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
        for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
            if (mct.NeutrinoSet()) {
                simb::MCNeutrino nuw = mct.GetNeutrino();
                fNeutrinoType.push_back(nuw.Nu().PdgCode());
                fInteractionType.push_back(nuw.InteractionType());
            }
        }
    }

    // Need a non-constant backtracker instance, for now.
    cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
    cheat::BackTrackerCore*             bt = const_cast<cheat::BackTrackerCore*>(const_bt);

    bool const dumpMCP = false;
    sim::ParticleList* partList = bt->GetParticleList();
    int nMotherless = 0;
    for ( auto const& mcp : *MCPHandle ) {
        if (mcp.Mother() > 0) {
            simb::MCParticle* TPCeve = bt->FindTPCEve(mcp.TrackId());
            if ( TPCeve->TrackId() == mcp.TrackId() ) ++nMotherless;
        }
        if (dumpMCP) {
            std::cout << "TrackID: " << mcp.TrackId() << " is PDG: " << mcp.PdgCode() <<
                " with mother " << mcp.Mother() << " produced by process " << mcp.Process()
                << std::endl;
        }
    }
    if (dumpMCP) {
        std::cout << "Number of children: " << partList->size() << std::endl;
        std::cout << "Number of motherless children: " << nMotherless << std::endl;
    }




    std::vector<art::Ptr<rec::Hit>> allhits;
    art::PtrMaker<rec::Hit> makeHitPtr(e,HitHandle.id());
    for (size_t iHit=0; iHit<HitHandle->size(); ++iHit) {
        art::Ptr<rec::Hit> aPtr = makeHitPtr(iHit);
        allhits.push_back(aPtr);
    }
    for (simb::MCParticle mcp : *MCPHandle) {
		if ( mcp.Mother() > 0 ) continue;
	    if ( abs(mcp.PdgCode()) == 13 ) bt->ParticleToHits(&mcp,allhits);
	}





/*  for ( rec::Cluster cluster : *RecoClusterHandle ) {
        std::vector<std::pair<simb::MCParticle*,float>> whatMatches;
        whatMatches = bt->ClusterToMCParticles(&cluster);
        std::cout << "\nCluster No. " << cluster.getIDNumber() << " is made of MCParticles: " << std::endl;
        for (auto itr = whatMatches.begin(); itr!=whatMatches.end(); ++itr) {
            std::cout << "G4 track number " << itr->first->TrackId() << "\thas PDG code " <<
                itr->first->PdgCode() << "\tand its mother is G4 track " << itr->first->Mother()
                << "\tand energy fraction " << 100*(itr->second) << "%\n";
        }
    }

    std::vector<art::Ptr<rec::Cluster>> clusterCol;
    art::PtrMaker<rec::Cluster> makeClusterPtr(e,RecoClusterHandle.id());
    for (size_t iCluster=0; iCluster<RecoClusterHandle->size(); ++iCluster ) {
        art::Ptr<rec::Cluster> aPtr = makeClusterPtr(iCluster);
        clusterCol.push_back(aPtr);
    }
    for ( simb::MCParticle mcp : *MCPHandle ) {
        if (mcp.Mother() == 0) {
            std::cout << "Particle PDG: " << mcp.PdgCode() << " matches to\n";
            std::vector<art::Ptr<rec::Cluster>> clusterList =
                bt->MCParticleToClusters(&mcp, clusterCol);
            for (art::Ptr<rec::Cluster> iCluster : clusterList)
                std::cout << iCluster->getIDNumber() << std::endl;
        }
    } */





    return;
}



DEFINE_ART_MODULE(gar::anatest)
