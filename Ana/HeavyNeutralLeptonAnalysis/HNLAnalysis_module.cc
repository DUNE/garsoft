////////////////////////////////////////////////////////////////////////////////
// Module for heavy neutral lepton analysis
// 
// File: HNLAnalysis_module.cc
//
// Original version by Georgios Christodoulou, November 2020
//
// Some code is borrowed from anatree_module.cc
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
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/TrackIoniz.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "TTree.h"
#include "TH2.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TRandom3.h"

#include "CLHEP/Random/RandFlat.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace gar {
  
  typedef std::pair<float, std::string> P;

  class HNLAnalysis : public art::EDAnalyzer {
  public:
    explicit HNLAnalysis(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    HNLAnalysis(HNLAnalysis const &) = delete;
    HNLAnalysis(HNLAnalysis &&) = delete;
    HNLAnalysis & operator = (HNLAnalysis const &) = delete;
    HNLAnalysis & operator = (HNLAnalysis &&) = delete;
    
    virtual void beginJob() override;
    
    // Required functions.
    void analyze(art::Event const & e) override;
    
  private:
    void ClearVectors();
    void FillVectors(art::Event const & e);
    
    // Helpers copied from anatree_module
    void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate, float& forwardIonVal, float& backwardIonVal);
    float processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate);
    // Helper method for processOneDirection
    static bool lessThan_byE(std::pair<float,float> a, std::pair<float,float> b) {return a.first < b.first;}    

    // Get parameterized TPC PID
    std::vector< std::pair<int, float> > processPIDInfo(float p);

    // Get the parameterized ECal PID - this from CAF.cpp
    int GetECalPIDParameterized(int pdg, float ptrue);
    
    // Position of TPC from geometry service
    double ftpccentre[3];
    // TPC length
    double ftpclength;
    // TPC radius
    double ftpcradius;
    // ECal geometry
    double fECALBarrelInnerRadius;
    double fECALBarrelOuterRadius;
    double fECALEndcapInnerRadius;
    double fECALEndcapOuterRadius;
    double fECALStartX;
    double fECALEndX;
    
    // Input data labels
    std::string fGeantLabel;
    //std::string fGeantInstanceCalo; ///< Instance name ECAL for sdp::CaloDeposit
    //std::string fGeantInstanceMuID; ///< Instance name MuID for sdp::CaloDeposit

    std::string fTrackLabel;
    std::string fClusterLabel;
    std::string fECALAssnLabel;

    // Truncation parameter for dE/dx (average this fraction lowest readings)
    float fIonizTruncate;      ///< Default=1.00;
    // MCParticle to MC vertex match
    float fMatchMCPtoVertDist; ///< Default=roundoff

    // The analysis tree
    TTree *fTree;
    // The neutrino truth tree
    TTree *fTruthTree;
    // Configuration tree
    TTree *fconfig;
    
    // Geometry service
    const gar::geo::GeometryCore* fGeo;
    // Backtracker service
    cheat::BackTrackerCore* fBack;
    // PDG database from ROOT
    TDatabasePDG* pdgInstance;

    // Map of PID TH2 per momentum value
    CLHEP::HepRandomEngine &fEngine;  ///< random engine
    std::unordered_map<int, TH2F*> m_pidinterp;

    // The analysis tree. Use Rtypes.h here, as these data get used by root global event info
    Int_t	fRun;
    Int_t	fSubRun;
    Int_t	fEvent;
    Float_t     fWeight;

    // Number of TPC tracks
    Int_t       fNTPCTracks;
    Int_t       fNTPCECalTracks;
    
    // This analysis starts from reco'd tracks
    std::vector<ULong64_t>			fTrackIDNumber;
    std::vector<Float_t>			fTrackX;
    std::vector<Float_t>			fTrackY;
    std::vector<Float_t>			fTrackZ;
    std::vector<Float_t>			fTrackPX;
    std::vector<Float_t>			fTrackPY;
    std::vector<Float_t>			fTrackPZ;
    std::vector<Float_t>                        fTrackP;
    std::vector<Int_t>				fTrackQ;
    std::vector<Float_t>			fTrackLen;
    std::vector<Float_t>			fTrackChi2;
    std::vector<Int_t>				fNTPCClustersOnTrack;
    std::vector<Double_t>			fTrackTime;
    std::vector<Float_t>                        fTrackAvgIon;
    std::vector<Float_t>                        fTrackEndX;
    std::vector<Float_t>                        fTrackEndY;
    std::vector<Float_t>                        fTrackEndZ;
    std::vector<Float_t>                        fTrackEndPX;
    std::vector<Float_t>                        fTrackEndPY;
    std::vector<Float_t>                        fTrackEndPZ;
    std::vector<Float_t>                        fTrackEndP;
    std::vector<Int_t>                          fTrackEndQ;
    std::vector<Int_t>                          fTrackPID;
    // Which match to qualified MCParticles
    std::vector<Int_t>                          fMCPTrackID;
    std::vector<Int_t>                          fMCPDG;
    std::vector<Int_t>                          fMCPDGMother;
    std::vector<Float_t>                        fMCPStartX;
    std::vector<Float_t>                        fMCPStartY;
    std::vector<Float_t>                        fMCPStartZ;
    std::vector<Float_t>                        fMCPStartPX;
    std::vector<Float_t>                        fMCPStartPY;
    std::vector<Float_t>                        fMCPStartPZ;
    std::vector<Float_t>                        fMCPStartP;
    std::vector<Float_t>                        fMCPStartEndP;
    std::vector<std::string>                    fMCPProc;
    std::vector<std::string>                    fMCPEndProc;
    std::vector<Float_t>                        fMCPTime;
    // MC particle enters the ECal
    std::vector<Int_t>                          fMCPInECal;
    // Matched vertex id
    //std::vector<Int_t>                          fMCPFromVertexID;

    // ECal variables
    std::vector<Int_t>                          fNECalClustersForTrack;
    std::vector<Int_t>                          fClusterTrackIDMatched;
    std::vector<Int_t>                          fClusterNhits;
    std::vector<Float_t>                        fClusterEnergy;
    std::vector<Float_t>                        fClusterTime;
    std::vector<Float_t>                        fClusterTimeDiffFirstLast;
    std::vector<Float_t>                        fClusterX;
    std::vector<Float_t>                        fClusterY;
    std::vector<Float_t>                        fClusterZ;
    std::vector<Float_t>                        fClusterTheta;
    std::vector<Float_t>                        fClusterPhi;
    std::vector<Int_t>                          fClusterPID;
    std::vector<Float_t>                        fClusterMainAxisX;
    std::vector<Float_t>                        fClusterMainAxisY;
    std::vector<Float_t>                        fClusterMainAxisZ;

    // Neutrino truth
    std::vector<Int_t>                          f_NeutrinoType;
    std::vector<Int_t>                          f_CCNC;
    std::vector<Int_t>                          f_Mode;
    std::vector<Int_t>                          f_InteractionType;
    std::vector<Float_t>                        f_Q2;
    std::vector<Float_t>                        f_W;
    std::vector<Float_t>                        f_X;
    std::vector<Float_t>                        f_Y;
    std::vector<Float_t>                        f_Theta;
    std::vector<Float_t>                        f_T;
    std::vector<Float_t>                        f_MCVertexX;
    std::vector<Float_t>                        f_MCVertexY;
    std::vector<Float_t>                        f_MCVertexZ;
    std::vector<Float_t>                        f_MCnuE;
    std::vector<Int_t>                          f_TgtPDG;
    std::vector<Int_t>                          f_MCVertexInGasTPC;

    // MC truth
    std::vector<Int_t>                          fNeutrinoType;
    std::vector<Int_t>                          fCCNC;
    std::vector<Int_t>                          fMode;
    std::vector<Int_t>                          fInteractionType;
    std::vector<Float_t>                        fQ2;
    std::vector<Float_t>                        fW;
    std::vector<Float_t>                        fX;
    std::vector<Float_t>                        fY;
    std::vector<Float_t>                        fTheta;
    std::vector<Float_t>                        fT;
    std::vector<Float_t>                        fMCVertexX;
    std::vector<Float_t>                        fMCVertexY;
    std::vector<Float_t>                        fMCVertexZ;
    std::vector<Float_t>                        fMCnuPx;
    std::vector<Float_t>                        fMCnuPy;
    std::vector<Float_t>                        fMCnuPz;
    std::vector<Float_t>                        fMCnuE;
    std::vector<Int_t>                          fTgtPDG;
    std::vector<Int_t>                          fMCVertexInGasTPC;
    //std::vector<Int_t>                          fMCVertexID;

    // MC HNL truth decay products
    std::vector<Int_t>                          fMCP_TrackID;
    std::vector<Int_t>                          fMCP_PDG;
    std::vector<Float_t>                        fMCP_X;
    std::vector<Float_t>                        fMCP_Y;
    std::vector<Float_t>                        fMCP_Z;
    std::vector<Float_t>                        fMCP_T;
    std::vector<Float_t>                        fMCP_PX;
    std::vector<Float_t>                        fMCP_PY;
    std::vector<Float_t>                        fMCP_PZ;
    std::vector<Float_t>                        fMCP_E;
    std::vector<Int_t>                          fMCP_InGasTPC;
    //std::vector<Int_t>                          fMCP_VertexID;

  };

}

//==============================================================================
// constructor
gar::HNLAnalysis::HNLAnalysis(fhicl::ParameterSet const & p) : EDAnalyzer(p), fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed")) {

  fGeantLabel        = p.get<std::string>("GEANTLabel",   "geant");
  //fGeantInstanceCalo = p.get<std::string>("GEANTInstanceCalo","ECAL");
  //fGeantInstanceMuID = p.get<std::string>("GEANTInstanceMuID","MuID");

  fTrackLabel        = p.get<std::string>("TrackLabel",   "track");
  fClusterLabel      = p.get<std::string>("ClusterLabel", "calocluster");
  fECALAssnLabel     = p.get<std::string>("ECALAssnLabel","trkecalassn");

  fIonizTruncate     = p.get<float>("IonizTruncate",      1.00);

  float MCPtoVertDefault = 10.0*std::numeric_limits<Float_t>::epsilon();
  fMatchMCPtoVertDist    = p.get<float>("MatchMCPtoVertDist",MCPtoVertDefault);

  pdgInstance = TDatabasePDG::Instance();
  
  consumesMany<std::vector<simb::MCTruth> >();
  consumes<std::vector<simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<rec::Track> >(fTrackLabel);
  consumes<std::vector<rec::Cluster> >(fClusterLabel);
  consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

  // Calorimetry related
  //art::InputTag ecalgeanttag(fGeantLabel, fGeantInstanceCalo);
  //consumes<std::vector<gar::sdp::CaloDeposit> >(ecalgeanttag);

  // Muon system related
  //if (fGeo->HasMuonDetector()) {
  //art::InputTag muidgeanttag(fGeantLabel, fGeantInstanceMuID);
  //consumes<std::vector<gar::sdp::CaloDeposit> >(muidgeanttag);
  //}


} // end constructor

//==============================================================================
void gar::HNLAnalysis::beginJob() {
  
  // Backtracker
  cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
  fBack = const_cast<cheat::BackTrackerCore*>(const_bt);

  // Geometry
  fGeo = gar::providerFrom<geo::GeometryGAr>();
  ftpccentre[0]          = fGeo->TPCXCent();
  ftpccentre[1]          = fGeo->TPCYCent();
  ftpccentre[2]          = fGeo->TPCZCent();
  ftpcradius             = fGeo->TPCRadius();
  ftpclength             = fGeo->TPCLength();
  fECALBarrelInnerRadius = fGeo->GetECALInnerBarrelRadius();
  fECALBarrelOuterRadius = fGeo->GetECALOuterBarrelRadius();
  fECALEndcapInnerRadius = fGeo->GetECALInnerEndcapRadius();
  fECALEndcapOuterRadius = fGeo->GetECALOuterEndcapRadius();
  fECALStartX            = fGeo->GetECALEndcapStartX();
  fECALEndX              = fGeo->GetECALEndcapOuterX();
  
  art::ServiceHandle<art::TFileService> tfs;

  // Save TPC/ECal geometry
  fconfig = tfs->make<TTree>("config","config");
  fconfig->Branch("tpccentre",                          &ftpccentre,                    "tpccentre[3]/D");
  fconfig->Branch("tpcradius",                          &ftpcradius,                    "tpcradius/D");
  fconfig->Branch("tpclength",                          &ftpclength,                    "tpclength/D");
  fconfig->Branch("ECALBarrelInnerRadius",              &fECALBarrelInnerRadius,        "ECALBarrelInnerRadius/D");
  fconfig->Branch("ECALBarrelOuterRadius",              &fECALBarrelOuterRadius,        "ECALBarrelOuterRadius/D");
  fconfig->Branch("ECALEndcapInnerRadius",              &fECALEndcapInnerRadius,        "ECALEndcapInnerRadius/D");
  fconfig->Branch("ECALEndcapOuterRadius",              &fECALEndcapOuterRadius,        "ECALEndcapOuterRadius/D");
  fconfig->Branch("ECALStartX",                         &fECALStartX,                   "ECALStartX/D");
  fconfig->Branch("ECALEndX",                           &fECALEndX,                     "ECALEndX/D");
  fconfig->Fill();

  // Analysis tree with reco 
  fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

  fTree->Branch("Run",					&fRun,				"Run/I");
  fTree->Branch("SubRun",				&fSubRun,			"SubRun/I");
  fTree->Branch("Event",				&fEvent,			"Event/I");
  fTree->Branch("Weight",                               &fWeight,                       "Weight/F");
  fTree->Branch("NTPCTracks",                           &fNTPCTracks,                   "NTPCTracks/I");
  fTree->Branch("NTPCECalTracks",                       &fNTPCECalTracks,               "NTPCECalTracks/I");

  fTree->Branch("TrackIDNumber",			&fTrackIDNumber);
  fTree->Branch("TrackX",				&fTrackX);
  fTree->Branch("TrackY", 				&fTrackY);
  fTree->Branch("TrackZ", 				&fTrackZ);
  fTree->Branch("TrackPX",				&fTrackPX);
  fTree->Branch("TrackPY",				&fTrackPY);
  fTree->Branch("TrackPZ",				&fTrackPZ);
  fTree->Branch("TrackP",                               &fTrackP);
  fTree->Branch("TrackQ", 				&fTrackQ);
  fTree->Branch("TrackLen",				&fTrackLen);
  fTree->Branch("TrackChi2",				&fTrackChi2);
  fTree->Branch("NTPCClustersOnTrack",	                &fNTPCClustersOnTrack);
  fTree->Branch("TrackTime",				&fTrackTime);
  fTree->Branch("TrackAvgIon",                          &fTrackAvgIon);
  fTree->Branch("TrackEndX",                            &fTrackEndX);
  fTree->Branch("TrackEndY",                            &fTrackEndY);
  fTree->Branch("TrackEndZ",                            &fTrackEndZ);
  fTree->Branch("TrackEndPX",                           &fTrackEndPX);
  fTree->Branch("TrackEndPY",                           &fTrackEndPY);
  fTree->Branch("TrackEndPZ",                           &fTrackEndPZ);
  fTree->Branch("TrackEndP",                            &fTrackEndP);
  fTree->Branch("TrackEndQ",                            &fTrackEndQ);
  fTree->Branch("TrackPID",                             &fTrackPID);

  fTree->Branch("NECalClustersForTrack",                &fNECalClustersForTrack);
  fTree->Branch("ClusterTrackIDMatched",                &fClusterTrackIDMatched);
  fTree->Branch("ClusterNhits",                         &fClusterNhits);
  fTree->Branch("ClusterEnergy",                        &fClusterEnergy);
  fTree->Branch("ClusterTime",                          &fClusterTime);
  fTree->Branch("ClusterTimeDiffFirstLast",             &fClusterTimeDiffFirstLast);
  fTree->Branch("ClusterX",                             &fClusterX);
  fTree->Branch("ClusterY",                             &fClusterY);
  fTree->Branch("ClusterZ",                             &fClusterZ);
  fTree->Branch("ClusterTheta",                         &fClusterTheta);
  fTree->Branch("ClusterPhi",                           &fClusterPhi);
  fTree->Branch("ClusterPID",                           &fClusterPID);
  fTree->Branch("ClusterMainAxisX",                     &fClusterMainAxisX);
  fTree->Branch("ClusterMainAxisY",                     &fClusterMainAxisY);
  fTree->Branch("ClusterMainAxisZ",                     &fClusterMainAxisZ);

  fTree->Branch("MCPTrackID",                           &fMCPTrackID);
  fTree->Branch("PDG",                                  &fMCPDG);
  fTree->Branch("PDGMother",                            &fMCPDGMother);
  fTree->Branch("MCPStartX",                            &fMCPStartX);
  fTree->Branch("MCPStartY",                            &fMCPStartY);
  fTree->Branch("MCPStartZ",                            &fMCPStartZ);
  fTree->Branch("MCPStartPX",                           &fMCPStartPX);
  fTree->Branch("MCPStartPY",                           &fMCPStartPY);
  fTree->Branch("MCPStartPZ",                           &fMCPStartPZ);
  fTree->Branch("MCPStartP",                            &fMCPStartP);
  fTree->Branch("MCPStartEndP",                         &fMCPStartEndP);
  fTree->Branch("MCPProc",                              &fMCPProc);
  fTree->Branch("MCPEndProc",                           &fMCPEndProc);
  fTree->Branch("MCPTime",                              &fMCPTime);
  fTree->Branch("MCPInECal",                            &fMCPInECal);
  //fTree->Branch("MCPFromVertexID",                      &fMCPFromVertexID);

  // Neutrino truth
  fTree->Branch("Nu_Type",                              &f_NeutrinoType);
  fTree->Branch("Nu_CCNC",                              &f_CCNC);
  fTree->Branch("Nu_Mode",                              &f_Mode);
  fTree->Branch("Nu_InterT",                            &f_InteractionType);
  fTree->Branch("Nu_Q2",                                &f_Q2);
  fTree->Branch("Nu_W",                                 &f_W);
  fTree->Branch("Nu_X",                                 &f_X);
  fTree->Branch("Nu_Y",                                 &f_Y);
  fTree->Branch("Nu_nuTheta",                           &f_Theta);
  fTree->Branch("Nu_nuVertX",                           &f_MCVertexX);
  fTree->Branch("Nu_nuVertY",                           &f_MCVertexY);
  fTree->Branch("Nu_nuVertZ",                           &f_MCVertexZ);
  fTree->Branch("Nu_E",                                 &f_MCnuE);
  fTree->Branch("Nu_TgtPDG",                            &f_TgtPDG);
  fTree->Branch("Nu_VertexInGasTPC",                    &f_MCVertexInGasTPC);

  // Truth tree
  fTruthTree = tfs->make<TTree>("GArTruthTree","GArTruthTree");
  fTruthTree->Branch("Run",                             &fRun,                          "Run/I");
  fTruthTree->Branch("SubRun",                          &fSubRun,                       "SubRun/I");
  fTruthTree->Branch("Event",                           &fEvent,                        "Event/I");
  fTruthTree->Branch("Weight",                          &fWeight,                       "Weight/F");

  fTruthTree->Branch("Nu_Type",                         &fNeutrinoType);
  fTruthTree->Branch("Nu_CCNC",                         &fCCNC);
  fTruthTree->Branch("Nu_Mode",                         &fMode);
  fTruthTree->Branch("Nu_InterT",                       &fInteractionType);
  fTruthTree->Branch("Nu_Q2",                           &fQ2);
  fTruthTree->Branch("Nu_W",                            &fW);
  fTruthTree->Branch("Nu_X",                            &fX);
  fTruthTree->Branch("Nu_Y",                            &fY);
  fTruthTree->Branch("Nu_Theta",                        &fTheta);
  fTruthTree->Branch("Nu_VertX",                        &fMCVertexX);
  fTruthTree->Branch("Nu_VertY",                        &fMCVertexY);
  fTruthTree->Branch("Nu_VertZ",                        &fMCVertexZ);
  fTruthTree->Branch("Nu_Px",                           &fMCnuPx);
  fTruthTree->Branch("Nu_Py",                           &fMCnuPy);
  fTruthTree->Branch("Nu_Pz",                           &fMCnuPz);
  fTruthTree->Branch("Nu_E",                            &fMCnuE);
  fTruthTree->Branch("Nu_TgtPDG",                       &fTgtPDG);
  fTruthTree->Branch("Nu_VertexInGasTPC",               &fMCVertexInGasTPC);
  //fTruthTree->Branch("MCVertexID",                      &fMCVertexID);

  fTruthTree->Branch("MCP_TrackID",                     &fMCP_TrackID);
  fTruthTree->Branch("MCP_PDG",                         &fMCP_PDG);
  fTruthTree->Branch("MCP_X",                           &fMCP_X);
  fTruthTree->Branch("MCP_Y",                           &fMCP_Y);
  fTruthTree->Branch("MCP_Z",                           &fMCP_Z);
  fTruthTree->Branch("MCP_T",                           &fMCP_T);
  fTruthTree->Branch("MCP_PX",                          &fMCP_PX);
  fTruthTree->Branch("MCP_PY",                          &fMCP_PY);
  fTruthTree->Branch("MCP_PZ",                          &fMCP_PZ);
  fTruthTree->Branch("MCP_E",                           &fMCP_E);
  fTruthTree->Branch("MCP_InGasTPC",                    &fMCP_InGasTPC);
  //fTruthTree->Branch("MCP_VertexID",                    &fMCP_VertexID);

  // Get histograms for TPC PID parameterization
  std::string filename = "${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
  TFile infile(filename.c_str(), "READ");

  m_pidinterp.clear();
  char str[11];
  for(int q = 0; q < 501; ++q){
    sprintf(str, "%d", q);
    std::string s = "pidmatrix";
    s.append(str);
    // read the 500 histograms one by one; each histogram is a
    // 6 by 6 matrix of probabilities for a given momentum value
    m_pidinterp.insert( std::make_pair(q, (TH2F*) infile.Get(s.c_str())->Clone("pidinterp")) );
  }

  return;

}  // End of :HNLAnalysis::beginJob

//==============================================================================
void gar::HNLAnalysis::analyze(art::Event const & event) {
  
  /// Is the backtracker good to go?
  //cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
  //fBack = const_cast<cheat::BackTrackerCore*>(const_bt);
  //if ( !fBack->HasTracks() || !fBack->HasClusters() ) {
  //throw cet::exception("HNLAnalysis") << " Not enough backtracker info"
  //<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  //}
  
  ClearVectors();
  FillVectors(event);

  return;

}

//==============================================================================
void gar::HNLAnalysis::ClearVectors() {

  fNTPCECalTracks       = 0;

  // clear out all our vectors  
  fTrackIDNumber.clear();
  fTrackX.clear();
  fTrackY.clear();
  fTrackZ.clear();
  fTrackPX.clear();
  fTrackPY.clear();
  fTrackPZ.clear();
  fTrackP.clear();
  fTrackQ.clear();
  fTrackLen.clear();
  fTrackChi2.clear();
  fNTPCClustersOnTrack.clear();
  fTrackTime.clear();
  fTrackAvgIon.clear();
  fTrackEndX.clear();
  fTrackEndY.clear();
  fTrackEndZ.clear();
  fTrackEndPX.clear();
  fTrackEndPY.clear();
  fTrackEndPZ.clear();
  fTrackEndP.clear();
  fTrackEndQ.clear();
  fTrackPID.clear();

  fNECalClustersForTrack.clear();
  fClusterTrackIDMatched.clear();
  fClusterNhits.clear();
  fClusterEnergy.clear();
  fClusterTime.clear();
  fClusterTimeDiffFirstLast.clear();
  fClusterX.clear();
  fClusterY.clear();
  fClusterZ.clear();
  fClusterTheta.clear();
  fClusterPhi.clear();
  fClusterPID.clear();
  fClusterMainAxisX.clear();
  fClusterMainAxisY.clear();
  fClusterMainAxisZ.clear();
  
  fMCPTrackID.clear();
  fMCPDG.clear();
  fMCPDGMother.clear();
  fMCPStartX.clear();
  fMCPStartY.clear();
  fMCPStartZ.clear();
  fMCPStartPX.clear();
  fMCPStartPY.clear();
  fMCPStartPZ.clear();
  fMCPStartP.clear();
  fMCPStartEndP.clear();
  fMCPProc.clear();
  fMCPEndProc.clear();
  fMCPTime.clear();
  fMCPInECal.clear();
  //fMCPFromVertexID.clear();

  f_NeutrinoType.clear();
  f_CCNC.clear();
  f_Mode.clear();
  f_InteractionType.clear();
  f_Q2.clear();
  f_W.clear();
  f_X.clear();
  f_Y.clear();
  f_Theta.clear();
  f_MCVertexX.clear();
  f_MCVertexY.clear();
  f_MCVertexZ.clear();
  f_MCnuE.clear();
  f_TgtPDG.clear();
  f_MCVertexInGasTPC.clear();

  fNeutrinoType.clear();
  fCCNC.clear();
  fMode.clear();
  fInteractionType.clear();
  fQ2.clear();
  fW.clear();
  fX.clear();
  fY.clear();
  fTheta.clear();
  fMCVertexX.clear();
  fMCVertexY.clear();
  fMCVertexZ.clear();
  fMCnuPx.clear();
  fMCnuPy.clear();
  fMCnuPz.clear();
  fMCnuE.clear();
  fTgtPDG.clear();
  fMCVertexInGasTPC.clear();
  //fMCVertexID.clear();

  fMCP_TrackID.clear();
  fMCP_PDG.clear();
  fMCP_X.clear();
  fMCP_Y.clear();
  fMCP_Z.clear();
  fMCP_T.clear();
  fMCP_PX.clear();
  fMCP_PY.clear();
  fMCP_PZ.clear();
  fMCP_E.clear();
  fMCP_InGasTPC.clear();
  //fMCP_VertexID.clear();
  
  return;

} // end :HNLAnalysis::ClearVectors

//==============================================================================
void gar::HNLAnalysis::FillVectors(art::Event const& event) {

  // =============  Get art handles ==========================================
  // Get handles for Tracks, clusters, associations between them
  auto TrackHandle = event.getHandle< std::vector<rec::Track> >(fTrackLabel);
  if (!TrackHandle) {
    throw cet::exception("HNLAnalysis") << " No rec::Track"
					<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }
  auto RecoClusterHandle = event.getHandle< std::vector<rec::Cluster> >(fClusterLabel);
  if (!RecoClusterHandle) {
    throw cet::exception("HNLAnalysis") << " No rec::Cluster branch."
					<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  auto TrackIonHandle = event.getHandle< std::vector<rec::TrackIoniz> >(fTrackLabel);
  if (!TrackIonHandle) {
    throw cet::exception("HNLAnalysis") << " No rec::TrackIoniz branch."
					<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
  findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,event,fTrackLabel);

  art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
  findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(RecoClusterHandle,event,fECALAssnLabel);

  // Get handles for MCinfo, also good for MCPTrajectory.  Want to get all MCTruths, regardless of generator label.

  auto mctruthHandles = event.getMany<std::vector<simb::MCTruth> >();
  if (mctruthHandles.size()!=1) {
    throw cet::exception("HNLAnalysis") << " Need just 1 simb::MCTruth"
					<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }
  
  auto MCPHandle = event.getHandle< std::vector<simb::MCParticle> >(fGeantLabel);
  if (!MCPHandle) {
    throw cet::exception("HNLAnalysis") << " No simb::MCParticle"
					<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  // =============  Start analysis =========================================
  fRun    = event.run();
  fSubRun = event.subRun();
  fEvent  = event.id().event();
  fWeight = 1.0;

  // Save all neutrino vertices
  //int nuvertexid = 0;
  for (size_t imchl = 0; imchl < mctruthHandles.size(); ++imchl) {
    for ( auto const& mct : (*mctruthHandles.at(imchl)) ) {
      if(mct.NeutrinoSet()){ // Standard neutrino
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
	fMCVertexX.push_back(nuw.Nu().EndX());
	fMCVertexY.push_back(nuw.Nu().EndY());
	fMCVertexZ.push_back(nuw.Nu().EndZ());
	fMCnuPx.push_back(nuw.Nu().Px());
	fMCnuPy.push_back(nuw.Nu().Py());
	fMCnuPz.push_back(nuw.Nu().Pz());
	fMCnuE.push_back( sqrt(nuw.Nu().Px()*nuw.Nu().Px() + nuw.Nu().Py()*nuw.Nu().Py() + nuw.Nu().Pz()*nuw.Nu().Pz()) );
	fTgtPDG.push_back(nuw.Target());
	TVector3 nupos(nuw.Nu().EndX(), nuw.Nu().EndY(), nuw.Nu().EndZ());
	fMCVertexInGasTPC.push_back(fGeo->PointInGArTPC(nupos));
	//fMCVertexID.push_back(nuvertexid++);
      }
      else{
	fNeutrinoType.push_back(2000020000);
        fCCNC.push_back(0);
        fMode.push_back(-1);
        fInteractionType.push_back(-1);
        fQ2.push_back(0);
        fW.push_back(0);
        fX.push_back(0);
        fY.push_back(0);
        fTheta.push_back(0);
        fMCVertexX.push_back(0);
        fMCVertexY.push_back(0);
        fMCVertexZ.push_back(0);
        fMCnuPx.push_back(0);
        fMCnuPy.push_back(0);
        fMCnuPz.push_back(0);
	fMCnuE.push_back(0);
        fTgtPDG.push_back(0);
        fMCVertexInGasTPC.push_back(1);
        //fMCVertexID.push_back(nuvertexid++);
      }

      int nparticles = mct.NParticles();
      for(int i = 0; i < nparticles; i++){
	simb::MCParticle mcpart = mct.GetParticle(i);

	if(mcpart.Weight() != fWeight) fWeight = mcpart.Weight();
	
	const TLorentzVector& position = mcpart.Position(0);
	const TLorentzVector& momentum = mcpart.Momentum(0);
	
	fMCP_PDG.push_back(mcpart.PdgCode());
	fMCP_TrackID.push_back(mcpart.TrackId());
	fMCP_X.push_back(position.X());
	fMCP_Y.push_back(position.Y());
	fMCP_Z.push_back(position.Z());
	fMCP_T.push_back(mcpart.T());
	fMCP_PX.push_back(momentum.Px());
	fMCP_PY.push_back(momentum.Py());
	fMCP_PZ.push_back(momentum.Pz());
	fMCP_E.push_back(mcpart.E());
	TVector3 nupos(position.X(), position.Y(), position.Z());
	fMCP_InGasTPC.push_back(fGeo->PointInGArTPC(nupos));
	//fMCP_VertexID.push_back(nuvertexid);	
      }
    }
  }
  
  // Fill the truth tree
  fTruthTree->Fill();

  // Make the usual map for the MC info
  typedef int TrkId;
  std::unordered_map<TrkId, Int_t> TrackIdToIndex;
  Int_t index = 0;
  for ( auto const& mcp : (*MCPHandle) ) {
    int TrackId = mcp.TrackId();
    TrackIdToIndex[TrackId] = index++;
  }

  // Loop over TPC tracks
  size_t iTrack = 0;
  for ( auto const& track : (*TrackHandle) ) {

    std::vector<std::pair<simb::MCParticle*,float>> MCsfromTrack;
    gar::rec::Track* nonconst_track = const_cast<gar::rec::Track*>(&track);
    MCsfromTrack = fBack->TrackToMCParticles(nonconst_track);
    int nMCsfromTrack = MCsfromTrack.size();
    if(nMCsfromTrack==0){
      continue;
    }
    
    simb::MCParticle theMCPart;
    for (int iMCfromTrack =0; iMCfromTrack<nMCsfromTrack; ++iMCfromTrack) {
      // Plausible MCParticle to make this track?
      theMCPart = *(MCsfromTrack[iMCfromTrack].first);
      TParticlePDG* particle = pdgInstance->GetParticle(theMCPart.PdgCode());
      if ( particle->Charge() == 0.0 ) continue;
      if ( particle->Stable() == 0   ) continue;
      //if ( particle->PdgCode() == 22 || particle->PdgCode() == 2112) continue;
      break;			// Take highest fraction candidate
    }
    
    // Does theMCParticle go out to the ECAL?
    int nTraj = theMCPart.NumberTrajectoryPoints();
    TVector3 doink;			
    bool whackThatECAL = false;
    for (int iTraj=0; iTraj<nTraj; ++iTraj) {
      doink.SetXYZ(theMCPart.Vx(iTraj),theMCPart.Vy(iTraj),theMCPart.Vz(iTraj));
      if ( fGeo->PointInECALBarrel(doink) || fGeo->PointInECALEndcap(doink) ) {
	whackThatECAL = true;
	break;
      }
    }

    fMCPInECal.push_back(whackThatECAL);
    
    // Which end of this track corresponds to the start of the reco track?
    const TLorentzVector& positionMCP = theMCPart.Position(0);
    float distStart = std::hypot(track.Vertex()[0] -positionMCP[0],
				 track.Vertex()[1] -positionMCP[1],
				 track.Vertex()[2] -positionMCP[2]);
    float distEnd   = std::hypot(track.End()[0]	   -positionMCP[0],
				 track.End()[1]	   -positionMCP[1],
				 track.End()[2]	   -positionMCP[2]);
    rec::TrackEnd kate = (distStart<distEnd) ? rec::TrackEndBeg : rec::TrackEndEnd;

    // To store parameterized TPC PID
    std::vector< std::pair<int, float> > pid;

    // Ready to record some info about the track & matching MC info
    fTrackIDNumber.push_back(track.getIDNumber());
    if(kate==rec::TrackEndBeg){
      fTrackX.push_back     (track.Vertex()[0]);
      fTrackY.push_back     (track.Vertex()[1]);
      fTrackZ.push_back     (track.Vertex()[2]);
      fTrackPX.push_back    (track.Momentum_beg()*track.VtxDir()[0]);
      fTrackPY.push_back    (track.Momentum_beg()*track.VtxDir()[1]);
      fTrackPZ.push_back    (track.Momentum_beg()*track.VtxDir()[2]);
      fTrackP.push_back     (track.Momentum_beg());
      fTrackQ.push_back     (track.ChargeBeg());
      fTrackLen.push_back   (track.LengthForward());
      fTrackChi2.push_back  (track.ChisqForward());

      fTrackEndX.push_back  (track.End()[0]);
      fTrackEndY.push_back  (track.End()[1]);
      fTrackEndZ.push_back  (track.End()[2]);
      fTrackEndPX.push_back (track.Momentum_end()*track.EndDir()[0]);
      fTrackEndPY.push_back (track.Momentum_end()*track.EndDir()[1]);
      fTrackEndPZ.push_back (track.Momentum_end()*track.EndDir()[2]);
      fTrackEndP.push_back  (track.Momentum_end());
      fTrackEndQ.push_back  (track.ChargeEnd());

      // Add the TPC PID information based on Tom's parametrization
      pid = processPIDInfo( track.Momentum_beg() );
    } 
    else{
      fTrackX.push_back     (track.End()[0]);
      fTrackY.push_back     (track.End()[1]);
      fTrackZ.push_back     (track.End()[2]);
      fTrackPX.push_back    (track.Momentum_end()*track.EndDir()[0]);
      fTrackPY.push_back    (track.Momentum_end()*track.EndDir()[1]);
      fTrackPZ.push_back    (track.Momentum_end()*track.EndDir()[2]);
      fTrackP.push_back     (track.Momentum_end());
      fTrackQ.push_back     (track.ChargeEnd());
      fTrackLen.push_back   (track.LengthBackward());		
      fTrackChi2.push_back  (track.ChisqBackward());

      fTrackEndX.push_back  (track.Vertex()[0]);
      fTrackEndY.push_back  (track.Vertex()[1]);
      fTrackEndZ.push_back  (track.Vertex()[2]);
      fTrackEndPX.push_back (track.Momentum_beg()*track.VtxDir()[0]);
      fTrackEndPY.push_back (track.Momentum_beg()*track.VtxDir()[1]);
      fTrackEndPZ.push_back (track.Momentum_beg()*track.VtxDir()[2]);
      fTrackEndP.push_back  (track.Momentum_beg());
      fTrackEndQ.push_back  (track.ChargeBeg());

      // Add the TPC PID information based on Tom's parametrization
      pid = processPIDInfo( track.Momentum_end() );
    }
    fNTPCClustersOnTrack.push_back(track.NHits());
    fTrackTime.push_back(track.Time());

    // Check TPC PID
    float pidhypo[5] = {5.0, 5.0, 5.0, 5.0, 5.0};
    for(size_t ipid = 0; ipid < pid.size(); ipid++){
      if(pid.at(ipid).first == 1000010020) continue;
      if     (pid.at(ipid).first == 13   && pid.at(ipid).second < pidhypo[0]) pidhypo[0] = pid.at(ipid).second;
      else if(pid.at(ipid).first == 211  && pid.at(ipid).second < pidhypo[1]) pidhypo[1] = pid.at(ipid).second;
      else if(pid.at(ipid).first == 321  && pid.at(ipid).second < pidhypo[2]) pidhypo[2] = pid.at(ipid).second;
      else if(pid.at(ipid).first == 2212 && pid.at(ipid).second < pidhypo[3]) pidhypo[3] = pid.at(ipid).second;
      else if(pid.at(ipid).first == 11   && pid.at(ipid).second < pidhypo[4]) pidhypo[4] = pid.at(ipid).second;     
    }

    for(int ihypo = 0; ihypo < 5; ihypo++){
      if(pidhypo[ihypo] == 5.0) pidhypo[ihypo] = 0.0;
    }

    // Backtracker is not perfect 
    int trpdg = abs(theMCPart.PdgCode());
    if(trpdg == 22) trpdg = 11;
    else if(trpdg == 2112) trpdg = 2212;

    CLHEP::RandFlat FlatRand(fEngine);
    if( (trpdg == 11   and pidhypo[4] == 1) ||
	(trpdg == 13   and pidhypo[0] == 1) ||
	(trpdg == 211  and pidhypo[1] == 1) ||
	(trpdg == 321  and pidhypo[2] == 1) ||
	(trpdg == 2212 and pidhypo[3] == 1) ) { 
      fTrackPID.push_back(trpdg);
    }
    else{
      int pidpdg = -1;
      for(int ihypo = 0; ihypo < 5; ihypo++){
	if(pidhypo[ihypo] == 1 || pidhypo[ihypo] == 0) continue;
	// Throw a random number between 0 and 1
	float random_number = FlatRand.fire();
	if(random_number < pidhypo[ihypo]){ 
	  pidpdg = ihypo;
	  break;
	}
      }
      
      if     (pidpdg == 0) fTrackPID.push_back(13);
      else if(pidpdg == 1) fTrackPID.push_back(211);
      else if(pidpdg == 2) fTrackPID.push_back(321);
      else if(pidpdg == 3) fTrackPID.push_back(2212);
      else if(pidpdg == 4) fTrackPID.push_back(11);
      else fTrackPID.push_back(0);
    }

    // Check if the track goes in the ECal and save the ECal info
    std::vector<ULong64_t> clusterID;
    size_t iCluster = 0;
    for ( auto const& cluster : (*RecoClusterHandle) ) {
      int nCALedTracks(0);
      if ( findManyCALTrackEnd->isValid() ) {
	nCALedTracks = findManyCALTrackEnd->at(iCluster).size();
      }

      for (int iCALedTrack=0; iCALedTrack<nCALedTracks; ++iCALedTrack) {
	rec::Track temptrack  = *(findManyCALTrackEnd->at(iCluster).at(iCALedTrack));
	if(temptrack.getIDNumber() == track.getIDNumber()){
	  clusterID.push_back(cluster.getIDNumber());
	}
      }
      iCluster++;
    }

    //std::sort(clusterID.begin(), clusterID.end());
    //auto iter = std::unique(clusterID.begin(), clusterID.end());
    //clusterID.erase(iter, clusterID.end());

    fNECalClustersForTrack.push_back(clusterID.size());
    if(clusterID.size() > 0){
      fNTPCECalTracks++;
      for ( auto const& cluster : (*RecoClusterHandle) ) {
	for(long unsigned int i = 0; i < clusterID.size(); i++){
	  if(clusterID.at(i) == cluster.getIDNumber()){
	    fClusterTrackIDMatched.push_back(track.getIDNumber());
	    fClusterNhits.push_back(cluster.CalorimeterHits().size());
	    fClusterEnergy.push_back(cluster.Energy());
	    fClusterTime.push_back(cluster.Time());
	    fClusterTimeDiffFirstLast.push_back(cluster.TimeDiffFirstLast());
	    fClusterX.push_back(cluster.Position()[0]);
	    fClusterY.push_back(cluster.Position()[1]);
	    fClusterZ.push_back(cluster.Position()[2]);
	    fClusterTheta.push_back(cluster.ITheta());
	    fClusterPhi.push_back(cluster.IPhi());
	    //fClusterPID.push_back(cluster.ParticleID());
	    fClusterPID.push_back(GetECalPIDParameterized( abs(theMCPart.PdgCode()), sqrt(theMCPart.EndPx()*theMCPart.EndPx() + theMCPart.EndPy()*theMCPart.EndPy() + theMCPart.EndPz()*theMCPart.EndPz())) );
	    fClusterMainAxisX.push_back(cluster.EigenVectors()[0]);
	    fClusterMainAxisY.push_back(cluster.EigenVectors()[1]);
	    fClusterMainAxisZ.push_back(cluster.EigenVectors()[2]);
	  }	
	}
      }
    }

    // Save dE/dx
    if (findIonization->isValid()) {
      // No calibration for now.  Someday this should all be in reco
      rec::TrackIoniz ionization = *(findIonization->at(iTrack));
      float avgIonF, avgIonB;
      processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
      if (kate==rec::TrackEndBeg) {
	fTrackAvgIon.push_back( avgIonF );
      }
      else{
	fTrackAvgIon.push_back( avgIonB );
      }
    } 
    else {
      // must push_back something so that fTrackAvgIon is of correct size.
      fTrackAvgIon.push_back( 0.0 );
    }
    iTrack++;

    fMCPDG.push_back(theMCPart.PdgCode());
    TrkId momTrkId = theMCPart.Mother();
    int momPDG = 0;
    if (momTrkId>0) {
      int momIndex = TrackIdToIndex[momTrkId];
      momPDG   = (*MCPHandle).at(momIndex).PdgCode();
    }
    fMCPDGMother.push_back(momPDG);
    
    const TLorentzVector& position = theMCPart.Position(0);
    fMCPStartX.push_back(position.X());
    fMCPStartY.push_back(position.Y());
    fMCPStartZ.push_back(position.Z());
    const TLorentzVector& momentum = theMCPart.Momentum(0);
    fMCPStartPX.push_back(momentum.Px());
    fMCPStartPY.push_back(momentum.Py());
    fMCPStartPZ.push_back(momentum.Pz());
    fMCPStartP.push_back(theMCPart.P());
    fMCPStartEndP.push_back(sqrt(theMCPart.EndPx()*theMCPart.EndPx() + theMCPart.EndPy()*theMCPart.EndPy() + theMCPart.EndPz()*theMCPart.EndPz()));

    fMCPProc.push_back(theMCPart.Process());
    fMCPEndProc.push_back(theMCPart.EndProcess());
    fMCPTime.push_back(theMCPart.T());
    fMCPTrackID.push_back(theMCPart.TrackId());

    // Find the vertex that this MCParticle is coming from
    bool vfound = false;
    for (size_t imchl = 0; imchl < mctruthHandles.size(); ++imchl) {
      for ( auto const& mct : (*mctruthHandles.at(imchl)) ) {
        if (mct.NeutrinoSet()) {
	  simb::MCNeutrino nuw = mct.GetNeutrino();
	  Float_t vertX = nuw.Nu().EndX();
          Float_t vertY = nuw.Nu().EndY();
          Float_t vertZ = nuw.Nu().EndZ();
          Float_t dist  = std::hypot(position.X()-vertX,position.Y()-vertY,position.Z()-vertZ);

	  if ( dist <= fMatchMCPtoVertDist ) { // First neutrino found
	    f_NeutrinoType.push_back(nuw.Nu().PdgCode());
	    f_CCNC.push_back(nuw.CCNC());
	    f_Mode.push_back(nuw.Mode());
	    f_InteractionType.push_back(nuw.InteractionType());
	    f_Q2.push_back(nuw.QSqr());
	    f_W.push_back(nuw.W());
	    f_X.push_back(nuw.X());
	    f_Y.push_back(nuw.Y());
	    f_Theta.push_back(nuw.Theta());
	    f_MCVertexX.push_back(nuw.Nu().EndX());
	    f_MCVertexY.push_back(nuw.Nu().EndY());
	    f_MCVertexZ.push_back(nuw.Nu().EndZ());
	    f_MCnuE.push_back( sqrt(nuw.Nu().Px()*nuw.Nu().Px() + nuw.Nu().Py()*nuw.Nu().Py() + nuw.Nu().Pz()*nuw.Nu().Pz()) );
	    f_TgtPDG.push_back(nuw.Target());
	    TVector3 nupos(nuw.Nu().EndX(), nuw.Nu().EndY(), nuw.Nu().EndZ());
	    f_MCVertexInGasTPC.push_back(fGeo->PointInGArTPC(nupos));
	    vfound = true;

	    break;
	  }
	}
      }
    }
    
    // Fill something if no neutrino is found 
    if(!vfound){
      f_NeutrinoType.push_back(0);
      f_CCNC.push_back(0);
      f_Mode.push_back(-1);
      f_InteractionType.push_back(-1);
      f_Q2.push_back(0);
      f_W.push_back(0);
      f_X.push_back(0);
      f_Y.push_back(0);
      f_Theta.push_back(0);
      f_MCVertexX.push_back(0);
      f_MCVertexY.push_back(0);
      f_MCVertexZ.push_back(0);
      f_MCnuE.push_back(0);
      f_TgtPDG.push_back(0);
      f_MCVertexInGasTPC.push_back(0);
    }

    /*
    // Find the vertex that this MCParticle is coming from
    nuvertexid = 0; int hnlvertexid = 0; int count = -1;
    for (size_t imchl = 0; imchl < mctruthHandles.size(); ++imchl) {
      for ( auto const& mct : (*mctruthHandles.at(imchl)) ) {
	if (mct.NeutrinoSet()) {
	  simb::MCNeutrino nuw = mct.GetNeutrino();
	  Float_t vertX = nuw.Nu().EndX();
	  Float_t vertY = nuw.Nu().EndY();
	  Float_t vertZ = nuw.Nu().EndZ();
	  Float_t dist  = std::hypot(position.X()-vertX,position.Y()-vertY,position.Z()-vertZ);
	  nuvertexid++;
	  if ( dist <= fMatchMCPtoVertDist ) {
	    count = nuvertexid;
	    break;
	  }
	}
	else{
	  int nparticles = mct.NParticles();
	  hnlvertexid++;
	  for(int i = 0; i < nparticles; i++){
	    simb::MCParticle mcpart = mct.GetParticle(i);
	    
	    if(mcpart.TrackId() == theMCPart.TrackId()){
	      count = hnlvertexid;
	      break;
	    }

	  }
	} // else

      }
    }

    fMCPFromVertexID.push_back(count);
    */
  } // end loop over TrackHandle

  fNTPCTracks = (Int_t)iTrack;
  if(iTrack > 0) fTree->Fill();
  
  return;

} // end HNLAnalysis::FillVectors

//==============================================================================
// Process ionization.  Eventually this moves into the reco code.
void gar::HNLAnalysis::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate, float& forwardIonVal, float& backwardIonVal) {

  // NO CALIBRATION SERVICE FOR NOW

  std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
  forwardIonVal = processOneDirection(SigData, ionizeTruncate);

  SigData = ion.getBAK_dSigdXs();
  backwardIonVal = processOneDirection(SigData, ionizeTruncate);

  return;

}

//==============================================================================
float gar::HNLAnalysis::processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {
  //==============================================================================

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
std::vector< std::pair<int, float> > gar::HNLAnalysis::processPIDInfo(float p){
  //==============================================================================

  std::vector<std::string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};
  std::vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
  std::vector< std::pair<int, float> > pid;
  pid.resize(6);

  int qclosest = 0;
  float dist = 100000000.;
  CLHEP::RandFlat FlatRand(fEngine);

  for(int q = 0; q < 501; ++q){
    //Check the title and the reco momentum take only the one that fits
    std::string fulltitle = m_pidinterp[q]->GetTitle();
    unsigned first = fulltitle.find("=");
    unsigned last = fulltitle.find("GeV");
    std::string substr = fulltitle.substr(first+1, last - first-1);
    float pidinterp_mom = std::atof(substr.c_str());
    //calculate the distance between the bin and mom, store the q the closest
    float disttemp = std::abs(pidinterp_mom - p);
    
    if( disttemp < dist ){
      dist = disttemp;
      qclosest = q;
    }
  } // closes the "pidmatrix" loop
  
  //Compute all the probabities for each type of true to reco
  //loop over the columns (true pid)
  for(int pidm = 0; pidm < 6; ++pidm){
    //loop over the columns (true pid)
    std::vector< std::pair<float, std::string> > v_prob;
    
    //loop over the rows (reco pid)
    for(int pidr = 0; pidr < 6; ++pidr){
      std::string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
      float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1, pidr+1);
      //Need to check random number value and prob value then associate the recopdg to the reco prob
      v_prob.push_back( std::make_pair(prob, recoparticlename) );
    }

    //Compute the pid from it
    if(v_prob.size() > 1){
      //Order the vector of prob
      std::sort(v_prob.begin(), v_prob.end());
      //Throw a random number between 0 and 1
      float random_number = FlatRand.fire();
      //Make cumulative sum to get the range
      std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});

      for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++){
	if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first ){
	  pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ), v_prob.at(ivec+1).first) );
	}
      }
    }
    else{
      pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ), v_prob.at(0).first) );
    }
  }

  //return a vector of pid and prob
  return pid;

}

//==============================================================================
int gar::HNLAnalysis::GetECalPIDParameterized(int pdg, float ptrue){
  //==============================================================================

  CLHEP::RandFlat FlatRand(fEngine);
  float ranNum = FlatRand.fire();

  if(abs(pdg) == 11 || pdg == 22) return pdg;
  else if(abs(pdg) == 2212) return pdg;
  else if(abs(pdg) == 13 || abs(pdg) == 211){
    Int_t pdg2 = 13;
    if(abs(pdg) == 13) pdg2 = 211;

    if(ptrue < 0.48) return pdg;
    if(ptrue >= 0.48 && ptrue < 0.75 && ranNum < 0.8) return pdg;
    if(ptrue >= 0.75 && ptrue < 0.90 && ranNum < 0.9) return pdg;
    if(ptrue >= 0.90 && ranNum < 0.95) return pdg; 

    if(ptrue >= 0.48 && ptrue < 0.75 && ranNum >= 0.8) return pdg2;
    if(ptrue >= 0.75 && ptrue < 0.90 && ranNum >= 0.9) return pdg2;
    if(ptrue >= 0.90 && ranNum >= 0.95) return pdg2;
  }
  else{// for any other case return a random PID
    if(ranNum < 0.25) return 11;
    else if(ranNum >= 0.25 && ranNum < 0.5) return 13;
    else if(ranNum >= 0.5 && ranNum < 0.75) return 211;
    else return 2212;
  }

  return pdg;

}

DEFINE_ART_MODULE(gar::HNLAnalysis)
