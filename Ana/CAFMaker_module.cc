////////////////////////////////////////////////////////////////////////////////
// Class:       CAFMaker
// Plugin Type: analyzer (art v2_11_02)
// File:        CAFMaker_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Additions from Leo Bellantoni, 2019-20
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
#include "SimulationDataProducts/GenieParticle.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/Vee.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "RawDataProducts/CaloRawDigit.h"

#include "SummaryDataProducts/POTSummary.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"

#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TH2.h"
#include "TFile.h"
#include "TVector3.h"

#include "CLHEP/Random/RandFlat.h"

#include <string>
#include <vector>
#include <unordered_map>

#include "StandardRecord/StandardRecord.h"

namespace gar {

  typedef std::pair<float, std::string> P;

  class CAFMaker : public art::EDAnalyzer {
  public:
    explicit CAFMaker(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CAFMaker(CAFMaker const &) = delete;
    CAFMaker(CAFMaker &&) = delete;
    CAFMaker & operator = (CAFMaker const &) = delete;
    CAFMaker & operator = (CAFMaker &&) = delete;

    virtual void beginJob() override;

    // Required functions.
    void analyze(art::Event const & e) override;
    void endRun(const art::Run& run) override;

  private:
    //Working Horse
    cheat::BackTrackerCore* BackTrack;
    void FillGeneratorMonteCarloInfo(art::Event const & e, caf::StandardRecord& rec);
    void FillRawInfo(art::Event const & e, caf::StandardRecord& rec);
    void FillRecoInfo(art::Event const & e, caf::StandardRecord& rec);
    void FillHighLevelRecoInfo(art::Event const & e, caf::StandardRecord& rec);

    //Helpers
    void processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
                               float& forwardIonVal, float& backwardIonVal);
    float processOneDirection(std::vector<std::pair<float,float>> SigData,
                              float ionizeTruncate);
    // Helper method for processOneDirection
    static bool lessThan_byE(std::pair<float,float> a, std::pair<float,float> b)
    {return a.first < b.first;}

    // Compute T for coherent pion analysis
    float computeT( simb::MCTruth theMCTruth );

    // Calculate the track PID based on reco momentum and Tom's parametrization
    std::vector< std::pair<int, float> > processPIDInfo( float p );

    // Position of TPC from geometry service; 1 S Boston Ave.
    float ItsInTulsa[3];
    float xTPC;
    float rTPC;

    // Input data labels
    std::vector<std::string> fGeneratorLabels;
    std::vector<std::string> fGENIEGeneratorLabels;
    std::string fGeantLabel; ///< module label for geant4 simulated hits
    std::string fInstanceLabelCalo; ///< Instance name for ECAL
    std::string fInstanceLabelMuID; ///< Instance name for MuID

    std::string fHitLabel; ///< module label for reco TPC hits rec::Hit
    std::string fTPCClusterLabel; ///< module label for TPC Clusters rec::TPCCluster
    std::string fTrackLabel; ///< module label for TPC Tracks rec:Track
    std::string fTrackTrajectoryLabel; ///< module label for TPC Track Trajectories rec:TrackTrajectory
    std::string fVertexLabel; ///< module label for vertexes rec:Vertex
    std::string fVeeLabel; ///< module label for conversion/decay vertexes rec:Vee

    std::string fRawCaloHitLabel; ///< module label for digitized calo hits raw::CaloRawDigit
    std::string fRawMuIDHitLabel;

    std::string fCaloHitLabel; ///< module label for reco calo hits rec::CaloHit
    std::string fMuIDHitLabel;

    std::string fClusterLabel; ///< module label for calo clusters rec::Cluster
    std::string fClusterMuIDLabel; ///< module label for calo clusters rec::Cluster in MuID
    std::string fPFLabel; ///< module label for reco particles rec::PFParticle
    std::string fECALAssnLabel; ///< module label for track-clusters associations

    std::string fPOTSummaryLabel;

    // Optionally keep/drop parts of the analysis tree
    bool  fWriteMCinfo;        ///< Info from MCTruth, GTruth     Default=true
    bool  fWriteMCPTrajectory; ///< Write MCP Trajectory                Default=true
    bool  fWriteMCPTrajMomenta; ///< Write Momenta associated with MCP Trajectory  Default=false
    bool  fWriteMCCaloInfo;    ///< Write MC info for calorimeter Default=true
    float fMatchMCPtoVertDist; ///< MCParticle to MC vertex match Default=roundoff

    bool  fWriteHits;          ///< Write info about TPC Hits     Default=false
    bool  fWriteTPCClusters;   ///< Write TPCClusters info        Default=true
    bool  fWriteTracks;        ///< Start/end X, P for tracks     Default=true
    bool  fWriteTrackTrajectories;        ///< Point traj of reco tracks     Default=false
    bool  fWriteVertices;      ///< Reco vertexes & their tracks  Default=true
    bool  fWriteVees;          ///< Reco vees & their tracks      Default=true

    bool  fWriteMuID;
    bool  fWriteCaloDigits;    ///< Raw digits for calorimetry.   Default=false
    bool  fWriteCaloHits;      ///< Write ECAL hits.              Default=true
    bool  fWriteCaloClusters;  ///< Write ECAL clusters.          Default=true
    bool  fWriteMatchedTracks; ///< Write ECAL-track Assns        Default=true

    // Truncation parameter for dE/dx (average this fraction of the lowest readings)
    float fIonizTruncate;      ///<                               Default=1.00;
    bool  fWriteCohInfo;       ///< MC level t for coherent pi+.  Default=false

    // the analysis tree
    TTree *fTree;
    TH1* fPOTHist;

    Float_t fTPC_X;      ///< center of TPC stored as per-event & compressed by root
    Float_t fTPC_Y;
    Float_t fTPC_Z;


    //Geometry
    const geo::GeometryCore* fGeo; ///< pointer to the geometry

    typedef int TrkId;
    std::unordered_map<TrkId, Int_t> TrackIdToIndex;

    //map of PID TH2 per momentum value
    CLHEP::HepRandomEngine &fEngine;  ///< random engine
    std::unordered_map<int, TH2F*> m_pidinterp;
  };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::CAFMaker::CAFMaker(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed")) {
  fGeo     = gar::providerFrom<geo::GeometryGAr>();

  bool usegenlabels =
    p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
  if (!usegenlabels) fGeneratorLabels.clear();

  bool usegeniegenlabels  =
    p.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
  if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

  //Sim Hits
  fGeantLabel        = p.get<std::string>("GEANTLabel","geant");
  fInstanceLabelCalo = p.get<std::string>("InstanceLabelCalo","ECAL");
  fInstanceLabelMuID = p.get<std::string>("InstanceLabelMuID","MuID");

  fHitLabel         = p.get<std::string>("HitLabel","hit");
  fTPCClusterLabel  = p.get<std::string>("TPCClusterLabel","tpccluster");
  fTrackLabel       = p.get<std::string>("TrackLabel","track");
  fTrackTrajectoryLabel  = p.get<std::string>("TrackTrajectoryLabel","track");
  fVertexLabel      = p.get<std::string>("VertexLabel","vertex");
  fVeeLabel         = p.get<std::string>("VeeLabel","veefinder1");

  //Calorimetric related ECAL/MuID
  fRawCaloHitLabel     = p.get<std::string>("RawCaloHitLabel","daqsipm");
  fRawMuIDHitLabel     = p.get<std::string>("RawMuIDHitLabel","daqsipmmuid");
  fCaloHitLabel        = p.get<std::string>("CaloHitLabel","sipmhit");
  fMuIDHitLabel        = p.get<std::string>("MuIDHitLabel","sipmhit");

  fClusterLabel     = p.get<std::string>("ClusterLabel","calocluster");
  fClusterMuIDLabel = p.get<std::string>("MuIDClusterLabel","caloclustermuid");
  fPFLabel          = p.get<std::string>("PFLabel","pandora");
  fECALAssnLabel    = p.get<std::string>("ECALAssnLabel","trkecalassn");

  fPOTSummaryLabel  = p.get<std::string>("POTSummaryLabel");

  // What to write
  fWriteMCinfo              = p.get<bool>("WriteMCinfo",        true);
  fWriteMCPTrajectory       = p.get<bool>("WriteMCPTrajectory", true);
  fWriteMCPTrajMomenta      = p.get<bool>("WriteMCPTrajMomenta",false);
  fWriteMCCaloInfo          = p.get<bool>("WriteMCCaloInfo",    true);
  float MCPtoVertDefault    = 10.0*std::numeric_limits<Float_t>::epsilon();
  fMatchMCPtoVertDist       = p.get<float>("MatchMCPtoVertDist",MCPtoVertDefault);

  fWriteHits                = p.get<bool>("WriteHits",         false);
  fWriteTPCClusters         = p.get<bool>("WriteTPCClusters",  true);
  fWriteTracks              = p.get<bool>("WriteTracks",       true);
  fWriteTrackTrajectories   = p.get<bool>("WriteTrackTrajectories", false);
  fWriteVertices            = p.get<bool>("WriteVertices",     true);
  fWriteVees                = p.get<bool>("WriteVees",         true);

  fWriteMuID                = p.get<bool>("WriteMuID",         true);
  fWriteCaloDigits          = p.get<bool>("WriteCaloDigits",   false);
  fWriteCaloHits            = p.get<bool>("WriteCaloHits",     true);
  fWriteCaloClusters        = p.get<bool>("WriteCaloClusters", true);
  fWriteMatchedTracks       = p.get<bool>("WriteMatchedTracks",true);

  fIonizTruncate            = p.get<float>("IonizTruncate",    0.70);
  fWriteCohInfo             = p.get<bool> ("WriteCohInfo",     false);

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

  consumesMany<std::vector<sdp::GenieParticle> >();
  //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<simb::MCParticle> >(fGeantLabel);

  //TPC related
  consumes<std::vector<sdp::EnergyDeposit> >(fGeantLabel);
  consumes<std::vector<rec::TPCCluster> >(fTPCClusterLabel);
  consumes<art::Assns<rec::Track, rec::TPCCluster> >(fTPCClusterLabel);
  consumes<std::vector<rec::Hit> >(fHitLabel);
  consumes<std::vector<rec::Track> >(fTrackLabel);
  consumes<std::vector<rec::TrackTrajectory> >(fTrackTrajectoryLabel);
  consumes<std::vector<rec::Vertex> >(fVertexLabel);
  consumes<art::Assns<rec::Track, rec::Vertex> >(fVertexLabel);
  consumes<std::vector<rec::Vee> >(fVeeLabel);
  consumes<art::Assns<rec::Track, rec::Vee> >(fVeeLabel);
  consumes<std::vector<rec::TrackIoniz>>(fTrackLabel);
  consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);

  //Calorimetry related
  art::InputTag ecalgeanttag(fGeantLabel, fInstanceLabelCalo);
  consumes<std::vector<gar::sdp::CaloDeposit> >(ecalgeanttag);

  art::InputTag ecalrawtag(fRawCaloHitLabel, fInstanceLabelCalo);
  consumes<std::vector<raw::CaloRawDigit> >(ecalrawtag);

  art::InputTag ecalhittag(fCaloHitLabel, fInstanceLabelCalo);
  consumes<std::vector<rec::CaloHit> >(ecalhittag);

  art::InputTag ecalclustertag(fClusterLabel, fInstanceLabelCalo);
  consumes<std::vector<rec::Cluster> >(ecalclustertag);
  consumes<art::Assns<rec::Cluster, rec::CaloHit>>(ecalclustertag);

  //Muon system related
  if (fGeo->HasMuonDetector() && fWriteMuID) {
    art::InputTag muidgeanttag(fGeantLabel, fInstanceLabelMuID);
    consumes<std::vector<gar::sdp::CaloDeposit> >(muidgeanttag);

    art::InputTag muidrawtag(fRawMuIDHitLabel, fInstanceLabelMuID);
    consumes<std::vector<raw::CaloRawDigit> >(muidrawtag);

    art::InputTag muidhittag(fCaloHitLabel, fInstanceLabelMuID);
    consumes<std::vector<rec::CaloHit> >(muidhittag);

    art::InputTag muidclustertag(fClusterMuIDLabel, fInstanceLabelMuID);
    consumes<std::vector<rec::Cluster> >(muidclustertag);
    consumes<art::Assns<rec::Cluster, rec::CaloHit>>(muidclustertag);
  }

  consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

  return;
} // end constructor



//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::beginJob() {

  fTPC_X = ItsInTulsa[0] = fGeo->TPCXCent();
  fTPC_Y = ItsInTulsa[1] = fGeo->TPCYCent();
  fTPC_Z = ItsInTulsa[2] = fGeo->TPCZCent();

  xTPC = fGeo->TPCLength() / 2.;
  rTPC = fGeo->TPCRadius();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("recTree", "recTree");

  // Create the branch. We will update the address before we write the tree
  caf::StandardRecord* rec = 0;
  fTree->Branch("rec", &rec);

  fPOTHist = fPOTSummaryLabel.empty() ? 0 : tfs->make<TH1F>("TotalPOT", "", 1, 0, 1);

  if (fWriteMCPTrajectory) {
    // Checking for parameter file validity only once to simplify
    // later code in CAFMaker::process etc.
    if (!fWriteMCinfo) {
      throw cet::exception("CAFMaker")
        << " fWriteMCPTrajectory, but !fWriteMCinfo."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    } 
  }

  if (fWriteMCCaloInfo) {
    if (!fWriteMCinfo) {
      throw cet::exception("CAFMaker")
        << " fWriteMCCaloInfo, but !fWriteMCinfo."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
  }

  // Reco'd verts & their track-ends
  if (fWriteVertices) {
    if (!fWriteTracks) {
      throw cet::exception("CAFMaker")
        << " fWriteVertices, but !fWriteTracks."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
  }

  // Reco'd vees & their track-ends
  if (fWriteVees) {
    if (!fWriteTracks) {
      throw cet::exception("CAFMaker")
        << " fWriteVees, but !fWriteTracks."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
  }

  if (fWriteMatchedTracks) {
    if (!fWriteTracks || !fWriteCaloClusters) {
      throw cet::exception("CAFMaker")
        << " fWriteMatchedTracks, but (!fWriteTracks || !fWriteCaloClusters)."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
  }

  std::string filename = "${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
  TFile infile(filename.c_str(), "READ");

  m_pidinterp.clear();
  char str[11];
  for (int q = 0; q < 501; ++q) {
    sprintf(str, "%d", q);
    std::string s = "pidmatrix";
    s.append(str);
    // read the 500 histograms one by one; each histogram is a
    // 6 by 6 matrix of probabilities for a given momentum value
    m_pidinterp.insert( std::make_pair(q, (TH2F*) infile.Get(s.c_str())->Clone("pidinterp")) );
  }

  return;
}  // End of :CAFMaker::beginJob



//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::analyze(art::Event const & e) {

  caf::StandardRecord rec;
  caf::StandardRecord* prec = &rec;
  fTree->SetBranchAddress("rec", &prec);

  rec.hdr.run    = e.run();
  rec.hdr.subrun = e.subRun();
  rec.hdr.event  = e.id().event();

  rec.hdr.TPC_X = fTPC_X;
  rec.hdr.TPC_Y = fTPC_Y;
  rec.hdr.TPC_Z = fTPC_Z;


  // Need a non-constant backtracker instance, for now, in anayze ot beginJob
  cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
  BackTrack = const_cast<cheat::BackTrackerCore*>(const_bt);

  //Fill generator and MC Information
  if (fWriteMCinfo) FillGeneratorMonteCarloInfo(e, rec);

  //Fill Raw Information
  if (fWriteCaloDigits) FillRawInfo(e, rec);

  //Fill Reco Information
  FillRecoInfo(e, rec);

  //Fill HL Reco Information
  FillHighLevelRecoInfo(e, rec);

  fTree->Fill();
  return;
}

//==============================================================================
void gar::CAFMaker::endRun(const art::Run& run)
{
  if(fPOTHist){
    art::InputTag itag(fPOTSummaryLabel);
    auto potsum = run.getHandle<gar::sumdata::POTSummary>(itag);
    fPOTHist->Fill(.5, potsum->totgoodpot);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::FillGeneratorMonteCarloInfo(art::Event const & e, caf::StandardRecord& rec) {

  // =============  Get art handles ==========================================
  // Get handles for MCinfo, also good for MCPTrajectory
  std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
  std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;

  if (fGeneratorLabels.size()<1) {
    mcthandlelist = e.getMany<std::vector<simb::MCTruth> >();    // get them all (even if there are none)
  } else {
    mcthandlelist.resize(fGeneratorLabels.size());
    for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
      // complain if we wanted a specific one but didn't find it
      art::InputTag itag(fGeneratorLabels.at(i));
      mcthandlelist.at(i) = e.getHandle< std::vector<simb::MCTruth> >(itag);
      if (!mcthandlelist.at(i)) {
        throw cet::exception("CAFMaker") << " No simb::MCTruth branch."
                                         << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }
  }

  if (fGENIEGeneratorLabels.size()<1) {
    gthandlelist = e.getMany<std::vector<simb::GTruth> >();  // get them all (even if there are none)
  } else {
    gthandlelist.resize(fGENIEGeneratorLabels.size());
    for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i) {
      gthandlelist.at(i) = e.getHandle<std::vector<simb::GTruth> >(fGENIEGeneratorLabels.at(i));
      // complain if we wanted a specific one but didn't find it
      if (!gthandlelist.at(i)) {
        throw cet::exception("CAFMaker") << " No simb::GTruth branch."
                                         << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }
  }

  // save MCTruth info
  for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
    for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
      if (mct.NeutrinoSet()) {
        simb::MCNeutrino nuw = mct.GetNeutrino();
        caf::SRNeutrino srnu;
        srnu.pdg = nuw.Nu().PdgCode();
        srnu.ccnc = nuw.CCNC();
        srnu.mode = nuw.Mode();
        srnu.interactionType = nuw.InteractionType();
        srnu.Q2 = nuw.QSqr();
        srnu.W = nuw.W();
        srnu.X = nuw.X();
        srnu.Y = nuw.Y();
        srnu.theta = nuw.Theta();
        if (fWriteCohInfo) {
          srnu.T = computeT(mct);
        }
        srnu.vtx = caf::SRVector3D(nuw.Nu().EndX(), nuw.Nu().EndY(), nuw.Nu().EndZ());
        srnu.p = caf::SRVector3D(nuw.Nu().Px(), nuw.Nu().Py(), nuw.Nu().Pz());

        rec.mc.nu.push_back(srnu);
      }  // end MC info from MCTruth
    }
  }

  rec.mc.nnu = rec.mc.nu.size();

  unsigned int nuidx = 0;
  // save GTruth info
  for (size_t igthl = 0; igthl < gthandlelist.size(); ++igthl) {
    for ( auto const& gt : (*gthandlelist.at(igthl)) ) {
      if(nuidx >= rec.mc.nu.size()) abort();

      // We're assuming we visit the gtruths in the same order (should
      // use an Assn?)
      rec.mc.nu[nuidx].gint = gt.fGint;
      rec.mc.nu[nuidx].tgtPDG = gt.ftgtPDG;
      rec.mc.nu[nuidx].weight = gt.fweight;
      rec.mc.nu[nuidx].gT = gt.fgT;

      ++nuidx;
    }
  }

  // Save the particle list from the GENIE event record

  auto gparthandlelist = e.getMany<std::vector<sdp::GenieParticle>>();

  for (size_t igphl = 0; igphl < gparthandlelist.size(); ++igphl) {
    for ( auto const& gpart : (*gparthandlelist.at(igphl)) ) {
      caf::SRGenieParticle srgpart;
      srgpart.intidx = gpart.InteractionIndex();
      srgpart.idx = gpart.Index();
      srgpart.pdg = gpart.Pdg();
      srgpart.status = gpart.Status();
      srgpart.name = gpart.Name();
      srgpart.firstmom = gpart.FirstMother();
      srgpart.lastmom = gpart.LastMother();
      srgpart.firstdaugh = gpart.FirstDaughter();
      srgpart.lastdaugh = gpart.LastDaughter();
      srgpart.p = caf::SRVector3D(gpart.Px(), gpart.Py(), gpart.Pz());
      srgpart.E = gpart.E();
      srgpart.mass = gpart.Mass();

      rec.mc.gpart.push_back(srgpart);
    }
  }

  rec.mc.ngpart = rec.mc.gpart.size();

  auto MCPHandle = e.getHandle< std::vector<simb::MCParticle> >(fGeantLabel);
  if (!MCPHandle) {
    throw cet::exception("CAFMaker") << " No simb::MCParticle branch."
                                     << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  // Save MCParticle info (post-GEANT; for pre-GEANT, maybe try MCTruth?)
  // Helps to have a map from the TrackId from generator to index number
  // into *MCPHandle, which is a vector<MCParticle>.  Not all TrackId value
  // appear in MCPHandle, although I think the missing ones are all from
  // the GEANT rather than the GENIE stage of GENIEGen.  (*MCPHandle).size()
  // is the same as eg fMCPStartX.size() by the time this entry is written.
  // TracIdToIndex isa class variable for accessibility by per-Track and
  // per-Cluster code
  Int_t index = 0;
  for ( auto const& mcp : (*MCPHandle) ) {
    int TrackId = mcp.TrackId();
    TrackIdToIndex[TrackId] = index++;
  }

  for ( auto const& mcp : (*MCPHandle) ) {
    caf::SRParticle srpart;
    srpart.trkid = mcp.TrackId();
    srpart.pdg = mcp.PdgCode();

    // If mcp.Mother() == 0, particle is from initial vertex;
    // Set MCMotherIndex to -1 in that case.
    TrkId momTrkId = mcp.Mother();
    Int_t momIndex = -1;
    int momPDG = 0;
    if (momTrkId>0) {
      //Check if it exists!
      if(TrackIdToIndex.find(momTrkId) != TrackIdToIndex.end()){
        momIndex = TrackIdToIndex[momTrkId];
        momPDG   = (*MCPHandle).at(momIndex).PdgCode();
      } else {
        MF_LOG_DEBUG("CAFMaker_module")
          << " mcp trkid " << mcp.TrackId()
          << " pdg code " << mcp.PdgCode()
          << " could not find mother trk id " << momTrkId
          << " in the TrackIdToIndex map"
          << " creating process is [ " << mcp.Process() << " ]";
      }
    }

    srpart.momidx = momIndex;
    srpart.momtrkid = mcp.Mother(); //directly trackid not index
    srpart.mompdg = momPDG;

    const TLorentzVector& position = mcp.Position(0);
    const TLorentzVector& momentum = mcp.Momentum(0);
    srpart.start.x = position.X();
    srpart.start.y = position.Y();
    srpart.start.z = position.Z();
    srpart.time = mcp.T();
    srpart.startp.x = momentum.Px();
    srpart.startp.y = momentum.Py();
    srpart.startp.z = momentum.Pz();

    const TLorentzVector& positionEnd = mcp.EndPosition();
    const TLorentzVector& momentumEnd = mcp.EndMomentum();
    srpart.end.x = positionEnd.X();
    srpart.end.y = positionEnd.Y();
    srpart.end.z = positionEnd.Z();
    srpart.endp.x = momentumEnd.Px();
    srpart.endp.y = momentumEnd.Py();
    srpart.endp.z = momentumEnd.Pz();
    srpart.proc = mcp.Process();
    srpart.endproc = mcp.EndProcess();

    rec.mc.part.push_back(srpart);
  }

  rec.mc.npart = rec.mc.part.size();

  // Get vertex for each MCParticle.  In principle, the mct.GetParticle()
  // method should produce all of them; filter out the ones with a non-stable
  // status code and the GENIE pseudo-particles (Pid>= 2000000000) and you
  // should have the contents of (*MCPHandle).  But you don't.  On the 1st
  // event I (Leo Bellantoni) tested, with 40 vertices and 351 primary
  // particles, particle 143, a proton, evidently decayed or something and is
  // not in (*MCPHandle).  So we use the following primitive earth technology.
  size_t nMCParticles = (*MCPHandle).size();
  size_t iMCParticle  = 0;
  for (; iMCParticle<nMCParticles; ++iMCParticle) {
    // Assign noprimary to start with
    rec.mc.part[iMCParticle].vtxidx = -1;
    // Do the primaries first
    if (rec.mc.part[iMCParticle].momidx!=-1) break;
    Float_t trackX = rec.mc.part[iMCParticle].start.x;
    Float_t trackY = rec.mc.part[iMCParticle].start.y;
    Float_t trackZ = rec.mc.part[iMCParticle].start.z;
    int vertexIndex = 0;
    for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
      for ( auto const& mct : (*mcthandlelist.at(imchl)) ) {
        if (mct.NeutrinoSet()) {
          simb::MCNeutrino nuw = mct.GetNeutrino();
          Float_t vertX = nuw.Nu().EndX();
          Float_t vertY = nuw.Nu().EndY();
          Float_t vertZ = nuw.Nu().EndZ();
          Float_t dist = std::hypot(trackX-vertX,trackY-vertY,trackZ-vertZ);
          if ( dist <= fMatchMCPtoVertDist ) {
            rec.mc.part[iMCParticle].vtxidx = vertexIndex;
            goto foundMCvert; // break out of all inner loops
          }
        }
        ++vertexIndex;
      }
    }
  foundMCvert:
    ;
  }

  // Now the secondaries.  As they are after the primaries, do not re-init iMCParticle
  for (; iMCParticle<nMCParticles; ++iMCParticle) {
    int momIndex = rec.mc.part[iMCParticle].momidx;
    int lastMCParticle = iMCParticle;
    while (momIndex != -1) {
      lastMCParticle = momIndex;
      momIndex       = rec.mc.part[momIndex].momidx;
    }
    rec.mc.part[iMCParticle].vtxidx = rec.mc.part[lastMCParticle].vtxidx;
  }

  if (fWriteMCPTrajectory) {
    // It's in the MCParticle table
    Int_t mcpIndex = 0;
    for ( auto const& mcp : (*MCPHandle) ) {
      const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
      const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
      //No charge don't store the trajectory
      if (definition==nullptr || definition->Charge() == 0) continue;
      //TrackID of the mcp to keep track to which mcp this trajectory is
      int trackId = mcp.TrackId();
      for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++) {
        float xTraj = mcp.Trajectory().X(iTraj);
        float yTraj = mcp.Trajectory().Y(iTraj);
        float zTraj = mcp.Trajectory().Z(iTraj);
        float rTraj = std::hypot( yTraj - ItsInTulsa[1], zTraj - ItsInTulsa[2]);

        if (abs(xTraj - ItsInTulsa[0]) > xTPC) continue;
        if (rTraj > rTPC) continue;

        caf::SRParticleTrajectoryPoint srtrajpt;
        srtrajpt.x = xTraj;
        srtrajpt.y = yTraj;
        srtrajpt.z = zTraj;
        srtrajpt.t = mcp.Trajectory().T(iTraj);
        if (fWriteMCPTrajMomenta) {
          srtrajpt.px = mcp.Trajectory().Px(iTraj);
          srtrajpt.py = mcp.Trajectory().Py(iTraj);
          srtrajpt.pz = mcp.Trajectory().Pz(iTraj);
        }
        srtrajpt.E = mcp.Trajectory().E(iTraj);

        rec.mc.part[mcpIndex].traj.push_back(srtrajpt);
      }
      rec.mc.part[mcpIndex].trkid = trackId;

      mcpIndex++;
    }
  }

  // Get handles for MCCaloInfo
  art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;//ecal
  art::Handle< std::vector<gar::sdp::CaloDeposit> > MuIDSimHitHandle;//ecal
  art::InputTag ecalgeanttag(fGeantLabel, fInstanceLabelCalo);
  art::InputTag muidgeanttag(fGeantLabel, fInstanceLabelMuID);
  if (fWriteMCCaloInfo) {
    SimHitHandle = e.getHandle< std::vector<gar::sdp::CaloDeposit> >(ecalgeanttag);
    if (!SimHitHandle) {
      throw cet::exception("CAFMaker") << " No gar::sdp::CaloDeposit branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    if(fWriteMuID) {
      if (fGeo->HasMuonDetector())
        {
          MuIDSimHitHandle = e.getHandle< std::vector<gar::sdp::CaloDeposit> >(muidgeanttag);
          if (!MuIDSimHitHandle) {
            throw cet::exception("CAFMaker") << " No gar::sdp::CaloDeposit branch."
                                             << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
          }
        }
    }

    // Save simulation ecal hit info
    for ( auto const& SimHit : (*SimHitHandle) ) {
      caf::SRSimHit srsimhit;
      srsimhit.x = SimHit.X();
      srsimhit.y = SimHit.Y();
      srsimhit.z = SimHit.Z();
      srsimhit.t = SimHit.Time();
      srsimhit.E = SimHit.Energy();
      srsimhit.trkid = SimHit.TrackID();
      srsimhit.cellid = SimHit.CellID();

      rec.mc.simhit.ecal.push_back(srsimhit);

      rec.mc.simhit.ecaltotE += SimHit.Energy();
    }

    rec.mc.simhit.necal = rec.mc.simhit.ecal.size();

    // Save simulation muon system hit info
    for ( auto const& SimHit : (*MuIDSimHitHandle) ) {
      caf::SRSimHit srsimhit;
      srsimhit.x = SimHit.X();
      srsimhit.y = SimHit.Y();
      srsimhit.z = SimHit.Z();
      srsimhit.t = SimHit.Time();
      srsimhit.E = SimHit.Energy();
      srsimhit.trkid = SimHit.TrackID();
      srsimhit.cellid = SimHit.CellID();

      rec.mc.simhit.muid.push_back(srsimhit);

      rec.mc.simhit.muidtotE += SimHit.Energy();
    }

    rec.mc.simhit.nmuid = rec.mc.simhit.muid.size();
  }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::FillRawInfo(art::Event const & e, caf::StandardRecord& rec) {


  // Get handle for CaloDigits
  art::InputTag ecalrawtag(fRawCaloHitLabel, fInstanceLabelCalo);
  auto RawHitHandle = e.getHandle<std::vector<gar::raw::CaloRawDigit> >(ecalrawtag);
  if (!RawHitHandle) {
    throw cet::exception("CAFMaker") << " No :raw::CaloRawDigit branch."
                                     << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  art::InputTag muidrawtag(fRawMuIDHitLabel, fInstanceLabelMuID);
  art::Handle< std::vector<gar::raw::CaloRawDigit> > MuIDRawHitHandle;
  if (fGeo->HasMuonDetector())
    {
      MuIDRawHitHandle = e.getHandle< std::vector<gar::raw::CaloRawDigit> >(muidrawtag);
      if (!MuIDRawHitHandle) {
        throw cet::exception("CAFMaker") << " No :raw::CaloRawDigit branch."
                                         << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }

  // save ecal raw digits info
  for ( auto const& DigiHit : (*RawHitHandle) ) {
    caf::SRDigit srdig;
    srdig.x = DigiHit.X();
    srdig.y = DigiHit.Y();
    srdig.z = DigiHit.Z();
    srdig.t = (DigiHit.Time().first + DigiHit.Time().second) / 2.0 ;
    srdig.adc = DigiHit.ADC().first;
    srdig.cellid = DigiHit.CellID();

    rec.dig.ecal.push_back(srdig);
  }
  rec.dig.necal = rec.dig.ecal.size();

  // save muon system raw digits info
  for ( auto const& DigiHit : (*MuIDRawHitHandle) ) {
    caf::SRDigit srdig;
    srdig.x = DigiHit.X();
    srdig.y = DigiHit.Y();
    srdig.z = DigiHit.Z();
    srdig.t = (DigiHit.Time().first + DigiHit.Time().second) / 2.0 ;
    srdig.adc = DigiHit.ADC().first;
    srdig.cellid = DigiHit.CellID();
  }
  rec.dig.nmuid = rec.dig.muid.size();
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::FillRecoInfo(art::Event const & e, caf::StandardRecord& rec) {

  // Get handle for TPC hit data
  if (fWriteHits) {
    auto HitHandle = e.getHandle< std::vector<rec::Hit> >(fHitLabel);
    if (!HitHandle) {
      throw cet::exception("CAFMaker") << " No rec::Hit branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // save hits in the TPC
    for ( auto const& Hit : (*HitHandle) ) {
      caf::SRTPCRecoHit srhit;
      srhit.x = Hit.Position()[0];
      srhit.y = Hit.Position()[1];
      srhit.z = Hit.Position()[2];
      srhit.sig = Hit.Signal();
      srhit.rms = Hit.RMS();
      srhit.chan = Hit.Channel();
      rec.hit.tpc.push_back(srhit);
    }
    rec.hit.ntpc = rec.hit.tpc.size();
  }

  // Get handle for CaloHits

  if (fWriteCaloHits) {

    art::InputTag ecalrecotag(fCaloHitLabel, fInstanceLabelCalo);
    auto RecoHitHandle = e.getHandle< std::vector<rec::CaloHit> >(ecalrecotag);
    if (!RecoHitHandle) {
      throw cet::exception("CAFMaker") << " No rec::CaloHit branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // save reco'd Calorimetry hits
    for ( auto const& Hit : (*RecoHitHandle) ) {
      caf::SRRecoHit srhit;
      srhit.id = Hit.getIDNumber();
      srhit.x = Hit.Position()[0];
      srhit.y = Hit.Position()[1];
      srhit.z = Hit.Position()[2];
      srhit.t = Hit.Time().first;
      srhit.E = Hit.Energy();
      srhit.cellid = Hit.CellID();
      rec.hit.ecal.push_back(srhit);

      rec.hit.ecaltotE += srhit.E;
    }
    rec.hit.necal = rec.hit.ecal.size();


    if(fWriteMuID) {
      art::InputTag muirecotag(fMuIDHitLabel, fInstanceLabelMuID);
      art::Handle< std::vector<rec::CaloHit> > MuIDRecoHitHandle;
      if (fGeo->HasMuonDetector())
        { 
          MuIDRecoHitHandle = e.getHandle< std::vector<rec::CaloHit> >(muirecotag);       
          if(!MuIDRecoHitHandle) {
            throw cet::exception("CAFMaker") << " No rec::CaloHit branch."
                                             << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
          }
        }
      for ( auto const& Hit : (*MuIDRecoHitHandle) ) {
        caf::SRRecoHit srhit;
        srhit.id = Hit.getIDNumber();
        srhit.x = Hit.Position()[0];
        srhit.y = Hit.Position()[1];
        srhit.z = Hit.Position()[2];
        srhit.t = Hit.Time().first;
        srhit.E = Hit.Energy();
        srhit.cellid = Hit.CellID();
        rec.hit.muid.push_back(srhit);

        rec.hit.muidtotE += srhit.E;
      }
      rec.hit.nmuid = rec.hit.muid.size();

    }

  }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::CAFMaker::FillHighLevelRecoInfo(art::Event const & e, caf::StandardRecord& rec) {

  // Get handle for TPCClusters
  art::Handle< std::vector<rec::TPCCluster> > TPCClusterHandle;
  if (fWriteTPCClusters) {
    TPCClusterHandle = e.getHandle< std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    if (!TPCClusterHandle) {
      throw cet::exception("CAFMaker") << " No rec::TPCCluster branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
  }

  // Get handles for Tracks and their ionizations; also Assn's to TPCClusters, TrackIoniz
  art::Handle< std::vector<rec::Track> > TrackHandle;
  art::Handle< std::vector<rec::TrackIoniz> > TrackIonHandle;
  art::Handle< std::vector<rec::TrackTrajectory> > TrackTrajHandle;
  art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
  art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
  if(fWriteTracks) {
    TrackHandle = e.getHandle< std::vector<rec::Track> >(fTrackLabel);
    if (!TrackHandle) {
      throw cet::exception("CAFMaker") << " No rec::Track branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    TrackIonHandle = e.getHandle< std::vector<rec::TrackIoniz> >(fTrackLabel);
    if (!TrackIonHandle) {
      throw cet::exception("CAFMaker") << " No rec::TrackIoniz branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);
    findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,e,fTrackLabel);

    if(fWriteTrackTrajectories) {
      TrackTrajHandle = e.getHandle< std::vector<rec::TrackTrajectory> >(fTrackTrajectoryLabel);
      if (!TrackTrajHandle) {
        throw cet::exception("CAFMaker") << " No rec::TrackTrajectory branch."
                                         << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }
  }

  // Get handle for Vertices; also Assn's to Tracks
  art::Handle< std::vector<rec::Vertex> > VertexHandle;
  art::FindManyP<rec::Track, rec::TrackEnd>* findManyTrackEnd = NULL;
  if (fWriteVertices) {
    VertexHandle = e.getHandle< std::vector<rec::Vertex> >(fVertexLabel);
    if (!VertexHandle) {
      throw cet::exception("CAFMaker") << " No rec::Vertex branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    findManyTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VertexHandle,e,fVertexLabel);
  }

  // Get handle for Vees; also Assn's to Tracks
  art::Handle< std::vector<rec::Vee> > VeeHandle;
  art::FindManyP<rec::Track, rec::TrackEnd>* findManyVeeTrackEnd = NULL;
  if (fWriteVees) {
    VeeHandle = e.getHandle< std::vector<rec::Vee> >(fVeeLabel);
    if (!VeeHandle) {
      throw cet::exception("CAFMaker") << " No rec::Vee branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    findManyVeeTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(VeeHandle,e,fVeeLabel);
  }

  // Get handle for CaloClusters; also Assn for matching tracks
  art::InputTag ecalclustertag(fClusterLabel, fInstanceLabelCalo);
  art::InputTag muidclustertag(fClusterMuIDLabel, fInstanceLabelMuID);
  art::Handle< std::vector<rec::Cluster> > RecoClusterHandle;
  art::Handle< std::vector<rec::Cluster> > RecoClusterMuIDHandle;
  art::FindManyP<rec::Track, rec::TrackEnd>* findManyCALTrackEnd = NULL;
  // art::FindManyP<gar::rec::CaloHit>* findManyClusterRecoHit = NULL;
  // art::FindManyP<gar::rec::CaloHit>* findManyClusterMuIDHit = NULL;

  if (fWriteCaloClusters) {
    RecoClusterHandle = e.getHandle< std::vector<rec::Cluster> >(ecalclustertag);
    if (!RecoClusterHandle) {
      throw cet::exception("CAFMaker") << " No rec::Cluster branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if(fWriteMuID) {
      if (fGeo->HasMuonDetector())
	{
	  RecoClusterMuIDHandle = e.getHandle< std::vector<rec::Cluster> >(muidclustertag);
	  if (!RecoClusterMuIDHandle) {
	    throw cet::exception("CAFMaker") << " No rec::Cluster (MuID) branch."
					     << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
	  }
	}
    }

    // findManyClusterRecoHit = new art::FindManyP<gar::rec::CaloHit>(RecoClusterHandle,e,ecalclustertag);

    // if (fGeo->HasMuonDetector() && fWriteMuID)
    //     findManyClusterMuIDHit = new art::FindManyP<gar::rec::CaloHit>(RecoClusterMuIDHandle,e,muidclustertag);

    if (fWriteTracks)
      findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(RecoClusterHandle,e,fECALAssnLabel);
  }

  // save clusters in the TPC. For some reason, can't get FindOneP<rec::Track> or
  // FindManyP<rec::Track> to work; seems the underlying Assn isn't found.  Have
  // to FindManyP<TPCCluster> instead and  iterate if (fWriteTracks).  :(
  if (fWriteTPCClusters) {
    for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
      caf::SRTPCCluster srclust;
      srclust.x = TPCCluster.Position()[0];
      srclust.y = TPCCluster.Position()[1];
      srclust.z = TPCCluster.Position()[2];
      srclust.sig = TPCCluster.Signal();
      srclust.rms = TPCCluster.RMS();
      const float* cov = TPCCluster.CovMatPacked();
      srclust.cov.xx = cov[0];
      srclust.cov.xy = cov[1];
      srclust.cov.xz = cov[2];
      srclust.cov.yy = cov[3];
      srclust.cov.yz = cov[4];
      srclust.cov.zz = cov[5];

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
      srclust.trkid = trackForThisTPCluster;

      rec.clust.tpc.push_back(srclust);
    }
    rec.clust.ntpc = rec.clust.tpc.size();
  }

  // save per-track info
  if (fWriteTracks) {
    size_t iTrack = 0;
    for ( auto const& track : (*TrackHandle) ) {
      caf::SRTrack srtrk;

      // track is a rec::Track, not a rec::TrackPar
      srtrk.id = track.getIDNumber();

      srtrk.start = caf::SRVector3D(track.Vertex()[0],
                                    track.Vertex()[1],
                                    track.Vertex()[2]);
      srtrk.startp = caf::SRVector3D(track.Momentum_beg()*track.VtxDir()[0],
                                     track.Momentum_beg()*track.VtxDir()[1],
                                     track.Momentum_beg()*track.VtxDir()[2]);
      srtrk.startq = track.ChargeBeg();

      srtrk.end = caf::SRVector3D(track.End()[0],
                                  track.End()[1],
                                  track.End()[2]);
      srtrk.endp = caf::SRVector3D(track.Momentum_end()*track.EndDir()[0],
                                   track.Momentum_end()*track.EndDir()[1],
                                   track.Momentum_end()*track.EndDir()[2]);
      srtrk.endq = track.ChargeEnd();

      srtrk.fwd.len = track.LengthForward();
      srtrk.bak.len = track.LengthBackward();
      srtrk.fwd.chi2 = track.ChisqForward();
      srtrk.bak.chi2 = track.ChisqBackward();
      srtrk.NTPCClustersOnTrack = track.NHits();

      //Add the PID information based on Tom's parametrization
      TVector3 momF(track.Momentum_beg()*track.VtxDir()[0], track.Momentum_beg()*track.VtxDir()[1], track.Momentum_beg()*track.VtxDir()[2]);
      TVector3 momB(track.Momentum_end()*track.EndDir()[0], track.Momentum_end()*track.EndDir()[1], track.Momentum_end()*track.EndDir()[2]);

      //Reconstructed momentum forward and backward
      float pF = momF.Mag();
      float pB = momB.Mag();
      std::vector< std::pair<int, float> > pidF = processPIDInfo( pF );
      std::vector< std::pair<int, float> > pidB = processPIDInfo( pB );

      //Fill the pid and its probability
      for(size_t ipid = 0; ipid < pidF.size(); ipid++) {
        caf::SRTrackPID srpid;
        srpid.id = pidF.at(ipid).first;
        srpid.prob = pidF.at(ipid).second;
        srtrk.fwd.pid.push_back(srpid);
      }
      for(size_t ipid = 0; ipid < pidB.size(); ipid++) {
        caf::SRTrackPID srpid;
        srpid.id = pidB.at(ipid).first;
        srpid.prob = pidB.at(ipid).second;
        srtrk.bak.pid.push_back(srpid);
      }

      // Matching MCParticle info
      std::vector<std::pair<simb::MCParticle*,float>> trakt;
      trakt = BackTrack->TrackToMCParticles( const_cast<rec::Track*>(&track) );
      srtrk.mcidx = -1;
      if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
        srtrk.mcidx = TrackIdToIndex[trakt[0].first->TrackId()];
      }

      if (srtrk.mcidx > -1) {
        srtrk.mcfrac = trakt[0].second;
      } else {
        srtrk.mcfrac = 0;
      }

      if (findIonization->isValid()) {
        // No calibration for now.  Someday this should all be in reco
        rec::TrackIoniz ionization = *(findIonization->at(iTrack));
        float avgIonF, avgIonB;
        processIonizationInfo(ionization, fIonizTruncate, avgIonF, avgIonB);
        srtrk.fwd.avgion = avgIonF;
        srtrk.bak.avgion = avgIonB;
      } else {
        srtrk.fwd.avgion = 0;
        srtrk.bak.avgion = 0;
      }
      iTrack++;

      rec.trk.push_back(srtrk);
    } // end loop over TrackHandle
  }       

  rec.ntrk = rec.trk.size();

  //TrackTrajectories
  if(fWriteTrackTrajectories) {
    size_t iTrackTraj = 0;
    for ( auto const& tracktraj : (*TrackTrajHandle) ) {
      caf::SRTrack& srtrk = rec.trk[iTrackTraj];

      std::vector<TVector3> temp = tracktraj.getFWDTrajectory();
      for(size_t i = 0; i < temp.size(); i++) {
        srtrk.fwd.traj.emplace_back(temp.at(i).X(),
                                    temp.at(i).Y(),
                                    temp.at(i).Z());
      }

      temp = tracktraj.getBAKTrajectory();
      for(size_t i = 0; i < temp.size(); i++) {
        srtrk.bak.traj.emplace_back(temp.at(i).X(),
                                    temp.at(i).Y(),
                                    temp.at(i).Z());
      }
      iTrackTraj++;
    }
  }

  // save Vertex and Track-Vertex association info
  if (fWriteVertices) {
    size_t iVertex = 0;
    for ( auto const& vertex : (*VertexHandle) ) {
      caf::SRVertex srvtx;
      srvtx.id = vertex.getIDNumber();
      srvtx.x = vertex.Position()[0];
      srvtx.y = vertex.Position()[1];
      srvtx.z = vertex.Position()[2];
      srvtx.t = vertex.Time();

      srvtx.ntrks = 0;
      if ( findManyTrackEnd->isValid() ) {
        srvtx.ntrks = findManyTrackEnd->at(iVertex).size();
      }

      int vertexCharge = 0;
      for (int iVertexedTrack=0; iVertexedTrack < srvtx.ntrks; ++iVertexedTrack) {
        // TODO can we just make a vector of tracks belonging to this vtx?
        caf::SRVertexAssn srassn;
        srassn.vtxid = vertex.getIDNumber();

        // Get this vertexed track.
        rec::Track track = *(findManyTrackEnd->at(iVertex).at(iVertexedTrack));
        srassn.trkid = track.getIDNumber();

        // Get the end of the track in the vertex.  It isn't that odd for the end
        // of the track not used in the vertex to be closer to the vertex than the
        // one actually used; you might have a very short stub track in a 3 track
        // vertex with small opening angles and the other 2 tracks might pull the
        // vertex some distance towards the far end of the stub track
        rec::TrackEnd fee = *(findManyTrackEnd->data(iVertex).at(iVertexedTrack));
        // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
        srassn.trkend = fee;

        rec.assn.vtx.push_back(srassn);

        if (fee==rec::TrackEndBeg) {
          vertexCharge += track.ChargeBeg();
        } else {
          vertexCharge += track.ChargeEnd();
        }
      }
      srvtx.Q = vertexCharge;
      ++iVertex;

      rec.vtx.push_back(srvtx);
    } // end loop over VertexHandle
  }

  rec.nvtx = rec.vtx.size();
  rec.assn.nvtx = rec.assn.vtx.size();

  // save Vee and Track-Vee association info
  if (fWriteVees) {
    size_t iVee = 0;
    for ( auto const& vee : (*VeeHandle) ) {
      caf::SRVee srvee;
      srvee.id = vee.getIDNumber();
      srvee.x = vee.Position()[0];
      srvee.y = vee.Position()[1];
      srvee.z = vee.Position()[2];
      srvee.t = vee.Time();
      srvee.Kpipi.p.x = vee.FourMomentum(0).X();
      srvee.Kpipi.p.y = vee.FourMomentum(0).Y();
      srvee.Kpipi.p.z = vee.FourMomentum(0).Z();
      srvee.Kpipi.E = vee.FourMomentum(0).E();
      srvee.Kpipi.m = vee.FourMomentum(0).M();
      srvee.Lppi.p.x = vee.FourMomentum(1).X();
      srvee.Lppi.p.y = vee.FourMomentum(1).Y();
      srvee.Lppi.p.z = vee.FourMomentum(1).Z();
      srvee.Lppi.E = vee.FourMomentum(1).E();
      srvee.Lppi.m = vee.FourMomentum(1).M();
      srvee.Lpip.p.x = vee.FourMomentum(2).X();
      srvee.Lpip.p.y = vee.FourMomentum(2).Y();
      srvee.Lpip.p.z = vee.FourMomentum(2).Z();
      srvee.Lpip.E = vee.FourMomentum(2).E();
      srvee.Lpip.m = vee.FourMomentum(2).M();

      rec.vee.push_back(srvee);

      int nVeeTracks = 0;
      if ( findManyVeeTrackEnd->isValid() ) {
        nVeeTracks = findManyVeeTrackEnd->at(iVee).size();
      }

      for (int iVeeTrack=0; iVeeTrack<nVeeTracks; ++iVeeTrack) {
        caf::SRVeeAssn srassn;
        srassn.veeid = vee.getIDNumber();

        // Get this vertexed track.
        rec::Track track = *(findManyVeeTrackEnd->at(iVee).at(iVeeTrack));
        srassn.trkid = track.getIDNumber();

        rec::TrackEnd fee = *(findManyVeeTrackEnd->data(iVee).at(iVeeTrack));
        // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
        srassn.trkend = fee;

        rec.assn.vee.push_back(srassn);
      }
      ++iVee;
    } // end loop over VeeHandle
  }

  rec.nvee = rec.vee.size();
  rec.assn.nvee = rec.assn.vee.size();

  // save Cluster info
  if (fWriteCaloClusters) {
    for ( auto const& cluster : (*RecoClusterHandle) ) {
      caf::SRECalCluster srclust;
      srclust.id = cluster.getIDNumber();
      srclust.nhits = cluster.CalorimeterHits().size();
      srclust.E = cluster.Energy();
      srclust.t = cluster.Time();
      srclust.TimeDiffFirstLast = cluster.TimeDiffFirstLast();
      srclust.x = cluster.Position()[0];
      srclust.y = cluster.Position()[1];
      srclust.z = cluster.Position()[2];
      srclust.theta = cluster.ITheta();
      srclust.phi = cluster.IPhi();
      srclust.pid = cluster.ParticleID();
      // srclust.Shape = cluster.Shape();
      srclust.mainAxis.x = cluster.EigenVectors()[0];
      srclust.mainAxis.y = cluster.EigenVectors()[1];
      srclust.mainAxis.z = cluster.EigenVectors()[2];

      // Matching MCParticle info
      std::vector<std::pair<simb::MCParticle*,float>> trakt;
      trakt = BackTrack->ClusterToMCParticles( const_cast<rec::Cluster*>(&cluster) );
      srclust.mcidx = -1;
      if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
        srclust.mcidx = TrackIdToIndex[trakt[0].first->TrackId()];
      }
      if (srclust.mcidx > -1) {
        srclust.mcfrac = trakt[0].second;
      } else {
        srclust.mcfrac = 0;
      }

      rec.clust.ecal.push_back(srclust);
    }

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      for ( auto const& cluster : (*RecoClusterMuIDHandle) ) {
        caf::SRMuIDCluster srclust;
        srclust.id = cluster.getIDNumber();
        srclust.nhits = cluster.CalorimeterHits().size();
        srclust.E = cluster.Energy();
        srclust.t = cluster.Time();
        srclust.TimeDiffFirstLast = cluster.TimeDiffFirstLast();
        srclust.x = cluster.Position()[0];
        srclust.y = cluster.Position()[1];
        srclust.z = cluster.Position()[2];
        srclust.theta = cluster.ITheta();
        srclust.phi = cluster.IPhi();
        srclust.pid = cluster.ParticleID();
        // srclust.Shape = cluster.Shape();
        srclust.mainAxis.x = cluster.EigenVectors()[0];
        srclust.mainAxis.y = cluster.EigenVectors()[1];
        srclust.mainAxis.z = cluster.EigenVectors()[2];

        // Matching MCParticle info
        std::vector<std::pair<simb::MCParticle*,float>> trakt;
        trakt = BackTrack->ClusterToMCParticles( const_cast<rec::Cluster*>(&cluster) );
        srclust.mcidx = -1;
        if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
          srclust.mcidx = TrackIdToIndex[trakt[0].first->TrackId()];
        }
        if (srclust.mcidx > -1) {
          srclust.mcfrac = trakt[0].second;
        } else {
          srclust.mcfrac = 0;
        }

        rec.clust.muid.push_back(srclust);
      }
    }
  }

  rec.clust.necal = rec.clust.ecal.size();
  rec.clust.nmuid = rec.clust.muid.size();

  // Write info for ECAL-matched tracks
  if (fWriteMatchedTracks) {
    size_t iCluster = 0;
    for ( auto const& cluster : (*RecoClusterHandle) ) {
      int nCALedTracks(0);
      if ( findManyCALTrackEnd->isValid() ) {
        nCALedTracks = findManyCALTrackEnd->at(iCluster).size();
      }
      for (int iCALedTrack=0; iCALedTrack<nCALedTracks; ++iCALedTrack) {
        caf::SRClusterAssn srassn;
        srassn.clustid = cluster.getIDNumber();
        rec::Track track  = *(findManyCALTrackEnd->at(iCluster).at(iCALedTrack));
        srassn.trkid = track.getIDNumber();

        rec::TrackEnd fee = *(findManyCALTrackEnd->data(iCluster).at(iCALedTrack));
        srassn.trkend = fee;    // The rec::TrackEnd (see Track.h) that extrapolated to cluster
        rec.assn.ecalclust.push_back(srassn);
      }
      iCluster++;
    }
  } // end branch on fWriteCaloInfo

  rec.assn.necalclust = rec.assn.ecalclust.size();

  return;
} // end :CAFMaker::FillVectors



//==============================================================================
//==============================================================================
//==============================================================================
// Process ionization.  Eventually this moves into the reco code.
void gar::CAFMaker::processIonizationInfo(rec::TrackIoniz& ion, float ionizeTruncate,
                                          float& forwardIonVal, float& backwardIonVal) {

  // NO CALIBRATION SERVICE FOR NOW

  std::vector<std::pair<float,float>> SigData = ion.getFWD_dSigdXs();
  forwardIonVal = processOneDirection(SigData, ionizeTruncate);

  SigData = ion.getBAK_dSigdXs();
  backwardIonVal = processOneDirection(SigData, ionizeTruncate);

  return;
}



float gar::CAFMaker::processOneDirection(std::vector<std::pair<float,float>> SigData, float ionizeTruncate) {

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
// Coherent pion analysis specific code
float gar::CAFMaker::computeT( simb::MCTruth theMCTruth ) {
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



//==============================================================================
//==============================================================================
//==============================================================================
std::vector< std::pair<int, float> > gar::CAFMaker::processPIDInfo( float p ) {

  std::vector<std::string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};
  std::vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
  std::vector< std::pair<int, float> > pid;
  pid.resize(6);

  int qclosest = 0;
  float dist = 100000000.;
  CLHEP::RandFlat FlatRand(fEngine);

  for (int q = 0; q < 501; ++q) {
    //Check the title and the reco momentum take only the one that fits
    std::string fulltitle = m_pidinterp[q]->GetTitle();
    unsigned first = fulltitle.find("=");
    unsigned last = fulltitle.find("GeV");
    std::string substr = fulltitle.substr(first+1, last - first-1);
    float pidinterp_mom = std::atof(substr.c_str());
    //calculate the distance between the bin and mom, store the q the closest
    float disttemp = std::abs(pidinterp_mom - p);

    if( disttemp < dist ) {
      dist = disttemp;
      qclosest = q;
    }
  } // closes the "pidmatrix" loop

    //Compute all the probabities for each type of true to reco
    //loop over the columns (true pid)
  for (int pidm = 0; pidm < 6; ++pidm) {

    //loop over the columns (true pid)
    std::vector< std::pair<float, std::string> > v_prob;

    //loop over the rows (reco pid)
    for (int pidr = 0; pidr < 6; ++pidr) {
      std::string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
      float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1, pidr+1);
      //Need to check random number value and prob value then associate the recopdg to the reco prob
      v_prob.push_back( std::make_pair(prob, recoparticlename) );
    }

    //Compute the pid from it
    if (v_prob.size() > 1) {
      //Order the vector of prob
      std::sort(v_prob.begin(), v_prob.end());
      //Throw a random number between 0 and 1
      float random_number = FlatRand.fire();
      //Make cumulative sum to get the range
      std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const P& _x, const P& _y){return P(_x.first + _y.first, _y.second);});

      for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++) {
        if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first ) {
          pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) ), v_prob.at(ivec+1).first) );
        }
      }
    } else {
      pid.push_back( std::make_pair(pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) ), v_prob.at(0).first) );
    }
  }

  //return a vector of pid and prob
  return pid;
}

DEFINE_ART_MODULE(gar::CAFMaker)

