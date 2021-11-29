////////////////////////////////////////////////////////////////////////////////
// Class:        MatchingPerformance
// Plugin Type: analyzer (art v2_11_02)
// File:        MatchingPerformance_module.cc
//
// Generated 25 Sep 2020 by Leo Bellantoni
// 
// For each reco track:
//     Use backtracker to get matching MCparticles.  Take charged & stable
//         matching MCParticle with largest energy fraction for the track
//     Only use MCParticles which go eventually into the ECAL?
//     Use backtracker info to determine "starting" end of the track.
//     Record chi-2 and NTPC clusters in the track (and other stuff)
//     List all the clusters that match back to the MCParticle, and all
//         that match to the track.  Tabulate the number and energy sum
//             For clusters on track & not MCParticle
//             For clusters on MCParticle & not on track
//             For clusters on both
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

#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "MCCheater/BackTracker.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "MCCheater/BackTracker.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "RecoAlg/TrackPropagator.h"


#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "Math/ProbFuncMathCore.h"

#include <string>
#include <vector>
#include <unordered_map>





namespace gar {

  class MatchingPerformance : public art::EDAnalyzer {
  public:
    explicit MatchingPerformance(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    MatchingPerformance(MatchingPerformance const &) = delete;
    MatchingPerformance(MatchingPerformance &&) = delete;
    MatchingPerformance & operator = (MatchingPerformance const &) = delete;
    MatchingPerformance & operator = (MatchingPerformance &&) = delete;

    virtual void beginJob() override;

    // Required functions.
    void analyze(art::Event const & e) override;



  private:
    void ClearVectors();
    bool FillVectors(art::Event const & e);



    // Position of TPC from geometry service; 1 S Boston Ave.
    float ItsInTulsa[3];

    // Input data labels
    std::string fGeantLabel;
    std::string fTrackLabel;
    std::string fClusterLabel;
    std::string fECALAssnLabel;
    std::string fInstanceLabelECAL; ///< Instance name for ECAL

    // the analysis tree
    TTree *fTree;

    // Geometry service
    const gar::geo::GeometryCore* fGeo;
    // Backtracker service
    cheat::BackTrackerCore* fBack;
    // PDG database from ROOT
    TDatabasePDG* pdgInstance;
    // Detector properties
    const detinfo::DetectorProperties* fDetProp;
    const detinfo::DetectorClocks*     fClocks;
    float                              fMaxXdisplacement;

    // Keepin an eye on it all
    int   fVerbosity;
    int   fClusterDirNhitCut;     ///< Do not use cluster direction unless you have this many hits or more
    float fMinPvalCut;            ///< Do not analyse track fits with this pVal or less
    TH1F* chargeFrac;             ///< Fraction total ionization from theMCPart in any given track
    float fMinIonFracCut;         ///< Do not analyse track fits with low ionization fraction



    // The analysis tree. Use Rtypes.h here, as these data get used by root
    // global event info
    Int_t    fRun;
    Int_t    fSubRun;
    Int_t    fEvent;

    // This analysis starts from reco'd tracks
    std::vector<ULong64_t>          fTrackIDNumber;
    std::vector<Float_t>            fTrackX;
    std::vector<Float_t>            fTrackY;
    std::vector<Float_t>            fTrackZ;
    std::vector<Float_t>            fTrackPX;
    std::vector<Float_t>            fTrackPY;
    std::vector<Float_t>            fTrackPZ;
    std::vector<Float_t>            fTrackPmag;
    std::vector<Int_t>              fTrackQ;
    std::vector<Float_t>            fTrackLen;
    std::vector<Float_t>            fTrackChi2;
    std::vector<Int_t>              fNTPCClustersOnTrack;
    std::vector<Float_t>            fTrackPval;
    std::vector<Double_t>           fTrackTime;

    // Which match to qualified MCParticles
    std::vector<Int_t>              fMCPDG;
    std::vector<Int_t>              fMCPDGMother;
    std::vector<Float_t>            fMCPStartX;
    std::vector<Float_t>            fMCPStartY;
    std::vector<Float_t>            fMCPStartZ;
    std::vector<Float_t>            fMCPStartPX;
    std::vector<Float_t>            fMCPStartPY;
    std::vector<Float_t>            fMCPStartPZ;
    std::vector<Float_t>            fMCPStartE;
    std::vector<std::string>        fMCPProc;
    std::vector<std::string>        fMCPEndProc;
    std::vector<Float_t>            fMCPTime;

    // Track - ECAL energy info
    std::vector<Float_t>            fEassocNoInTruth;
    std::vector<Int_t>              fNassocNoInTruth;
    std::vector<Float_t>            fEunassocInTruth;
    std::vector<Int_t>              fNunassocInTruth;
    std::vector<Float_t>            fEisassocInTruth;
    std::vector<Int_t>              fNisassocInTruth;

    // The matching variables here; they will have one entry for each
    // cluster, for each "nice" track end; all other vectors here will
    // have just one entry for each track
    std::vector<Float_t>            fRadClusTrackTru;
    std::vector<Float_t>            fRadClusTrackFal;
    std::vector<Float_t>            fXClusTrackOver25Tru;
    std::vector<Float_t>            fXClusTrackOver25Fal;
    std::vector<Float_t>            fXClusTrackUnder25Tru;
    std::vector<Float_t>            fXClusTrackUnder25Fal;
    std::vector<Float_t>            fRPhiClusTrackTru;
    std::vector<Float_t>            fRPhiClusTrackFal;
    std::vector<Float_t>            fDotClustTrackTru;
    std::vector<Float_t>            fDotClustTrackFal;
    std::vector<Float_t>            fDistClustTrackTru;
    std::vector<Float_t>            fDistClustTrackFal;
  };
}





//==============================================================================
//==============================================================================
//==============================================================================
gar::MatchingPerformance::MatchingPerformance(fhicl::ParameterSet const & p) : EDAnalyzer(p) {

  fGeantLabel        = p.get<std::string>("GEANTLabel",       "geant");
  fTrackLabel        = p.get<std::string>("TrackLabel",       "track");
  fClusterLabel      = p.get<std::string>("ClusterLabel",     "calocluster");
  fECALAssnLabel     = p.get<std::string>("ECALAssnLabel",    "trkecalassn");
  fInstanceLabelECAL = p.get<std::string>("InstanceLabelCalo","ECAL");
  fVerbosity         = p.get<int>        ("Verbosity",         0);
  fClusterDirNhitCut = p.get<int>        ("ClusterDirNhitCut", 5);
  fMinPvalCut        = p.get<float>      ("MinPvalCut",        0.0);
  fMinIonFracCut     = p.get<float>      ("MinIonFracCut",     0.50);

  pdgInstance = TDatabasePDG::Instance();

  consumesMany<std::vector<simb::MCTruth> >();
  consumes<std::vector<simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<rec::Track> >(fTrackLabel);
  consumes<std::vector<rec::Cluster> >(fClusterLabel);
  consumes<art::Assns<rec::Cluster, rec::Track>>(fECALAssnLabel);

  return;
} // end constructor





//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::beginJob() {

  fGeo = gar::providerFrom<geo::GeometryGAr>();
  ItsInTulsa[0] = fGeo->TPCXCent();        // 1 S Boston Ave
  ItsInTulsa[1] = fGeo->TPCYCent();
  ItsInTulsa[2] = fGeo->TPCZCent();



  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("GArPerfTree","GArPerfTree");

  if (fVerbosity>0) {
    chargeFrac = tfs->make<TH1F>("chargeFrac", "Fraction of MCPart's ionization in track",
                                 101,0.00,1.01);
  }



  fTree->Branch("Run",                    &fRun,              "Run/I");
  fTree->Branch("SubRun",                 &fSubRun,           "SubRun/I");
  fTree->Branch("Event",                  &fEvent,            "Event/I");

  fTree->Branch("TrackIDNumber",          &fTrackIDNumber);
  fTree->Branch("TrackX",                 &fTrackX);
  fTree->Branch("TrackY",                 &fTrackY);
  fTree->Branch("TrackZ",                 &fTrackZ);
  fTree->Branch("TrackPX",                &fTrackPX);
  fTree->Branch("TrackPY",                &fTrackPY);
  fTree->Branch("TrackPZ",                &fTrackPZ);
  fTree->Branch("TrackPmag",              &fTrackPmag);
  fTree->Branch("TrackQ",                 &fTrackQ);
  fTree->Branch("TrackLen",               &fTrackLen);
  fTree->Branch("TrackChi2",              &fTrackChi2);
  fTree->Branch("NTPCClustersOnTrack",    &fNTPCClustersOnTrack);
  fTree->Branch("TrackPval",              &fTrackPval);
  fTree->Branch("TrackTime",              &fTrackTime);

  fTree->Branch("PDG",                    &fMCPDG);
  fTree->Branch("PDGMother",              &fMCPDGMother);
  fTree->Branch("MCPStartX",              &fMCPStartX);
  fTree->Branch("MCPStartY",              &fMCPStartY);
  fTree->Branch("MCPStartZ",              &fMCPStartZ);
  fTree->Branch("MCPStartPX",             &fMCPStartPX);
  fTree->Branch("MCPStartPY",             &fMCPStartPY);
  fTree->Branch("MCPStartPZ",             &fMCPStartPZ);
  fTree->Branch("MCPTime",                &fMCPTime);

  fTree->Branch("EassocNoInTruth",        &fEassocNoInTruth);
  fTree->Branch("NassocNoInTruth",        &fNassocNoInTruth);
  fTree->Branch("EunassocInTruth",        &fEunassocInTruth);
  fTree->Branch("NunassocInTruth",        &fNunassocInTruth);
  fTree->Branch("EisassocInTruth",        &fEisassocInTruth);
  fTree->Branch("NisassocInTruth",        &fNisassocInTruth);

  fTree->Branch("RadClusTrackTru",        &fRadClusTrackTru);
  fTree->Branch("RadClusTrackFal",        &fRadClusTrackFal);
  fTree->Branch("XClusTrackOver25Tru",    &fXClusTrackOver25Tru);
  fTree->Branch("XClusTrackOver25Fal",    &fXClusTrackOver25Fal);
  fTree->Branch("XClusTrackUnder25Tru",   &fXClusTrackUnder25Tru);
  fTree->Branch("XClusTrackUnder25Fal",   &fXClusTrackUnder25Fal);
  fTree->Branch("RPhiClusTrackTru",       &fRPhiClusTrackTru);
  fTree->Branch("RPhiClusTrackFal",       &fRPhiClusTrackFal);
  fTree->Branch("DotClustTrackTru",       &fDotClustTrackTru);
  fTree->Branch("DotClustTrackFal",       &fDotClustTrackFal);
  fTree->Branch("DistClustTrackTru",      &fDistClustTrackTru);
  fTree->Branch("DistClustTrackFal",      &fDistClustTrackFal);
  return;
}  // End of :MatchingPerformance::beginJob





//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::analyze(art::Event const & event) {

  // Need this here because the clock service is invoked at preProcessEvent.
  fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
  fClocks  = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();
  fMaxXdisplacement =
    fDetProp->DriftVelocity(fDetProp->Efield(),fDetProp->Temperature())
    *fClocks->SpillLength();

  // Is the backtracker good to go?  Need to fill fBack on per-event basis
  // not in beginJob()
  cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
  fBack = const_cast<cheat::BackTrackerCore*>(const_bt);
  if ( !fBack->HasTracks() || !fBack->HasClusters() ) {
    throw cet::exception("MatchingPerformance") << " Not enough backtracker info"
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  ClearVectors();
  if ( FillVectors(event) ) {
    fTree->Fill();
  }

  return;
}





//==============================================================================
//==============================================================================
//==============================================================================
void gar::MatchingPerformance::ClearVectors() {
  // clear out all our vectors

  fTrackIDNumber.clear();
  fTrackX.clear();
  fTrackY.clear();
  fTrackZ.clear();
  fTrackPX.clear();
  fTrackPY.clear();
  fTrackPZ.clear();
  fTrackPmag.clear();
  fTrackQ.clear();
  fTrackLen.clear();
  fTrackChi2.clear();
  fNTPCClustersOnTrack.clear();
  fTrackPval.clear();
  fTrackTime.clear();

  fMCPDG.clear();
  fMCPDGMother.clear();
  fMCPStartX.clear();
  fMCPStartY.clear();
  fMCPStartZ.clear();
  fMCPStartPX.clear();
  fMCPStartPY.clear();
  fMCPStartPZ.clear();
  fMCPTime.clear();

  fEassocNoInTruth.clear();
  fNassocNoInTruth.clear();
  fEunassocInTruth.clear();
  fNunassocInTruth.clear();
  fEisassocInTruth.clear();
  fNisassocInTruth.clear();

  fRadClusTrackTru.clear();
  fRadClusTrackFal.clear();
  fXClusTrackOver25Tru.clear();
  fXClusTrackOver25Fal.clear();
  fXClusTrackUnder25Tru.clear();
  fXClusTrackUnder25Fal.clear();
  fRPhiClusTrackTru.clear();
  fRPhiClusTrackFal.clear();
  fDotClustTrackTru.clear();
  fDotClustTrackFal.clear();
  fDistClustTrackTru.clear();
  fDistClustTrackFal.clear();

  return;
} // end :MatchingPerformance::ClearVectors





//==============================================================================
//==============================================================================
//==============================================================================
bool gar::MatchingPerformance::FillVectors(art::Event const& event) {

  // =============  Get art handles ==========================================
  // Get handles for Tracks, clusters, associations between them
  auto TrackHandle = event.getHandle< std::vector<rec::Track> >(fTrackLabel);
  if (!TrackHandle) {
    throw cet::exception("MatchingPerformance") << " No rec::Track"
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  art::InputTag ecalclustertag(fClusterLabel, fInstanceLabelECAL);
  auto ClusterHandle = event.getHandle< std::vector<rec::Cluster> >(ecalclustertag);
  if (!ClusterHandle) {
    throw cet::exception("MatchingPerformance") << " No rec::Cluster branch."
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  art::FindMany<rec::Cluster, rec::TrackEnd>* findManyTrackEndCAL = NULL;
  art::InputTag ecalassntag(fECALAssnLabel, fInstanceLabelECAL);
  findManyTrackEndCAL = new art::FindMany<rec::Cluster, rec::TrackEnd>
    (TrackHandle,event,ecalassntag);
  if ( !findManyTrackEndCAL->isValid() ) {
    throw cet::exception("MatchingPerformance") << " Bad TrackEnd-ECAL Assn."
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }
 
  // Get handles for MCinfo, also good for MCPTrajectory.  Want to get all
  // MCTruths, regardless of generator label.

  auto mctruthHandles = event.getMany<std::vector<simb::MCTruth> >();

  if (mctruthHandles.size()!=1) {
    throw cet::exception("MatchingPerformance") << " Need just 1 simb::MCTruth"
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }

  auto MCPHandle = event.getHandle< std::vector<simb::MCParticle> >(fGeantLabel);
  if (!MCPHandle) {
    throw cet::exception("MatchingPerformance") << " No simb::MCParticle"
                                                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }



  // Need something to look at!
  if ( TrackHandle->size()==0 || ClusterHandle->size()==0 ) return false;

  // Make the usual map for the MC info
  typedef int TrkId;
  std::unordered_map<TrkId, Int_t> TrackIdToIndex;
  Int_t index = 0;
  for ( auto const& mcp : (*MCPHandle) ) {
    int TrackId = mcp.TrackId();
    TrackIdToIndex[TrackId] = index++;
  }

  // Will need a vector of art::Ptrs to all the ECAL clusters for the
  // backtracker calls.
  std::vector<art::Ptr<rec::Cluster>> ClusterGrabber;
  art::PtrMaker<rec::Cluster> makeClusterPtr(event,ClusterHandle.id());
  for (size_t iCluster=0; iCluster<ClusterHandle->size(); ++iCluster) {
    art::Ptr<rec::Cluster> aPtr = makeClusterPtr(iCluster);
    ClusterGrabber.push_back(aPtr);
  }
  // Will need a vector<art::Ptr<rec::Track>> of tracks also
  std::vector<art::Ptr<rec::Track>> TrackGrabber;
  art::PtrMaker<rec::Track> makeTrackPtr(event,TrackHandle.id());
  for (size_t iTrack=0; iTrack<TrackHandle->size(); ++iTrack) {
    art::Ptr<rec::Track> aPtr = makeTrackPtr(iTrack);
    TrackGrabber.push_back(aPtr);
  }





  // =============  Pull art's handles =======================================
  fRun    = event.run();
  fSubRun = event.subRun();
  fEvent  = event.id().event();



  // Require tracks match to right kind of MC first
  struct niceNice {rec::Track recoTrk; simb::MCParticle simiTrk;};
  std::vector<niceNice> niceTracks;
  for ( auto const& track : (*TrackHandle) ) {

    /// Want to be able to know what you've matched back to at MCTruth
    std::vector<std::pair<simb::MCParticle*,float>> MCsfromTrack;
    gar::rec::Track* nonconst_track = const_cast<gar::rec::Track*>(&track);
    MCsfromTrack = fBack->TrackToMCParticles(nonconst_track);
    int nMCsfromTrack = MCsfromTrack.size();
    if (nMCsfromTrack!=1) continue;

    simb::MCParticle theMCPart  = *(MCsfromTrack[0].first);
    TParticlePDG* particle = pdgInstance->GetParticle(theMCPart.PdgCode());
    if (particle==NULL         ) continue;   // What causes this err? If anything.
    if (particle->Charge()==0.0) continue;
    if (particle->Stable()==0  ) continue;

    // Does theMCParticle go out to the ECAL?
    int nTraj = theMCPart.NumberTrajectoryPoints();
    TVector3 doink;            bool whackThatECAL = false;
    for (int iTraj=0; iTraj<nTraj; ++iTraj) {
      doink.SetXYZ(theMCPart.Vx(iTraj),theMCPart.Vy(iTraj),theMCPart.Vz(iTraj));
      if ( fGeo->PointInECALBarrel(doink) || fGeo->PointInECALEndcap(doink) ) {
        whackThatECAL = true;
        break;
      }
    }
    if (!whackThatECAL) continue;

    // This track is nice
    niceNice tmp = {track, theMCPart};
    niceTracks.push_back(tmp);
  }



  // Anti-stubbiness cut.
  std::vector<niceNice>::iterator iNiceTrk=niceTracks.begin();
  for (; iNiceTrk!=niceTracks.end(); ++iNiceTrk) {
    simb::MCParticle theMCPart = iNiceTrk->simiTrk;
    std::vector<art::Ptr<rec::Track>> tracksFromSameMCP =
      fBack->MCParticleToTracks(&theMCPart,TrackGrabber);

    // Get ionization fraction for every track with contribution from this MCP
    std::vector<float> ionzFromSameMCP;    ionzFromSameMCP.resize(tracksFromSameMCP.size());
    for (size_t iTrkFromSame=0; iTrkFromSame<tracksFromSameMCP.size(); ++iTrkFromSame) {
      rec::Track tmpTrk = *(tracksFromSameMCP[iTrkFromSame]);
      std::vector<art::Ptr<rec::Hit>> hitList = fBack->TrackToHits(&tmpTrk);
      float thisTracksI = 0;
      for (auto& ptrHit : hitList) thisTracksI += ptrHit->Signal();
      ionzFromSameMCP[iTrkFromSame] = thisTracksI;
    }
    float ionzTotSameMCP = 0;
    for (size_t iIonSame=0; iIonSame<ionzFromSameMCP.size(); ++iIonSame) 
      ionzTotSameMCP += ionzFromSameMCP[iIonSame];
    for (size_t iIonSame=0; iIonSame<ionzFromSameMCP.size(); ++iIonSame) {
      ionzFromSameMCP[iIonSame] /= ionzTotSameMCP;
      if (fVerbosity>0) chargeFrac->Fill(ionzFromSameMCP[iIonSame]);
    }

    // Cut on ion fraction in niceTracks
    for (size_t iTrkFromSame=0; iTrkFromSame<tracksFromSameMCP.size(); ++iTrkFromSame) {
      if ( *(tracksFromSameMCP[iTrkFromSame])==iNiceTrk->recoTrk
           && ionzFromSameMCP[iTrkFromSame] <fMinIonFracCut) {
        niceTracks.erase(iNiceTrk);
        --iNiceTrk;
        break;
      }
    }
  }


  iNiceTrk = niceTracks.begin();
  std::vector<niceNice>::iterator jNiceTrk = iNiceTrk +1;
  for (; iNiceTrk!=niceTracks.end(); ++iNiceTrk) {
    for (; jNiceTrk!=niceTracks.end(); ++jNiceTrk) {
      if ( iNiceTrk->simiTrk.TrackId() == jNiceTrk->simiTrk.TrackId() ){
        niceTracks.erase(jNiceTrk);
        --jNiceTrk;
      }
    }
  }



  for (size_t iNiceTrk_local=0; iNiceTrk_local<niceTracks.size(); ++iNiceTrk_local) {
    rec::Track       track      = niceTracks[iNiceTrk_local].recoTrk;
    simb::MCParticle theMCPart  = niceTracks[iNiceTrk_local].simiTrk;



    // Examine matched clusters on downstream ends of the track unless 
    // MC vertex in gas - don't test the matching quality then.
    // Note: rec::TrackEndBeg is at the front of the deque here!
    std::deque<rec::TrackEnd> endList = {rec::TrackEndBeg,rec::TrackEndEnd};

    TVector3 positionMCP = theMCPart.Position(0).Vect();
    if (fGeo->PointInGArTPC(positionMCP)) {
      // Pick the end where the particle ended.
      // Don't include the drift distance in matching.
      float distStart = std::hypot(track.Vertex()[1] -positionMCP[1],
                                   track.Vertex()[2] -positionMCP[2]);
      float distEnd   = std::hypot(track.End()[1]    -positionMCP[1],
                                   track.End()[2]    -positionMCP[2]);
      if ( distStart < distEnd ) {
        endList.pop_front();
      } else {
        endList.pop_back();
      }
    }

    // Having selected now the end we wil extrapolate,
    std::deque<rec::TrackEnd>::iterator iEnd = endList.begin();
    // Record info about the track at extrapolated end - need trackPar, 
    // TrackEnd later and they must be in MPD coordinates.
    fTrackIDNumber.push_back(track.getIDNumber());
    Float_t saveChi2;
    float trackPar[5];      float trackEnd[3];

    if ( (*iEnd)==rec::TrackEndBeg ) {
      for (size_t i=0; i<5; ++i) trackPar[i] = track.TrackParBeg()[i];
      for (size_t i=0; i<3; ++i) trackEnd[i] = track.Vertex()[i];
      fTrackX.push_back   ( track.Vertex()[0] -ItsInTulsa[0] );
      fTrackY.push_back   ( track.Vertex()[1] -ItsInTulsa[1] );
      fTrackZ.push_back   ( track.Vertex()[2] -ItsInTulsa[2] );
      fTrackPX.push_back  (-track.Momentum_beg()*track.VtxDir()[0] );
      fTrackPY.push_back  (-track.Momentum_beg()*track.VtxDir()[1] );
      fTrackPZ.push_back  (-track.Momentum_beg()*track.VtxDir()[2] );
      fTrackPmag.push_back( track.Momentum_beg() );
      fTrackQ.push_back   ( track.ChargeEnd() );
      fTrackLen.push_back ( track.LengthBackward() );
      saveChi2 = track.ChisqBackward();
      fTrackChi2.push_back( saveChi2 );
    } else {
      for (size_t i=0; i<5; ++i) trackPar[i] = track.TrackParEnd()[i];
      for (size_t i=0; i<3; ++i) trackEnd[i] = track.End()[i];
      fTrackX.push_back   ( track.End()[0] -ItsInTulsa[0] );
      fTrackY.push_back   ( track.End()[1] -ItsInTulsa[1] );
      fTrackZ.push_back   ( track.End()[2] -ItsInTulsa[2] );
      fTrackPX.push_back  (-track.Momentum_end()*track.EndDir()[0] );
      fTrackPY.push_back  (-track.Momentum_end()*track.EndDir()[1] );
      fTrackPZ.push_back  (-track.Momentum_end()*track.EndDir()[2] );
      fTrackPmag.push_back( track.Momentum_end() );
      fTrackQ.push_back   ( track.ChargeBeg() );
      fTrackLen.push_back ( track.LengthForward() );        
      saveChi2 = track.ChisqForward();
      fTrackChi2.push_back( saveChi2);

    }

    // Apply the pVal cut
    Int_t nHits = track.NHits();
    fNTPCClustersOnTrack.push_back(nHits);
    Float_t pVal = ROOT::Math::chisquared_cdf_c(saveChi2,nHits-5);
    if (pVal <= fMinPvalCut) return false;
    fTrackPval.push_back(pVal);
    fTrackTime.push_back(track.Time());

    trackPar[0] -=ItsInTulsa[1];        trackPar[1] -=ItsInTulsa[2];
    for (int i=0; i<3; ++i) trackEnd[i] -= ItsInTulsa[i];

    // Record some info about the matching MC
    fMCPDG.push_back(theMCPart.PdgCode());
    TrkId momTrkId = theMCPart.Mother();
    int momPDG = 0;
    if (momTrkId>0) {
      int momIndex = TrackIdToIndex[momTrkId];
      momPDG   = (*MCPHandle).at(momIndex).PdgCode();
    }
    fMCPDGMother.push_back(momPDG);

    fMCPStartX.push_back(positionMCP.X());
    fMCPStartY.push_back(positionMCP.Y());
    fMCPStartZ.push_back(positionMCP.Z());
    const TLorentzVector& momentumMCP = theMCPart.Momentum(0);
    fMCPStartPX.push_back(momentumMCP.Px());
    fMCPStartPY.push_back(momentumMCP.Py());
    fMCPStartPZ.push_back(momentumMCP.Pz());
    fMCPTime.push_back(theMCPart.T());    



    // Now, look for the ECAL-track matching info
    std::vector<rec::Cluster> clustersOnMCP;
    for (rec::Cluster cluster : *ClusterHandle) {
      if (fBack->MCParticleCreatedCluster(&theMCPart,&cluster)) {
        clustersOnMCP.push_back(cluster);
      }
      if (fBack->ClusterCreatedMCParticle(&theMCPart,&cluster)) {
        clustersOnMCP.push_back(cluster);
      }
    }

    // Because of how IsForebearOf works, the ionization in the cluster
    // which is from theMCPart makes theMCPart both a parent and descendent
    // of the MCParticles in the cluster.  So we must remove duplicate 
    // clustersOnMCP, which means sorting them, which means we need a 
    // lambda function, etc. 
    std::sort(clustersOnMCP.begin(),clustersOnMCP.end(),
              [](rec::Cluster a,rec::Cluster b){return a.getIDNumber() > b.getIDNumber();});
    std::vector<rec::Cluster>::iterator itr =
      std::unique(clustersOnMCP.begin(),clustersOnMCP.end());
    clustersOnMCP.resize(std::distance(clustersOnMCP.begin(),itr));




    float eAssocNoInTruth = 0.0;            int nAssocNoInTruth = 0;
    float eUnassocInTruth = 0.0;            int nUnassocInTruth = 0;
    float eIsassocInTruth = 0.0;            int nIsassocInTruth = 0;
    size_t   nCALedClusts =   0;
    std::vector<rec::Cluster  const*> clustersFromAssn;
    std::vector<rec::TrackEnd const*> trackEndsFromAssn;

    // The art::FindMany object findManyTrackEndCAL is indexed over 
    // all the tracks in *TrackHandle (which is bigger than niceTracks).
    size_t iTrack=0;
    for (; iTrack<TrackHandle->size(); ++iTrack) {
      if ( (*TrackHandle)[iTrack] == track ) break;
    }

    // nCALedClusts.size() == trackEndsFromAssn.size() and is the same
    // number regardless of std::deque<rec::TrackEnd>::iterator iEnd but
    // the code reads better here
    findManyTrackEndCAL->get(iTrack, clustersFromAssn,trackEndsFromAssn);
    nCALedClusts = clustersFromAssn.size();

    for ( auto const& cluster : (*ClusterHandle) ) {

      // recompute the matching variables, similar to (if not exactly
      // the same as) in Reco/TPCECALAssociation_module.cc  They will be:
      // inECALBarrel, distRadially, Over25, cutQuantity, dotSee
      float radius = 1.0/trackPar[2];
      float zCent = trackPar[1] - radius*sin(trackPar[3]);
      float yCent = trackPar[0] + radius*cos(trackPar[3]);
      float xClus = cluster.Position()[0] -ItsInTulsa[0];
      float yClus = cluster.Position()[1] -ItsInTulsa[1];
      float zClus = cluster.Position()[2] -ItsInTulsa[2];
      float rClus = std::hypot(zClus,yClus);
      bool  inECALBarrel = fGeo->PointInECALBarrel(cluster.Position());

      float distRadially = std::hypot(zClus-zCent,yClus-yCent) -abs(radius);
      distRadially = abs(distRadially);

      float retXYZ1[3];    float retXYZ2[3];    TVector3 trackXYZ;
      float cutQuantity = -1.0;      bool over25 = false;

      if (inECALBarrel) {
        int errcode = util::TrackPropagator::PropagateToCylinder(
                                                                 trackPar,trackEnd,rClus, 0.0, 0.0, retXYZ1,retXYZ2);
        if ( errcode==0 ) {
          float extrapXerr;
          float transDist1 = std::hypot(retXYZ1[2]-zClus,retXYZ1[1]-yClus);
          float transDist2 = std::hypot(retXYZ2[2]-zClus,retXYZ2[1]-yClus);
          bool looneyExtrap;
          if (transDist1<transDist2) {
            trackXYZ.SetXYZ(retXYZ1[0],retXYZ1[1],retXYZ1[2]);
          } else {
            trackXYZ.SetXYZ(retXYZ2[0],retXYZ2[1],retXYZ2[2]);
          }
          trackXYZ += TVector3(ItsInTulsa);
          looneyExtrap =  !fGeo->PointInECALBarrel(trackXYZ)
            && !fGeo->PointInECALEndcap(trackXYZ);
          trackXYZ -= TVector3(ItsInTulsa);
          if (!looneyExtrap) {
            extrapXerr = trackXYZ.X() -xClus;
            float expected_mean = 0;
            if (trackEnd[0]<-25) {
              expected_mean = +fMaxXdisplacement/2.0;
              over25 = true;
            }
            if (trackEnd[0]>+25) {
              expected_mean = -fMaxXdisplacement/2.0;
              over25 = true;
            }
            cutQuantity = abs(extrapXerr -expected_mean);
          }
        }
      } else {
        // In an endcap.  How many radians in a fMaxXdisplacement?
        float radiansInDrift = trackPar[2]*fMaxXdisplacement
          / tan(trackPar[4]);
        if ( abs(radiansInDrift) >= 2.0*M_PI ) goto angleCut;
        int errcode = util::TrackPropagator::PropagateToX(
                                                          trackPar,trackEnd, xClus, retXYZ1);
        if ( errcode==0 ) {
          trackXYZ.SetXYZ(retXYZ1[0],retXYZ1[1],retXYZ1[2]);
          trackXYZ += TVector3(ItsInTulsa);
          bool looneyExtrap =  !fGeo->PointInECALBarrel(trackXYZ)
            && !fGeo->PointInECALEndcap(trackXYZ);
          trackXYZ -= TVector3(ItsInTulsa);
          if (!looneyExtrap) {
                            
            float angClus  = std::atan2(yClus-yCent,zClus-zCent);
            float angXtrap = std::atan2(retXYZ1[1] -yCent,retXYZ1[2] -zCent) -angClus;
            // angXtrap can indeed be outside of -PI to +PI
            if (angXtrap > +M_PI) angXtrap -= 2.0*M_PI;
            if (angXtrap < -M_PI) angXtrap += 2.0*M_PI;
            cutQuantity = abs(radius*angXtrap);
          }
        }
      }

    angleCut:
      float trackDir[3];
      float dotSee = -2.0;
      int nCells = cluster.CalorimeterHits().size();
      if (nCells >= fClusterDirNhitCut) {
        float xEval = trackXYZ.x();  // Ya not the cluster X
        int errcode = util::TrackPropagator::DirectionX(
                                                        trackPar,trackEnd, xEval,trackDir);
        if (errcode==0) {
          TVector3 clusterDir;
          clusterDir.SetXYZ(cluster.EigenVectors()[0],
                            cluster.EigenVectors()[1],cluster.EigenVectors()[2]);
          dotSee = clusterDir.Dot(trackDir);
        }
      }

      // Look at distance from extrapolated track to cluster.
      float dist3dXtrap = std::hypot(
                                     trackXYZ.x()-xClus, trackXYZ.y()-yClus, trackXYZ.z()-zClus);



      bool clusterOnTrack = false;
      for (size_t iCALedClust=0; iCALedClust<nCALedClusts; ++iCALedClust) {
        if (*trackEndsFromAssn[iCALedClust] != (*iEnd) ) continue;
        if (*clustersFromAssn[iCALedClust] == cluster) {
          clusterOnTrack = true;
          break;
        }
      }

      bool clusterOnMCP = false;
      for (size_t iClusOnMCP = 0; iClusOnMCP<clustersOnMCP.size(); ++iClusOnMCP) {
        if ( clustersOnMCP[iClusOnMCP] == cluster ) {
          clusterOnMCP = true;
          break;
        }
      }



 
      if (!clusterOnMCP ) {
        fRadClusTrackFal.push_back(distRadially);
        if (inECALBarrel) {
          if (over25) {
            fXClusTrackOver25Fal.push_back(cutQuantity);
          } else {
            fXClusTrackUnder25Fal.push_back(cutQuantity);
          }
        } else {
          fRPhiClusTrackFal.push_back(cutQuantity);
        }
        fDotClustTrackFal.push_back(dotSee);
        fDistClustTrackFal.push_back(dist3dXtrap);

        if ( clusterOnTrack ) {
          eAssocNoInTruth += cluster.Energy();
          nAssocNoInTruth += cluster.CalorimeterHits().size();
        }
      } else {
        fRadClusTrackTru.push_back(distRadially);
        if (inECALBarrel) {
          if (over25) {
            fXClusTrackOver25Tru.push_back(cutQuantity);
          } else {
            fXClusTrackUnder25Tru.push_back(cutQuantity);
          }
        } else {
          fRPhiClusTrackTru.push_back(cutQuantity);
        }
        fDotClustTrackTru.push_back(dotSee);
        fDistClustTrackTru.push_back(dist3dXtrap);

        if (!clusterOnTrack && clusterOnMCP ) {
          eUnassocInTruth += cluster.Energy();
          nUnassocInTruth += cluster.CalorimeterHits().size();
        } else {
          eIsassocInTruth += cluster.Energy();
          nIsassocInTruth += cluster.CalorimeterHits().size();
        }
      }
    }

    fEassocNoInTruth.push_back(eAssocNoInTruth);
    fNassocNoInTruth.push_back(nAssocNoInTruth);
    fEunassocInTruth.push_back(eUnassocInTruth);
    fNunassocInTruth.push_back(nUnassocInTruth);
    fEisassocInTruth.push_back(eIsassocInTruth);
    fNisassocInTruth.push_back(nIsassocInTruth);

  } // end loop over TrackHandle

  return true;
} // end MatchingPerformance::FillVectors


DEFINE_ART_MODULE(gar::MatchingPerformance)
