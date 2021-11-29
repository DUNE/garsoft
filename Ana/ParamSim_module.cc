////////////////////////////////////////////////////////////////////////////////
// Class:       ParamSim
// Plugin Type: analyzer
// File:        ParamSim_module.cc
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

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

#include "TFile.h"
#include "TTree.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"

#include <string>
#include <vector>
#include <unordered_map>

namespace gar {

  //Helper Class
  class CAFHelper
  {
  public:
    /* C'tor */
    CAFHelper(const geo::GeometryCore *fGeo, CLHEP::HepRandomEngine& engine);

    /* D'tor */
    ~CAFHelper();

    /* Check if MCP point is in fiducial */
    bool PointInFiducial(const TVector3& point);

    /* Check if MCP point is in TPC Volume (can be outside Fid) */
    bool PointInTPC(const TVector3& point);

    /* Check if MCP point is in the ECAL */
    bool PointInCalo(const TVector3& point);

    /* Check if MCP point is in between the TPC and the ECAL */
    bool PointStopBetween(const TVector3& point);

    /* Check if MCP point is not in the fiducial and not in the ECAL */
    bool isThroughCalo(const TVector3& point);

    /* Check if MCP decayed in calo */
    bool hasDecayedInCalo(const TVector3& point);

    /* Check if MCP is in the Barrel region */
    bool isBarrel(const TVector3& point);

    /* Check if MCP is in the Endcap region */
    bool isEndcap(const TVector3& point);

    float GetRamdomNumber();

    float GaussianSmearing(const float mean, const float sigma);

  private:
    void PrintParameters();

    CLHEP::HepRandomEngine& fEngine;

    //For the fiducial volume
    const double fTPCFidRadius = 222.5;
    const double fTPCFidLength = 215.;

    double fTPCRadius;
    double fTPCLength;
    double fECALBarrelInnerRadius;
    double fECALBarrelOuterRadius;
    double fECALEndcapInnerRadius;
    double fECALEndcapOuterRadius;
    double fECALStartX;
    double fECALEndX;
  };

  //==============================================================================
  inline float CAFHelper::GetRamdomNumber() {
    CLHEP::RandFlat FlatRand(fEngine);
    return FlatRand.fire();;
  }

  //==============================================================================
  inline float CAFHelper::GaussianSmearing(const float mean, const float sigma) {
    CLHEP::RandGauss GausRand(fEngine);
    return GausRand.fire(mean, sigma);;
  }

  //==============================================================================
  CAFHelper::CAFHelper(const geo::GeometryCore *fGeo, CLHEP::HepRandomEngine& engine)
    : fEngine(engine)
  {
    fTPCRadius = fGeo->TPCRadius();
    fTPCLength = fGeo->TPCLength()/2.;
    fECALBarrelInnerRadius = fGeo->GetECALInnerBarrelRadius();
    fECALBarrelOuterRadius = fGeo->GetECALOuterBarrelRadius();
    fECALEndcapInnerRadius = fGeo->GetECALInnerEndcapRadius();
    fECALEndcapOuterRadius = fGeo->GetECALOuterEndcapRadius();
    fECALStartX = fGeo->GetECALEndcapStartX();
    fECALEndX = fGeo->GetECALEndcapOuterX();

    PrintParameters();
  }

  //==============================================================================
  CAFHelper::~CAFHelper()
  {

  }

  //==============================================================================
  void CAFHelper::PrintParameters()
  {
    std::cout << " ==== CAFHelper Parameters ==== " << std::endl;
    std::cout << "TPC Radius " << fTPCRadius << " cm" << std::endl;
    std::cout << "TPC Length " << fTPCLength << " cm" << std::endl;
    std::cout << "ECAL Barrel Inner Radius " << fECALBarrelInnerRadius << " cm" << std::endl;
    std::cout << "ECAL Barrel Outer Radius " << fECALBarrelOuterRadius << " cm" << std::endl;
    std::cout << "ECAL Endcap Inner Radius  " << fECALEndcapInnerRadius << " cm" << std::endl;
    std::cout << "ECAL Endcap Outer Radius  " << fECALEndcapOuterRadius << " cm" << std::endl;
    std::cout << "ECAL Endcap Start X " << fECALStartX << " cm" << std::endl;
    std::cout << "ECAL Endcap End X " << fECALEndX << " cm" << std::endl;
    std::cout << " ==== CAFHelper Parameters ==== " << std::endl;
  }

  //==============================================================================
  inline bool CAFHelper::PointInFiducial(const TVector3& point)
  {
    //TPC Fiducial volume defined as
    //R < 260 cm
    //|X| < 250 cm
    bool isInFiducial = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > fTPCFidRadius ) isInFiducial = false;
    if( r_point < fTPCFidRadius && std::abs(point.X()) > fTPCFidLength ) isInFiducial = false;

    return isInFiducial;
  }

  //==============================================================================
  inline bool CAFHelper::PointInTPC(const TVector3& point)
  {
    //TPC volume defined as
    //R < 260 cm
    //|X| < 250 cm
    if(PointInFiducial(point)) return true;
    bool isInTPC = true;

    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    if( r_point > fTPCRadius ) isInTPC = false;
    if( r_point < fTPCRadius && std::abs(point.X()) > fTPCLength ) isInTPC = false;

    return isInTPC;
  }

  //==============================================================================
  inline bool CAFHelper::PointInCalo(const TVector3& point)
  {
    //Barrel Radius 278 cm
    //Endcap starts at 364 cm
    bool isInCalo = false;
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    //in the Barrel
    if( r_point > fECALBarrelInnerRadius && r_point < fECALBarrelOuterRadius && std::abs(point.X()) < fECALStartX ) isInCalo = true;
    //in the Endcap
    if( r_point < fECALEndcapOuterRadius && std::abs(point.X()) > fECALStartX && std::abs(point.X()) < fECALEndX ) isInCalo = true;

    return isInCalo;
  }

  //==============================================================================
  inline bool CAFHelper::PointStopBetween(const TVector3& point)
  {
    //Barrel Radius 278 cm
    //Endcap starts at 364 cm
    bool isStopBetween = false;
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    //in the Barrel
    if( r_point < fECALBarrelInnerRadius && r_point > fTPCRadius && std::abs(point.X()) < fTPCLength ) isStopBetween = true;
    //in the Endcap
    if( r_point < fECALEndcapOuterRadius && std::abs(point.X()) > fTPCLength && std::abs(point.X()) < fECALStartX ) isStopBetween = true;

    return isStopBetween;
  }

  //==============================================================================
  inline bool CAFHelper::isThroughCalo(const TVector3& point)
  {
    return !PointInTPC(point) && !PointStopBetween(point) && !PointInCalo(point);
  }

  //==============================================================================
  inline bool CAFHelper::hasDecayedInCalo(const TVector3& point)
  {
    return PointInCalo(point);
  }

  //==============================================================================
  inline bool CAFHelper::isBarrel(const TVector3& point)
  {
    bool isBarrel = false;
    float theta = std::atan(fECALBarrelInnerRadius / std::abs(fECALStartX) ); //angle for barrel/endcap transition
    float r_point = std::sqrt( point.Y()*point.Y() + point.Z()*point.Z() );
    float theta_point = std::atan(r_point / std::abs(point.X()) ); //angle for barrel/endcap transition for the point

    if( theta_point > theta ) isBarrel = true;
    return isBarrel;
  }

  //==============================================================================
  inline bool CAFHelper::isEndcap(const TVector3& point)
  {
    bool isEndcap = false;
    if( !isBarrel(point) ) isEndcap = true;
    return isEndcap;
  }

  //==============================================================================
  class ParamSim : public art::EDAnalyzer {
  public:
    explicit ParamSim(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParamSim(ParamSim const &) = delete;
    ParamSim(ParamSim &&) = delete;
    ParamSim & operator = (ParamSim const &) = delete;
    ParamSim & operator = (ParamSim &&) = delete;

    virtual void beginJob() override;

    // Required functions.
    void analyze(art::Event const & e) override;

  private:
    void ClearVectors();
    void SaveGtruthMCtruth(art::Event const & e);
    void TreatMCParticles(art::Event const & e);
    void TreatTPCVisible(const float ecaltime);
    void TreatTPCNotVisible(const float ecaltime);

    //Helpers
    void FillCommonVariables(const int pdg, const TLorentzVector& momentum, const TLorentzVector& position, const TLorentzVector& positionEnd, const int mctrackid, const int mothertrackid, const int motherpdg, const float ptrue, const float angle, const std::string mcp_process, const std::string mcp_endprocess, const float time);
    void ComputeTrkLength(const simb::MCParticle &mcp);
    void DoRangeCalculation(float &preco, float &angle_reco);
    void DoGluckSternCalculation(float &preco, float &angle_reco);
    void TPCParticleIdentification();
    void TreatNeutrons(const float ecaltime);
    void TreatPhotons(const float ecaltime);
    void TreatOthers(const float ecaltime);
    bool CheckVectorSize();

    //TTree
    TTree *fTree;
    std::unordered_map<int, TH2F*> m_pidinterp;

    //Geometry
    const geo::GeometryCore* fGeo; ///< pointer to the geometry
    CAFHelper* fHelper;
    CLHEP::HepRandomEngine              &fEngine;  ///< random engine

    //fcl parameters
    std::string fGeneratorLabel;
    std::string fGeantLabel; ///< module label for geant4 simulated hits
    bool fCorrect4origin;

    //Fixed parameters
    float fOrigin[3];
    const std::vector<int> neutrinos = {12, 14, 16};
    //gamma, neutron, pi0, k0L, k0S, k0, delta0
    const std::vector<int> pdg_neutral = {22, 2112, 111, 130, 310, 311, 2114};
    //pion, muon, proton, kaon, deuteron, electron
    const std::vector<int> pdg_charged = {211, 13, 2212, 321, 1000010020, 11};
    const float gastpc_len = 2.; // new track length cut in cm based on Thomas' study of low energy protons
    const float gastpc_B = 0.5; // B field strength in Tesla
    const float gastpc_padPitch = 0.1; // 1 mm. Actual pad pitch varies, which is going to be impossible to implement
    const float gastpc_X0 = 1300.; // cm = 13m radiation length
    //Resolution for short tracks //TODO check this numbers!
    const float sigmaP_short = 0.1; //in GeV
    // point resolution
    const float sigma_x = 0.1;
    //as function of KE
    //0 -> 50 MeV ~ 20%
    //> 50 MeV ~ 40%
    const float NeutronECAL_detEff[2] = {0.2, 0.4};
    const float sigmaNeutronECAL_first = 0.11;
    //TODO fraction of rescatters
    // float sigmaNeutronECAL_rescatter = 0.26;
    //ECAL energy resolution sigmaE/E
    const float ECAL_stock = 0.06; //in %
    const float ECAL_const = 0.02;
    TF1 *fRes;
    //ECAL sampling fraction
    // double sampling_frac = 4.32;
    //ECAL nlayers
    const int nLayers = 60;
    //ECAL MIP resolution (based on AHCAL numbers)
    const double ECAL_MIP_Res = 0.23;
    //MIP2GeV conversion factor
    const double MIP2GeV_factor = 0.814 / 1000;
    //float ECAL_pi0_resolution = 0.13; //sigmaE/E in between at rest (17%) and high energy (~few %)
    const float fECALTimeResolution = 1.; // 1 ns time resolution
    TParticlePDG *neutron = TDatabasePDG::Instance()->GetParticle(2112);
    const float neutron_mass = neutron->Mass(); //in GeV

    //CAF variables
    //Event-wise values
    unsigned int fRun, fEvent, fSubRun;
    //Generator values
    std::vector<int> mode, ccnc, ntype, gint, weight, tgtpdg, gt_t, intert, detected;
    std::vector<double> q2, w, y, x, theta, t, mctime, mcnupx, mcnupy, mcnupz, vertx, verty, vertz;
    //MC Particle Values, with motherid added
    unsigned int _nFSP;
    std::vector<int> mctrkid, motherid, pdgmother, truepdg, _MCPStartX, _MCPStartY, _MCPStartZ, _MCPEndX, _MCPEndY, _MCPEndZ;
    std::vector<std::string> _MCProc, _MCEndProc;
    std::vector<double> trkLen, trkLenPerp, truep, truepx, truepy, truepz, _angle;
    //Reco values
    std::vector<int> recopid, recopidecal;
    std::vector<double> prob_arr, anglereco, _preco, erecon, etime;
    //Geometry
    std::vector<unsigned int> isFidStart, isTPCStart, isCaloStart, isInBetweenStart, isThroughCaloStart;
    std::vector<unsigned int> isFidEnd, isTPCEnd, isCaloEnd, isInBetweenEnd, isThroughCaloEnd;
    std::vector<unsigned int> isBarrelStart, isEndcapStart, isBarrelEnd, isEndcapEnd;
  };

  //==============================================================================
  ParamSim::ParamSim(fhicl::ParameterSet const & p)
    : EDAnalyzer(p),
      fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed"))
  {
    fGeo     = gar::providerFrom<geo::GeometryGAr>();

    fGeneratorLabel    = p.get<std::string>("GeneratorLabel","genie");
    fGeantLabel        = p.get<std::string>("GEANTLabel","geant");
    fCorrect4origin    = p.get<bool>("Correct4Origin", false);

    consumes<std::vector<simb::MCTruth> >(fGeneratorLabel);
    consumes<std::vector<simb::GTruth> >(fGeneratorLabel);
    consumes<std::vector<simb::MCParticle> >(fGeantLabel);
  }

  //==============================================================================
  void ParamSim::beginJob()
  {
    fOrigin[0] = fGeo->GetOriginX();
    fOrigin[1] = fGeo->GetOriginY();
    fOrigin[2] = fGeo->GetOriginZ();

    fRes = new TF1("fRes", "TMath::Sqrt ( [0]*[0]/x + [1]*[1] )", 3);
    fRes->FixParameter(0, ECAL_stock);
    fRes->FixParameter(1, ECAL_const);
    fHelper  = new CAFHelper(fGeo, fEngine);

    // read the PID parametrization ntuple from T. Junk
    TString filename = "${DUNE_PARDATA_DIR}/MPD/dedxPID/dedxpidmatrices8kevcm.root";
    TFile pidfile(filename, "READ");

    m_pidinterp.clear();
    char str[11];
    for (int q = 0; q < 501; ++q)
      {
        sprintf(str, "%d", q);
        std::string s = "pidmatrix";
        s.append(str);
        // read the 500 histograms one by one; each histogram is a
        // 6 by 6 matrix of probabilities for a given momentum value
        m_pidinterp.insert( std::make_pair(q, (TH2F*) pidfile.Get(s.c_str())->Clone("pidinterp")) );
      }

    // pidfile.Close();

    art::ServiceHandle<art::TFileService> tfs;
    fTree = tfs->make<TTree>("caf","caf tree");

    //Event number, //Run number, //Sub Run number
    fTree->Branch("Event", &fEvent);
    fTree->Branch("Run", &fRun);
    fTree->Branch("SubRun", &fSubRun);

    //MC Truth info (Generator)
    fTree->Branch("mode", &mode);
    fTree->Branch("q2", &q2);
    fTree->Branch("w", &w);
    fTree->Branch("y", &y);
    fTree->Branch("x", &x);
    fTree->Branch("theta", &theta);
    fTree->Branch("t", &t);
    fTree->Branch("ntype", &ntype);
    fTree->Branch("ccnc", &ccnc);
    fTree->Branch("gint", &gint);
    fTree->Branch("tgtpdg", &tgtpdg);
    fTree->Branch("weight", &weight);
    fTree->Branch("gt_t", &gt_t);
    fTree->Branch("intert", &intert);
    fTree->Branch("mcnupx", &mcnupx);
    fTree->Branch("mcnupy", &mcnupy);
    fTree->Branch("mcnupz", &mcnupz);
    fTree->Branch("vertx", &vertx);
    fTree->Branch("verty", &verty);
    fTree->Branch("vertz", &vertz);

    //Number of final state particle (primaries)
    fTree->Branch("nFSP", &_nFSP);
    //MC Particle info
    fTree->Branch("detected", &detected);
    fTree->Branch("pdgmother", &pdgmother);
    fTree->Branch("MCPTime", &mctime);
    fTree->Branch("MCPStartX", &_MCPStartX);
    fTree->Branch("MCPStartY", &_MCPStartY);
    fTree->Branch("MCPStartZ", &_MCPStartZ);
    fTree->Branch("motherid", &motherid);
    fTree->Branch("mctrkid", &mctrkid);
    fTree->Branch("truepx", &truepx);
    fTree->Branch("truepy", &truepy);
    fTree->Branch("truepz", &truepz);
    fTree->Branch("MCPEndX", &_MCPEndX);
    fTree->Branch("MCPEndY", &_MCPEndY);
    fTree->Branch("MCPEndZ", &_MCPEndZ);
    fTree->Branch("MCProc", &_MCProc);
    fTree->Branch("MCEndProc", &_MCEndProc);
    fTree->Branch("angle", &_angle);
    fTree->Branch("truep", &truep);
    fTree->Branch("truepdg", &truepdg);
    //Reco info
    fTree->Branch("recopid", &recopid);
    fTree->Branch("recopidecal", &recopidecal);
    fTree->Branch("trkLen", &trkLen);
    fTree->Branch("trkLenPerp", &trkLenPerp);
    fTree->Branch("preco", &_preco);
    fTree->Branch("anglereco", &anglereco);
    fTree->Branch("erecon", &erecon);
    fTree->Branch("etime", &etime);
    fTree->Branch("prob_arr", &prob_arr);
    //Geometry
    fTree->Branch("isFidStart", &isFidStart);
    fTree->Branch("isTPCStart", &isTPCStart);
    fTree->Branch("isCaloStart", &isCaloStart);
    fTree->Branch("isInBetweenStart", &isInBetweenStart);
    fTree->Branch("isThroughCaloStart", &isThroughCaloStart);
    fTree->Branch("isBarrelStart", &isBarrelStart);
    fTree->Branch("isEndcapStart", &isEndcapStart);

    fTree->Branch("isFidEnd", &isFidEnd);
    fTree->Branch("isTPCEnd", &isTPCEnd);
    fTree->Branch("isCaloEnd", &isCaloEnd);
    fTree->Branch("isThroughCaloEnd", &isThroughCaloEnd);
    fTree->Branch("isInBetweenEnd", &isInBetweenEnd);
    fTree->Branch("isBarrelEnd", &isBarrelEnd);
    fTree->Branch("isEndcapEnd", &isEndcapEnd);
  }

  //==============================================================================
  void ParamSim::analyze(art::Event const & e)
  {
    ClearVectors();
    fRun    = e.run();
    fSubRun = e.subRun();
    fEvent  = e.id().event();

    SaveGtruthMCtruth(e);
    TreatMCParticles(e);

    //Checks
    if(CheckVectorSize()) {
      fTree->Fill();
    } else {
      std::cerr << "Event " << fEvent << std::endl;
      std::cerr << "Number of FSP " << _nFSP << std::endl;

      std::cerr << "Size of pdgmother " << pdgmother.size() << std::endl;
      std::cerr << "Size of truepdg " << truepdg.size() << std::endl;
      std::cerr << "Size of mctime " << mctime.size() << std::endl;
      std::cerr << "Size of mctrkid " << mctrkid.size() << std::endl;
      std::cerr << "Size of motherid " << motherid.size() << std::endl;
      std::cerr << "Size of _MCPStartX " << _MCPStartX.size() << std::endl;
      std::cerr << "Size of _MCPStartY " << _MCPStartY.size() << std::endl;
      std::cerr << "Size of _MCPStartZ " << _MCPStartZ.size() << std::endl;
      std::cerr << "Size of _MCPEndX " << _MCPEndX.size() << std::endl;
      std::cerr << "Size of _MCPEndY " << _MCPEndY.size() << std::endl;
      std::cerr << "Size of _MCPEndZ " << _MCPEndZ.size() << std::endl;
      std::cerr << "Size of _MCProc " << _MCProc.size() << std::endl;
      std::cerr << "Size of _MCEndProc " << _MCEndProc.size() << std::endl;
      std::cerr << "Size of trkLen " << trkLen.size() << std::endl;
      std::cerr << "Size of trkLenPerp " << trkLenPerp.size() << std::endl;
      std::cerr << "Size of truep " << truep.size() << std::endl;
      std::cerr << "Size of truepx " << truepx.size() << std::endl;
      std::cerr << "Size of truepy " << truepy.size() << std::endl;
      std::cerr << "Size of truepz " << truepz.size() << std::endl;
      std::cerr << "Size of _angle " << _angle.size() << std::endl;

      //Reco values
      std::cerr << "Size of detected " << detected.size() << std::endl;
      std::cerr << "Size of recopid " << recopid.size() << std::endl;
      std::cerr << "Size of recopidecal " << recopidecal.size() << std::endl;
      // std::cerr << "Size of prob_arr " << prob_arr.size() << std::endl;
      std::cerr << "Size of anglereco " << anglereco.size() << std::endl;
      std::cerr << "Size of _preco " << _preco.size() << std::endl;
      std::cerr << "Size of erecon " << erecon.size() << std::endl;
      std::cerr << "Size of etime " << etime.size() << std::endl;

      //Geometry
      std::cerr << "Size of isFidStart " << isFidStart.size() << std::endl;
      std::cerr << "Size of isTPCStart " << isTPCStart.size() << std::endl;
      std::cerr << "Size of isCaloStart " << isCaloStart.size() << std::endl;
      std::cerr << "Size of isThroughCaloStart " << isThroughCaloStart.size() << std::endl;
      std::cerr << "Size of isInBetweenStart " << isInBetweenStart.size() << std::endl;
      std::cerr << "Size of isBarrelStart " << isBarrelStart.size() << std::endl;
      std::cerr << "Size of isEndcapStart " << isEndcapStart.size() << std::endl;

      std::cerr << "Size of isFidEnd " << isFidEnd.size() << std::endl;
      std::cerr << "Size of isTPCEnd " << isTPCEnd.size() << std::endl;
      std::cerr << "Size of isCaloEnd " << isCaloEnd.size() << std::endl;
      std::cerr << "Size of isThroughCaloEnd " << isThroughCaloEnd.size() << std::endl;
      std::cerr << "Size of isInBetweenEnd " << isInBetweenEnd.size() << std::endl;
      std::cerr << "Size of isBarrelEnd " << isBarrelEnd.size() << std::endl;
      std::cerr << "Size of isEndcapEnd " << isEndcapEnd.size() << std::endl;

      std::cerr << "Event with wrong vector sizes... skipped" << std::endl;
    }

  }

  //==============================================================================
  void ParamSim::SaveGtruthMCtruth(art::Event const & e)
  {
    auto MCTHandle = e.getHandle< std::vector<simb::MCTruth> >(fGeneratorLabel);
    if (!MCTHandle) {
      throw cet::exception("ParamSim") << " No simb::MCTruth branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    auto GTHandle = e.getHandle< std::vector<simb::GTruth> >(fGeneratorLabel);
    if (!GTHandle) {
      throw cet::exception("ParamSim") << " No simb::GTruth branch."
                                       << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // save MCTruth info
    for ( auto const& mct : (*MCTHandle) ) {
      if (mct.NeutrinoSet()) {
        simb::MCNeutrino nuw = mct.GetNeutrino();
        ccnc.push_back(nuw.CCNC());
        ntype.push_back(nuw.Nu().PdgCode());
        q2.push_back(nuw.QSqr());
        w.push_back(nuw.W());
        y.push_back(nuw.Y());
        x.push_back(nuw.X());
        theta.push_back(nuw.Theta());
        mode.push_back(nuw.Mode());
        intert.push_back(nuw.InteractionType());
        if(fCorrect4origin){
          vertx.push_back(nuw.Nu().EndX() - fOrigin[0]);
          verty.push_back(nuw.Nu().EndY() - fOrigin[1]);
          vertz.push_back(nuw.Nu().EndZ() - fOrigin[2]);
        } else {
          vertx.push_back(nuw.Nu().EndX());
          verty.push_back(nuw.Nu().EndY());
          vertz.push_back(nuw.Nu().EndZ());
        }
        mcnupx.push_back(nuw.Nu().Px());
        mcnupy.push_back(nuw.Nu().Py());
        mcnupz.push_back(nuw.Nu().Pz());
      }  // end MC info from MCTruth
    }

    // save GTruth info
    for ( auto const& gt : (*GTHandle) ) {
      gint.push_back(gt.fGint);
      tgtpdg.push_back(gt.ftgtPDG);
      weight.push_back(gt.fweight);
      gt_t.push_back(gt.fgT);
    }
  }

  //==============================================================================
  void ParamSim::TreatMCParticles(art::Event const & e)
  {
    auto MCPHandle = e.getHandle< std::vector<simb::MCParticle> >(fGeantLabel);
    if (!MCPHandle) {
      throw cet::exception("anatree") << " No simb::MCParticle branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    std::unordered_map<int, int> TrackIdToIndex;
    int index = 0;
    for ( auto const& mcp : (*MCPHandle) ) {
      int TrackId = mcp.TrackId();
      TrackIdToIndex[TrackId] = index++;
    }

    //Loop over the mcp
    for ( auto const& mcp : (*MCPHandle) ) {

      const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
      const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
      if (definition == nullptr) continue;

      const std::string mcp_process = mcp.Process();
      const std::string mcp_endprocess = mcp.EndProcess();
      const int mctrackid = mcp.TrackId();
      const int mothertrackid = mcp.Mother();
      const int pdg = mcp.PdgCode();
      const TLorentzVector& position = mcp.Position(0);
      const TVector3 spoint(position.X() - fOrigin[0], position.Y() - fOrigin[1], position.Z() - fOrigin[2]);
      const TLorentzVector& positionEnd = mcp.EndPosition();
      const TVector3 epoint(positionEnd.X() - fOrigin[0], positionEnd.Y() - fOrigin[1], positionEnd.Z() - fOrigin[2]);
      const TLorentzVector& momentum = mcp.Momentum(0);
      const TVector3 mom(momentum.X(), momentum.Y(), momentum.Z());
      const float ptrue = mom.Mag();
      const float angle  = atan(mom.X() / mom.Z());
      const float time = mcp.T();
      int motherpdg = 0;

      //Find out mother pdg
      if(TrackIdToIndex.find(mothertrackid) != TrackIdToIndex.end()) {
        motherpdg = (*MCPHandle).at(TrackIdToIndex[mothertrackid]).PdgCode();
      }

      //need to ignore neutrals for this - put the value to 0
      auto result = std::find(pdg_neutral.begin(), pdg_neutral.end(), abs(pdg));
      bool isNeutral = (result != pdg_neutral.end()) ? true : false;

      if( isNeutral )
        {
          trkLen.push_back(-1);
          trkLenPerp.push_back(-1);
        }
      else {
        ComputeTrkLength(mcp);
      }

      float ecaltime = fHelper->GaussianSmearing(time, fECALTimeResolution);

      //Store common mcp properties
      FillCommonVariables(pdg, momentum, position, positionEnd, mctrackid, mothertrackid, motherpdg, ptrue, angle, mcp_process, mcp_endprocess, time);

      //Treat visible particles in the TPC
      if( trkLen.at(_nFSP) > gastpc_len ) {
        TreatTPCVisible(ecaltime);
      }
      else
        {
          TreatTPCNotVisible(ecaltime);
        }// end trkLen.at(_nFSP) < gastpc_len

      _nFSP++;
    }//end loop mcp
  }

  //==============================================================================
  void ParamSim::FillCommonVariables(const int pdg, const TLorentzVector& momentum, const TLorentzVector& position, const TLorentzVector& positionEnd, const int mctrackid, const int mothertrackid, const int motherpdg, const float ptrue, const float angle, const std::string mcp_process, const std::string mcp_endprocess, const float time)
  {
    const TVector3 spoint(position.X() - fOrigin[0], position.Y() - fOrigin[1], position.Z() - fOrigin[2]);
    const TVector3 epoint(positionEnd.X() - fOrigin[0], positionEnd.Y() - fOrigin[1], positionEnd.Z() - fOrigin[2]);

    truepdg.push_back(pdg);
    truepx.push_back(momentum.X());
    truepy.push_back(momentum.Y());
    truepz.push_back(momentum.Z());

    if(fCorrect4origin){
      _MCPStartX.push_back(position.X() - fOrigin[0]);
      _MCPStartY.push_back(position.Y() - fOrigin[1]);
      _MCPStartZ.push_back(position.Z() - fOrigin[2]);
      _MCPEndX.push_back(positionEnd.X() - fOrigin[0]);
      _MCPEndY.push_back(positionEnd.Y() - fOrigin[1]);
      _MCPEndZ.push_back(positionEnd.Z() - fOrigin[2]);
    } else {
      _MCPStartX.push_back(position.X());
      _MCPStartY.push_back(position.Y());
      _MCPStartZ.push_back(position.Z());
      _MCPEndX.push_back(positionEnd.X());
      _MCPEndY.push_back(positionEnd.Y());
      _MCPEndZ.push_back(positionEnd.Z());
    }

    // save the true momentum
    truep.push_back(ptrue);
    // save the true angle
    _angle.push_back(angle);
    //Save MC process
    _MCProc.push_back(mcp_process);
    _MCEndProc.push_back(mcp_endprocess);
    mctime.push_back(time);
    mctrkid.push_back(mctrackid);
    motherid.push_back(mothertrackid);
    pdgmother.push_back(motherpdg);

    isFidStart.push_back(fHelper->PointInFiducial(spoint));
    isTPCStart.push_back(fHelper->PointInTPC(spoint));
    isCaloStart.push_back(fHelper->PointInCalo(spoint));
    isThroughCaloStart.push_back(fHelper->isThroughCalo(spoint));
    isInBetweenStart.push_back(fHelper->PointStopBetween(spoint));
    isBarrelStart.push_back(fHelper->isBarrel(spoint));
    isEndcapStart.push_back(fHelper->isEndcap(spoint));
    //Check where endpoint of mcp is
    isFidEnd.push_back(fHelper->PointInFiducial(epoint));
    isTPCEnd.push_back(fHelper->PointInTPC(epoint));
    isCaloEnd.push_back(fHelper->PointInCalo(epoint));
    isThroughCaloEnd.push_back(fHelper->isThroughCalo(epoint));
    isInBetweenEnd.push_back(fHelper->PointStopBetween(epoint));
    isBarrelEnd.push_back(fHelper->isBarrel(epoint));
    isEndcapEnd.push_back(fHelper->isEndcap(epoint));
  }

  //==============================================================================
  void ParamSim::ComputeTrkLength(const simb::MCParticle &mcp)
  {
    //start track length
    //***************************************************************************************************************/

    // calculate the total and the transverse track lengths and restrict the
    // tracklength to be above the gas TPC track length threshold
    double tracklen = 0.;
    double tracklen_perp = 0.;

    //CAREFUL No offset for the trajectory points (origin for them is the TPC?)??????
    //TODO check if the mcp point is within the TPC volume! Skip for mcp in the ECAL (showers)
    //TODO Link showers to original mcp?
    for(size_t itraj = 1; itraj < mcp.Trajectory().size(); itraj++)
      {
        float xTraj = mcp.Trajectory().X(itraj);
        float yTraj = mcp.Trajectory().Y(itraj);
        float zTraj = mcp.Trajectory().Z(itraj);

        //Traj point+1
        TVector3 point(xTraj - fOrigin[0], yTraj - fOrigin[1], zTraj - fOrigin[2]);

        //point is not in the TPC anymore - stop traj loop
        if(not fHelper->PointInTPC(point))
          continue;

        // find the length of the track by getting the distance between each hit
        TVector3 diff(xTraj - mcp.Trajectory().X(itraj - 1), yTraj - mcp.Trajectory().Y(itraj - 1), zTraj - mcp.Trajectory().Z(itraj - 1));
        // perp length
        TVector2 tracklen_perp_vec(zTraj - mcp.Trajectory().Z(itraj - 1), yTraj - mcp.Trajectory().Y(itraj - 1));
        // Summing up
        tracklen += diff.Mag();
        tracklen_perp += tracklen_perp_vec.Mod();
      }

    trkLen.push_back(tracklen);
    trkLenPerp.push_back(tracklen_perp);
  }

  //==============================================================================
  void ParamSim::DoRangeCalculation(float &preco, float &angle_reco)
  {
    // calculate number of trackpoints
    float nHits = round (trkLen.at(_nFSP) / gastpc_padPitch);
    // angular resolution first term
    float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / trkLen.at(_nFSP)*trkLen.at(_nFSP)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
    // scattering term in Gluckstern formula
    float sigma_angle_2 = (0.015*0.015 / (3. * truep.at(_nFSP) * truep.at(_nFSP))) * (trkLen.at(_nFSP)/gastpc_X0);
    // angular resolution from the two terms above
    float sigma_angle_short = sqrt(sigma_angle_1 + sigma_angle_2);

    //reconstructed angle
    angle_reco = fHelper->GaussianSmearing(_angle.at(_nFSP), sigma_angle_short);
    //reconstructed momentum
    preco = fHelper->GaussianSmearing( truep.at(_nFSP), sigmaP_short );
  }

  //==============================================================================
  void ParamSim::DoGluckSternCalculation(float &preco, float &angle_reco)
  {
    //Case where the endpoint is not in the TPC, should be able to use the Gluckstern formula
    // calculate number of trackpoints
    float nHits = round (trkLen.at(_nFSP) / gastpc_padPitch);
    // measurement term in Gluckstern formula
    float fracSig_meas = sqrt(720./(nHits+4)) * ((0.01*gastpc_padPitch*truep.at(_nFSP)) / (0.3 * gastpc_B * 0.0001 *trkLenPerp.at(_nFSP)*trkLenPerp.at(_nFSP)));
    // multiple Coulomb scattering term in Gluckstern formula
    float fracSig_MCS = (0.052*sqrt(1.43)) / (gastpc_B * sqrt(gastpc_X0*trkLenPerp.at(_nFSP)*0.0001));
    // momentum resoltion from the two terms above
    float sigmaP = truep.at(_nFSP) * sqrt( fracSig_meas*fracSig_meas + fracSig_MCS*fracSig_MCS );
    // now Gaussian smear the true momentum using the momentum resolution
    preco = fHelper->GaussianSmearing( truep.at(_nFSP), sigmaP );

    // measurement term in the Gluckstern formula for calculating the
    // angular resolution
    float sigma_angle_1 = ((sigma_x * sigma_x * 0.0001) / trkLen.at(_nFSP)*trkLen.at(_nFSP)*0.0001) * (12*(nHits-1))/(nHits*(nHits+1));
    // scattering term in Gluckstern formula
    float sigma_angle_2 = (0.015*0.015 / (3. * truep.at(_nFSP) * truep.at(_nFSP))) * (trkLen.at(_nFSP)/gastpc_X0);
    // angular resolution from the two terms above
    float sigma_angle = sqrt(sigma_angle_1 + sigma_angle_2);
    // now Gaussian smear the true angle using the angular resolution
    angle_reco = fHelper->GaussianSmearing(_angle.at(_nFSP), sigma_angle);
  }

  //==============================================================================
  void ParamSim::TreatTPCVisible(const float ecaltime)
  {
    //start tpc
    //***************************************************************************************************************/
    int pdg = truepdg.at(_nFSP);

    if ( std::find(pdg_charged.begin(), pdg_charged.end(), std::abs(pdg)) != pdg_charged.end() )
      {
        TVector3 epoint(_MCPEndX.at(_nFSP), _MCPEndY.at(_nFSP), _MCPEndZ.at(_nFSP));
        float ptrue = truep.at(_nFSP);

        //Use range instead of Gluckstern for stopping tracks
        //TODO is that correct? What if it is a scatter in the TPC? Need to check if daughter is same particle
        float preco = 0.;
        float angle_reco = 0.;

        //Case for range, the end point of the mcp is in the TPC, does not reach the ecal
        if( fHelper->PointInTPC(epoint) )
          {
            DoRangeCalculation(preco, angle_reco);

            if(preco > 0)
              _preco.push_back(preco);
            else
              _preco.push_back(-1);
            anglereco.push_back(angle_reco);
            erecon.push_back(-1);
            recopidecal.push_back(-1);
            etime.push_back(-1);
            detected.push_back(-1);
          }
        else
          {
            DoGluckSternCalculation(preco, angle_reco);

            // save reconstructed momentum and angle to cafanatree
            if(preco > 0)
              _preco.push_back(preco);
            else
              _preco.push_back(-1);
            anglereco.push_back(angle_reco);

            //Reaches the ECAL and stops there
            if( fHelper->PointInCalo(epoint) )
              {
                //Need energy measurement in ecal
                TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(pdg));
                if(nullptr == part)
                  {
                    //deuteron
                    if( pdg == 1000010020 )
                      {
                        float mass = 1.8756;//in GeV mass deuteron
                        float etrue = std::sqrt(ptrue*ptrue + mass*mass) - mass;
                        float ECAL_resolution = fRes->Eval(etrue)*etrue;
                        float ereco = fHelper->GaussianSmearing(etrue, ECAL_resolution);
                        erecon.push_back((ereco > 0) ? ereco : 0.);
                        recopidecal.push_back(-1);
                        detected.push_back(1);
                        etime.push_back(ecaltime);
                      }
                    else
                      {
                        erecon.push_back(-1);
                        recopidecal.push_back(-1);
                        detected.push_back(0);
                        etime.push_back(-1);
                      }
                  }
                else
                  {
                    //by default should be tagged as an electron as it has a track,
                    //otherwise tag as gamma if not track -> need mis association rate, and use dE/dX in Scintillator?
                    //separation between e and mu/pi should be around 100%
                    //separation mu/pi -> based on Chris study with only the ECAL (no Muon ID detector)
                    //separation with p and mu/pi/e ?? high energy -> confusion with mu/pi, low energy confusion with e
                    //using E/p to ID?
                    float mass = part->Mass();//in GeV
                    float etrue = std::sqrt(ptrue*ptrue + mass*mass) - mass;
                    float ECAL_resolution = fRes->Eval(etrue)*etrue;
                    float ereco = fHelper->GaussianSmearing(etrue, ECAL_resolution);
                    erecon.push_back((ereco > 0) ? ereco : 0.);
                    detected.push_back(1);
                    etime.push_back(ecaltime);

                    //Electron
                    if( abs(pdg) == 11 ){
                      recopidecal.push_back(11);
                    }
                    else if( abs(pdg) == 13 || abs(pdg) == 211 )
                      {
                        //Muons and Pions
                        //ptrue < 480 MeV/c 100% separation
                        //80% from 480 to 750
                        //90% up to 750 to 900
                        //95% over 900
                        float random_number = fHelper->GetRamdomNumber();

                        if(ptrue < 0.48) {
                          recopidecal.push_back(abs(pdg));//100% efficiency by range
                        }
                        else if(ptrue >= 0.48 && ptrue < 0.75)
                          {
                            //case muon
                            if(abs(pdg) == 13)
                              {
                                if(random_number > (1 - 0.8)) {
                                  recopidecal.push_back(13);
                                }
                                else{
                                  recopidecal.push_back(211);
                                }
                              }
                            //case pion
                            if(abs(pdg) == 211)
                              {
                                if(random_number > (1 - 0.8)) {
                                  recopidecal.push_back(211);
                                }
                                else{
                                  recopidecal.push_back(13);
                                }
                              }
                          }
                        else if(ptrue >= 0.75 && ptrue < 0.9)
                          {
                            //case muon
                            if(abs(pdg) == 13){
                              if(random_number > (1 - 0.9)) {
                                recopidecal.push_back(13);
                              }
                              else{
                                recopidecal.push_back(211);
                              }
                            }
                            //case pion
                            if(abs(pdg) == 211) {
                              if(random_number > (1 - 0.9)) {
                                recopidecal.push_back(211);
                              }
                              else{
                                recopidecal.push_back(13);
                              }
                            }
                          }
                        else
                          {
                            //case muon
                            if(abs(pdg) == 13){
                              if(random_number > (1 - 0.95)) {
                                recopidecal.push_back(13);
                              }
                              else{
                                recopidecal.push_back(211);
                              }
                            }
                            //case pion
                            if(abs(pdg) == 211){
                              if(random_number > (1 - 0.95)) {
                                recopidecal.push_back(211);
                              }
                              else{
                                recopidecal.push_back(13);
                              }
                            }
                          }
                      }
                    else if( abs(pdg) == 2212 )
                      {
                        recopidecal.push_back(2212);//TODO for p/pi separation
                      }
                    else {
                      recopidecal.push_back(-1);
                    }
                  }
              }
            else if( fHelper->isThroughCalo(epoint) )
              {
                //Case the endpoint is outside the CALO -> it went through the ECAL (mu/pi/p possible)
                //the ECAL will see 60 MIPs on average
                double Evis = (double)nLayers; //in MIP
                //Smearing to account for Gaussian detector noise (Landau negligible)
                Evis = fHelper->GaussianSmearing(Evis, ECAL_MIP_Res);
                //1 MIP = 0.814 MeV
                double Erec = Evis * MIP2GeV_factor;
                erecon.push_back((Erec > 0) ? Erec : 0.);
                etime.push_back(ecaltime);
                detected.push_back(1);

                //Muon/Pions/Protons are reco as Muons (without MuID detector)
                if( abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 ) {
                  recopidecal.push_back(13);
                }
                else{
                  recopidecal.push_back(-1);
                }
              }
            else
              {
                //Does not reach the ECAL???
                erecon.push_back(-1);
                recopidecal.push_back(-1);
                etime.push_back(-1);
                detected.push_back(0);
              }
          } //end endpoint is not in TPC

        TPCParticleIdentification();

      } // end is charged mcp
    else
      {
        //not in the pdglist of particles but visible in TPC?
        auto found = std::find(pdg_charged.begin(), pdg_charged.end(), abs(pdg));
        if(found == pdg_charged.end())
          {
            detected.push_back(0);
            etime.push_back(-1);
            erecon.push_back(-1);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            recopid.push_back(-1);
            recopidecal.push_back(-1);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
      }

    //end tpc
    //***************************************************************************************************************/
  }

  void ParamSim::TPCParticleIdentification()
  {
    //start pid
    //***************************************************************************************************************/

    int pdg = truepdg.at(_nFSP);

    for (int pidm = 0; pidm < 6; ++pidm)
      {
        if ( abs(pdg) == pdg_charged.at(pidm) )
          {
            float p = _preco.at(_nFSP);
            std::vector<double> vec;
            std::vector<std::string> pnamelist     = {"#pi", "#mu", "p", "K", "d", "e"};
            std::vector<std::string> recopnamelist = {"#pi", "#mu", "p", "K", "d", "e"};

            int qclosest = 0;
            float dist = 100000000.;

            for (int q = 0; q < 501; ++q)
              {
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

                //loop over the columns (true pid)
            std::vector< std::pair<float, std::string> > v_prob;
            //get true particle name
            std::string trueparticlename = m_pidinterp[qclosest]->GetXaxis()->GetBinLabel(pidm+1);
            // std::cout << trueparticlename << std::endl;
            if ( trueparticlename == pnamelist[pidm] )
              {
                //loop over the rows (reco pid)
                for (int pidr = 0; pidr < 6; ++pidr)
                  {
                    std::string recoparticlename = m_pidinterp[qclosest]->GetYaxis()->GetBinLabel(pidr+1);
                    // std::cout << recoparticlename << std::endl;
                    if (recoparticlename == recopnamelist[pidr])
                      {
                        float prob = m_pidinterp[qclosest]->GetBinContent(pidm+1,pidr+1);
                        prob_arr.push_back(prob);
                        //Need to check random number value and prob value then associate the recopdg to the reco prob
                        v_prob.push_back( std::make_pair(prob, recoparticlename) );
                      }
                  }

                int pid = -1;
                if(v_prob.size() > 1)
                  {
                    //Order the vector of prob
                    std::sort(v_prob.begin(), v_prob.end());
                    //Throw a random number between 0 and 1
                    float random_number = fHelper->GetRamdomNumber();
                    //Make cumulative sum to get the range
                    std::partial_sum(v_prob.begin(), v_prob.end(), v_prob.begin(), [](const std::pair<float, std::string>& _x, const std::pair<float, std::string>& _y){return std::pair<float, std::string>(_x.first + _y.first, _y.second);});
                    for(size_t ivec = 0; ivec < v_prob.size()-1; ivec++)
                      {
                        if( random_number < v_prob.at(ivec+1).first && random_number >= v_prob.at(ivec).first )
                          {
                            pid = pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(ivec+1).second) ) );
                          }
                      }
                  }
                else
                  {
                    pid = pdg_charged.at( std::distance( recopnamelist.begin(), std::find(recopnamelist.begin(), recopnamelist.end(), v_prob.at(0).second) ) );
                  }

                recopid.push_back( pid );
              } // closes the if statement
          } // end if charged in pdg list
      } // end loop pidm

        //end pid
        //***************************************************************************************************************/
  }

  //==============================================================================
  void ParamSim::TreatTPCNotVisible(const float ecaltime)
  {
    int pdg = truepdg.at(_nFSP);

    //Case neutrinos
    if(std::find(neutrinos.begin(), neutrinos.end(), std::abs(pdg)) != neutrinos.end())
      {
        detected.push_back(0);
        etime.push_back(0.);
        erecon.push_back(0.);
        recopidecal.push_back(-1);
        _preco.push_back(-1);
        anglereco.push_back(-1);
        recopid.push_back(-1);
        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }
    else if( std::abs(pdg) == 2112 )
      {
        TreatNeutrons(ecaltime);
      }
    else if(std::abs(pdg) == 111)
      {
        //Pi0 case
        erecon.push_back(-1);
        recopid.push_back(-1);
        detected.push_back(0);
        recopidecal.push_back(-1);
        etime.push_back(-1);
        _preco.push_back(-1);
        anglereco.push_back(-1);

        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }
    else if(std::abs(pdg) == 22)
      {
        TreatPhotons(ecaltime);
      }
    else
      {
        TreatOthers(ecaltime);
      }//end case charged but not visible in TPC
  }

  //==============================================================================
  void ParamSim::TreatNeutrons(const float ecaltime)
  {
    TVector3 epoint(_MCPEndX.at(_nFSP), _MCPEndY.at(_nFSP), _MCPEndZ.at(_nFSP));
    //start neutrons
    //***************************************************************************************************************/

    if(fHelper->PointInCalo(epoint)) //needs to stop in the ECAL
      {
        //check if it can be detected by the ECAL
        //Assumes 40% efficiency to detect
        float random_number = fHelper->GetRamdomNumber();
        float true_KE = std::sqrt(truep.at(_nFSP)*truep.at(_nFSP) + neutron_mass*neutron_mass) - neutron_mass;
        // float true_KE = ptrue*ptrue / (2*neutron_mass); // in GeV
        int index = (true_KE >= 0.05) ? 1 : 0;

        if(random_number > (1 - NeutronECAL_detEff[index]) && true_KE > 0.003)//Threshold of 3 MeV
          {
            //TODO random is first interaction or rescatter and smear accordingly to Chris's study
            //Detected in the ECAL
            recopid.push_back(-1); //reco pid set to 0?
            detected.push_back(1);
            float eres = sigmaNeutronECAL_first * true_KE;
            float ereco = fHelper->GaussianSmearing( true_KE, eres );
            erecon.push_back(ereco > 0 ? ereco : 0.);
            etime.push_back(ecaltime);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            recopidecal.push_back(2112);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
        else
          {
            //neutron not detected
            detected.push_back(0);
            recopid.push_back(-1);
            erecon.push_back(-1);
            etime.push_back(-1);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            recopidecal.push_back(-1);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
      } //endpoint is in ECAL
    else
      {
        //Endpoint is not in calo (TPC/isInBetween or outside Calo)
        detected.push_back(0);
        recopid.push_back(-1);
        erecon.push_back(-1);
        etime.push_back(-1);
        _preco.push_back(-1);
        anglereco.push_back(-1);
        recopidecal.push_back(-1);
        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }

    //End neutrons
    //***************************************************************************************************************/
  }

  //==============================================================================
  void ParamSim::TreatPhotons(const float ecaltime)
  {
    TVector3 epoint(_MCPEndX.at(_nFSP), _MCPEndY.at(_nFSP), _MCPEndZ.at(_nFSP));
    //start gammas
    //***************************************************************************************************************/

    if( pdgmother.at(_nFSP) != 111 )
      {
        //Endpoint is in the ECAL
        if(fHelper->PointInCalo(epoint))
          {
            //if they hit the ECAL and smear their energy
            float ECAL_resolution = fRes->Eval(truep.at(_nFSP))*truep.at(_nFSP);
            float ereco = fHelper->GaussianSmearing(truep.at(_nFSP), ECAL_resolution);
            erecon.push_back( (ereco > 0) ? ereco : 0. );
            recopid.push_back(-1);
            detected.push_back(1);
            etime.push_back(ecaltime);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            //reach the ECAL, should be tagged as gamma
            recopidecal.push_back(22);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
        else if(fHelper->PointInTPC(epoint) || fHelper->PointStopBetween(epoint) || fHelper->isThroughCalo(epoint))
          {
            erecon.push_back(-1);
            recopid.push_back(-1);
            detected.push_back(0);
            etime.push_back(-1);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            //converted so not seen in ECAL
            recopidecal.push_back(-1);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
      }
    else
      {
        //From a pi0
        if(fHelper->PointInCalo(epoint))
          {
            //if they hit the ECAL and smear their energy
            float ECAL_resolution = fRes->Eval(truep.at(_nFSP))*truep.at(_nFSP);
            float ereco = fHelper->GaussianSmearing(truep.at(_nFSP), ECAL_resolution);
            erecon.push_back((ereco > 0) ? ereco : 0.);
            recopid.push_back(-1);
            detected.push_back(1);
            etime.push_back(ecaltime);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            //reaches the ecal
            recopidecal.push_back(22);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
        else if(fHelper->PointInTPC(epoint) || fHelper->PointStopBetween(epoint) || fHelper->isThroughCalo(epoint))
          {
            //from pi0 and converted in TPC or stopped between TPC/ECAL
            erecon.push_back(-1);
            recopid.push_back(-1);
            detected.push_back(0);
            etime.push_back(-1);
            _preco.push_back(-1);
            anglereco.push_back(-1);
            //converted not seen by ecal
            recopidecal.push_back(-1);
            for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
          }
      }

    //end gammas
    //***************************************************************************************************************/
  }

  //==============================================================================
  void ParamSim::TreatOthers(const float ecaltime)
  {
    int pdg = truepdg.at(_nFSP);
    TVector3 epoint(_MCPEndX.at(_nFSP), _MCPEndY.at(_nFSP), _MCPEndZ.at(_nFSP));
    //Case for particles that stop or go through ECAL (problematic particles with no track length????)
    //Not visible in the TPC and not neutron or gamma or pi0 (otherwise it has been already done above)
    if(fHelper->PointInCalo(epoint))
      {
        detected.push_back(1);
        etime.push_back(ecaltime);

        TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(abs(pdg));
        float mass = 0.;
        if(nullptr != part)
          mass = part->Mass();//in GeV

        float etrue = std::sqrt(truep.at(_nFSP)*truep.at(_nFSP) + mass*mass) - mass;
        float ECAL_resolution = fRes->Eval(etrue)*etrue;
        float ereco = fHelper->GaussianSmearing(etrue, ECAL_resolution);
        erecon.push_back((ereco > 0) ? ereco : 0.);

        //Electron
        if( abs(pdg) == 11 ){
          recopidecal.push_back(11);
        }
        else if( abs(pdg) == 13 || abs(pdg) == 211 )
          {
            //Muons and Pions
            //ptrue < 480 MeV/c 100% separation
            //80% from 480 to 750
            //90% up to 750 to 900
            //95% over 900
            float random_number = fHelper->GetRamdomNumber();

            if(truep.at(_nFSP) < 0.48) {
              recopidecal.push_back(abs(pdg));//100% efficiency by range
            }
            else if(truep.at(_nFSP) >= 0.48 && truep.at(_nFSP) < 0.75)
              {
                //case muon
                if(abs(pdg) == 13)
                  {
                    if(random_number > (1 - 0.8)) {
                      recopidecal.push_back(13);
                    }
                    else{
                      recopidecal.push_back(211);
                    }
                  }

                //case pion
                if(abs(pdg) == 211)
                  {
                    if(random_number > (1 - 0.8)) {
                      recopidecal.push_back(211);
                    }
                    else{
                      recopidecal.push_back(13);
                    }
                  }
              }
            else if(truep.at(_nFSP) >= 0.75 && truep.at(_nFSP) < 0.9)
              {
                //case muon
                if(abs(pdg) == 13){
                  if(random_number > (1 - 0.9)) {
                    recopidecal.push_back(13);
                  }
                  else{
                    recopidecal.push_back(211);
                  }
                }
                //case pion
                if(abs(pdg) == 211) {
                  if(random_number > (1 - 0.9)) {
                    recopidecal.push_back(211);
                  }
                  else{
                    recopidecal.push_back(13);
                  }
                }
              }
            else
              {
                //case muon
                if(abs(pdg) == 13){
                  if(random_number > (1 - 0.95)) {
                    recopidecal.push_back(13);
                  }
                  else{
                    recopidecal.push_back(211);
                  }
                }
                //case pion
                if(abs(pdg) == 211){
                  if(random_number > (1 - 0.95)) {
                    recopidecal.push_back(211);
                  }
                  else{
                    recopidecal.push_back(13);
                  }
                }
              }
          }
        else if( abs(pdg) == 2212 )
          {
            recopidecal.push_back(2212);//TODO for p/pi separation
          }
        else {
          recopidecal.push_back(-1);
        }

        _preco.push_back(-1);
        anglereco.push_back(-1);
        recopid.push_back(-1);
        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }
    else if (fHelper->isThroughCalo(epoint))
      {
        //Case the endpoint is outside the CALO -> it went through the ECAL (mu/pi/p possible)
        //the ECAL will see 60 MIPs on average
        double Evis = (double)nLayers; //in MIP
        //Smearing to account for Gaussian detector noise (Landau negligible)
        Evis = fHelper->GaussianSmearing(Evis, ECAL_MIP_Res);
        //1 MIP = 0.814 MeV
        double Erec = Evis * MIP2GeV_factor;
        erecon.push_back((Erec > 0) ? Erec : 0.);
        etime.push_back(ecaltime);
        detected.push_back(1);

        //Muon/Pions/Protons are reco as Muons (without MuID detector)
        if( abs(pdg) == 13 || abs(pdg) == 211 || abs(pdg) == 2212 ) {
          recopidecal.push_back(13);
        }
        else {
          recopidecal.push_back(-1);
        }

        _preco.push_back(-1);
        anglereco.push_back(-1);
        recopid.push_back(-1);
        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }
    else if(fHelper->PointInTPC(epoint) || fHelper->PointStopBetween(epoint))
      {
        detected.push_back(0);
        etime.push_back(-1);
        erecon.push_back(-1);
        _preco.push_back(-1);
        anglereco.push_back(-1);
        recopid.push_back(-1);
        recopidecal.push_back(-1);
        for (int pidr = 0; pidr < 6; ++pidr) prob_arr.push_back(-1);
      }
  }

  //==============================================================================
  void ParamSim::ClearVectors()
  {
    //Generator values
    mode.clear();
    ccnc.clear();
    ntype.clear();
    gint.clear();
    weight.clear();
    tgtpdg.clear();
    gt_t.clear();
    intert.clear();
    q2.clear();
    w.clear();
    y.clear();
    x.clear();
    theta.clear();
    t.clear();
    mcnupx.clear();
    mcnupy.clear();
    mcnupz.clear();
    vertx.clear();
    verty.clear();
    vertz.clear();

    //MC Particle values
    _nFSP = 0;
    detected.clear();
    pdgmother.clear();
    truepdg.clear();
    mctime.clear();
    mctrkid.clear();
    motherid.clear();
    _MCPStartX.clear();
    _MCPStartY.clear();
    _MCPStartZ.clear();
    _MCPEndX.clear();
    _MCPEndY.clear();
    _MCPEndZ.clear();
    _MCProc.clear();
    _MCEndProc.clear();
    trkLen.clear();
    trkLenPerp.clear();
    truep.clear();
    truepx.clear();
    truepy.clear();
    truepz.clear();
    _angle.clear();

    //Reco values
    recopid.clear();
    recopidecal.clear();
    prob_arr.clear();
    anglereco.clear();
    _preco.clear();
    erecon.clear();
    etime.clear();

    //Geometry
    isFidStart.clear();
    isTPCStart.clear();
    isCaloStart.clear();
    isThroughCaloStart.clear();
    isInBetweenStart.clear();
    isBarrelStart.clear();
    isEndcapStart.clear();

    isFidEnd.clear();
    isTPCEnd.clear();
    isCaloEnd.clear();
    isThroughCaloEnd.clear();
    isInBetweenEnd.clear();
    isBarrelEnd.clear();
    isEndcapEnd.clear();
  }

  //==============================================================================
  bool ParamSim::CheckVectorSize()
  {
    bool isOK = true;

    if(_nFSP != pdgmother.size()) isOK = false;
    if(_nFSP != truepdg.size()) isOK = false;
    if(_nFSP != mctime.size()) isOK = false;
    if(_nFSP != mctrkid.size()) isOK = false;
    if(_nFSP != motherid.size()) isOK = false;
    if(_nFSP != _MCPStartX.size()) isOK = false;
    if(_nFSP != _MCPStartY.size()) isOK = false;
    if(_nFSP != _MCPStartZ.size()) isOK = false;
    if(_nFSP != _MCPEndX.size()) isOK = false;
    if(_nFSP != _MCPEndY.size()) isOK = false;
    if(_nFSP != _MCPEndZ.size()) isOK = false;
    if(_nFSP != _MCProc.size()) isOK = false;
    if(_nFSP != _MCEndProc.size()) isOK = false;
    if(_nFSP != trkLen.size()) isOK = false;
    if(_nFSP != trkLenPerp.size()) isOK = false;
    if(_nFSP != truep.size()) isOK = false;
    if(_nFSP != truepx.size()) isOK = false;
    if(_nFSP != truepy.size()) isOK = false;
    if(_nFSP != truepz.size()) isOK = false;
    if(_nFSP != _angle.size()) isOK = false;

    //Reco values
    if(_nFSP != recopid.size()) isOK = false;
    if(_nFSP != detected.size()) isOK = false;
    if(_nFSP != recopidecal.size()) isOK = false;
    // if(_nFSP != prob_arr.size()) isOK = false;
    if(_nFSP != anglereco.size()) isOK = false;
    if(_nFSP != _preco.size()) isOK = false;
    if(_nFSP != erecon.size()) isOK = false;
    if(_nFSP != etime.size()) isOK = false;

    //Geometry
    if(_nFSP != isFidStart.size()) isOK = false;
    if(_nFSP != isTPCStart.size()) isOK = false;
    if(_nFSP != isCaloStart.size()) isOK = false;
    if(_nFSP != isThroughCaloStart.size()) isOK = false;
    if(_nFSP != isInBetweenStart.size()) isOK = false;
    if(_nFSP != isBarrelStart.size()) isOK = false;
    if(_nFSP != isEndcapStart.size()) isOK = false;

    if(_nFSP != isFidEnd.size()) isOK = false;
    if(_nFSP != isTPCEnd.size()) isOK = false;
    if(_nFSP != isCaloEnd.size()) isOK = false;
    if(_nFSP != isThroughCaloEnd.size()) isOK = false;
    if(_nFSP != isInBetweenEnd.size()) isOK = false;
    if(_nFSP != isBarrelEnd.size()) isOK = false;
    if(_nFSP != isEndcapEnd.size()) isOK = false;

    return isOK;
  }

} //end namespace gar

namespace gar {
  DEFINE_ART_MODULE(ParamSim)
}
