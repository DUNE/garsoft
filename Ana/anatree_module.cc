////////////////////////////////////////////////////////////////////////////////
// Class:       anatree
// Plugin Type: analyzer (art v2_11_02)
// File:        anatree_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Additions from Leo Bellantoni, 2019-21
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

#include "SummaryDataProducts/RunData.h"
#include "SummaryDataProducts/POTSummary.h"

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

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"
#include "Geometry/BitFieldCoder.h"

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



namespace gar {

  typedef std::pair<float, std::string> P;

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

    virtual void endRun(art::Run const& run) override;
    virtual void beginJob() override;

    // Required functions.
    void analyze(art::Event const & e) override;

  private:
    void ClearVectors();

    //Working Horse
    cheat::BackTrackerCore* BackTrack;
    void FillGeneratorMonteCarloInfo(art::Event const & e);
    void FillRawInfo(art::Event const & e);
    void FillRecoInfo(art::Event const & e);
    void FillHighLevelRecoInfo(art::Event const & e);

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
    std::string fPOTtag;

    std::string fInstanceLabelCalo; ///< Instance name for ECAL
    std::string fInstanceLabelMuID; ///< Instance name for MuID
        
    std::string fGeantLabel; ///< module label for geant4 simulated hits

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

    //Geometry
    const geo::GeometryCore* fGeo; ///< pointer to the geometry
    std::string fECALEncoding;
    std::string fMuIDEncoding;
    gar::geo::BitFieldCoder *fFieldDecoder_ECAL;
    gar::geo::BitFieldCoder *fFieldDecoder_MuID;

    typedef int TrkId;
    std::unordered_map<TrkId, Int_t> TrackIdToIndex;


    // global event info
    Int_t   fEvent;      ///< number of the event being processed
    Int_t   fRun;        ///< number of the run being processed
    Int_t   fSubRun;     ///< number of the sub-run being processed
    Float_t fTPC_X;      ///< center of TPC stored as per-event & compressed by root
    Float_t fTPC_Y;
    Float_t fTPC_Z;
    Int_t   fTotalPOT;
    Int_t   fNSpills;

    // MCTruth data.
    // GENIE kinematics computed ignoring Fermi momentum and the off-shellness
    // of the bound nucleon.  Well, that's what the documentation says.  Might
    // not be true!  But to get the right kinematics here, t has to be computed.
    // Mode and InteractionType given in simb::MCNeutrino.h of the nusimdata
    // product.
    // Use Rtypes.h here, as these data get used by root

    std::vector<Int_t>              fNeutrinoType;
    std::vector<Int_t>              fCCNC;
    std::vector<Int_t>              fMode;
    std::vector<Int_t>              fInteractionType;
    std::vector<Float_t>            fQ2;
    std::vector<Float_t>            fW;
    std::vector<Float_t>            fX;
    std::vector<Float_t>            fY;
    std::vector<Float_t>            fTheta;
    std::vector<Float_t>            fT;
    std::vector<Float_t>            fMCVertexX;
    std::vector<Float_t>            fMCVertexY;
    std::vector<Float_t>            fMCVertexZ;
    std::vector<Float_t>            fMCnuPx;
    std::vector<Float_t>            fMCnuPy;
    std::vector<Float_t>            fMCnuPz;

    // GTruth data
    std::vector<Int_t>              fGint;
    std::vector<Int_t>              fTgtPDG;
    std::vector<Float_t>            fWeight;
    std::vector<Float_t>            fgT;

    //GENIE particle list from event record
    std::vector<Int_t>              fnGPart;
    std::vector<Int_t>              fGPartIntIdx;
    std::vector<Int_t>              fGPartIdx;
    std::vector<Int_t>              fGPartPdg;
    std::vector<Int_t>              fGPartStatus;
    std::vector<Int_t>              fGPartFirstMom;
    std::vector<Int_t>              fGPartLastMom;
    std::vector<Int_t>              fGPartFirstDaugh;
    std::vector<Int_t>              fGPartLastDaugh;
    std::vector<std::string>        fGPartName;
    std::vector<Float_t>            fGPartPx;
    std::vector<Float_t>            fGPartPy;
    std::vector<Float_t>            fGPartPz;
    std::vector<Float_t>            fGPartE;
    std::vector<Float_t>            fGPartMass;

    // MCParticle data
    std::vector<Int_t>              fMCTrkID;
    std::vector<Int_t>              fMCPDG;
    std::vector<Int_t>              fMCMotherIndex;
    std::vector<Int_t>              fMCMotherTrkID;
    std::vector<Int_t>              fMCPDGMother;
    std::vector<Float_t>            fMCPStartX;
    std::vector<Float_t>            fMCPStartY;
    std::vector<Float_t>            fMCPStartZ;
    std::vector<Float_t>            fMCPTime;
    std::vector<Float_t>            fMCPStartPX;
    std::vector<Float_t>            fMCPStartPY;
    std::vector<Float_t>            fMCPStartPZ;
    std::vector<Float_t>            fMCPEndX;
    std::vector<Float_t>            fMCPEndY;
    std::vector<Float_t>            fMCPEndZ;
    std::vector<Float_t>            fMCPEndPX;
    std::vector<Float_t>            fMCPEndPY;
    std::vector<Float_t>            fMCPEndPZ;
    std::vector<std::string>        fMCPProc;
    std::vector<std::string>        fMCPEndProc;
    std::vector<Int_t>              fMCPVertIndex;

    // track trajectory of MCP
    std::vector<Float_t>            fTrajMCPX;
    std::vector<Float_t>            fTrajMCPY;
    std::vector<Float_t>            fTrajMCPZ;
    std::vector<Float_t>            fTrajMCPT;
    std::vector<Float_t>            fTrajMCPE;
    std::vector<Int_t>              fTrajMCPIndex;
    std::vector<Int_t>              fTrajMCPTrackID;
    std::vector<Float_t>            fTrajMCPPX;
    std::vector<Float_t>            fTrajMCPPY;
    std::vector<Float_t>            fTrajMCPPZ;

    // sim calo hit data
    UInt_t                          fSimnHits;
    std::vector<Float_t>            fSimHitX;
    std::vector<Float_t>            fSimHitY;
    std::vector<Float_t>            fSimHitZ;
    std::vector<Float_t>            fSimHitTime;
    std::vector<Float_t>            fSimHitEnergy;
    std::vector<Int_t>              fSimHitTrackID;
    std::vector<ULong64_t>          fSimHitCellID;     // ULong64_t is size_t on 64 bit machines
    Float_t                         fSimEnergySum;

    //Muon system sim hits
    UInt_t                          fSimnHits_MuID;
    std::vector<Float_t>            fSimHitX_MuID;
    std::vector<Float_t>            fSimHitY_MuID;
    std::vector<Float_t>            fSimHitZ_MuID;
    std::vector<Float_t>            fSimHitTime_MuID;
    std::vector<Float_t>            fSimHitEnergy_MuID;
    std::vector<Int_t>              fSimHitTrackID_MuID;
    std::vector<ULong64_t>          fSimHitCellID_MuID;     // ULong64_t is size_t on 64 bit machines
    Float_t                         fSimEnergySum_MuID;

    // Hit data
    std::vector<Float_t>            fHitX;
    std::vector<Float_t>            fHitY;
    std::vector<Float_t>            fHitZ;
    std::vector<Float_t>            fHitSig;
    std::vector<Float_t>            fHitRMS;
    std::vector<UInt_t>             fHitChan;

    // TPCCluster data
    std::vector<Float_t>            fTPCClusterX;
    std::vector<Float_t>            fTPCClusterY;
    std::vector<Float_t>            fTPCClusterZ;
    std::vector<Float_t>            fTPCClusterSig;
    std::vector<Float_t>            fTPCClusterRMS;
    std::vector<ULong64_t>          fTPCClusterTrkIDNumber;
    std::vector<Float_t>            fTPCClusterCovXX;
    std::vector<Float_t>            fTPCClusterCovXY;
    std::vector<Float_t>            fTPCClusterCovXZ;
    std::vector<Float_t>            fTPCClusterCovYY;
    std::vector<Float_t>            fTPCClusterCovYZ;
    std::vector<Float_t>            fTPCClusterCovZZ;
    std::vector<Int_t>              fTPCClusterMCindex;      // Branch index (NOT the GEANT track ID) of MCParticle
    std::vector<Float_t>            fTPCClusterMCfrac;       // that best matches & fraction of ionization therefrom

    // track data
    std::vector<ULong64_t>          fTrackIDNumber;
    std::vector<Float_t>            fTrackStartX;
    std::vector<Float_t>            fTrackStartY;
    std::vector<Float_t>            fTrackStartZ;
    std::vector<Float_t>            fTrackStartPX;
    std::vector<Float_t>            fTrackStartPY;
    std::vector<Float_t>            fTrackStartPZ;
    std::vector<Int_t>              fTrackStartQ;

    std::vector<Float_t>            fTrackEndX;
    std::vector<Float_t>            fTrackEndY;
    std::vector<Float_t>            fTrackEndZ;
    std::vector<Float_t>            fTrackEndPX;
    std::vector<Float_t>            fTrackEndPY;
    std::vector<Float_t>            fTrackEndPZ;
    std::vector<Int_t>              fTrackEndQ;

    std::vector<Float_t>            fTrackLenF;         // from forward fit, from the Beg end to the End end
    std::vector<Float_t>            fTrackLenB;
    std::vector<Float_t>            fTrackChi2F;        // from forward fit, from the Beg end to the End end
    std::vector<Float_t>            fTrackChi2B;
    std::vector<Int_t>              fNTPCClustersOnTrack;
    std::vector<Float_t>            fTrackAvgIonF;      // from forward fit, from the Beg end to the End end
    std::vector<Float_t>            fTrackAvgIonB;

    std::vector<Int_t>              fTrackPIDF;
    std::vector<Float_t>            fTrackPIDProbF;
    std::vector<Int_t>              fTrackPIDB;
    std::vector<Float_t>            fTrackPIDProbB;
    std::vector<Int_t>              fTrackMCindex;      // Branch index (NOT the GEANT track ID) of MCParticle
    std::vector<Float_t>            fTrackMCfrac;       // that best matches & fraction of ionization therefrom

    //TrackTrajectory
    std::vector<Float_t>            fTrackTrajectoryFWDX; //forward
    std::vector<Float_t>            fTrackTrajectoryFWDY; //forward
    std::vector<Float_t>            fTrackTrajectoryFWDZ; //forward
    std::vector<Int_t>              fTrackTrajectoryFWDID;

    std::vector<Float_t>            fTrackTrajectoryBWDX; //backward
    std::vector<Float_t>            fTrackTrajectoryBWDY; //backward
    std::vector<Float_t>            fTrackTrajectoryBWDZ; //backward
    std::vector<Int_t>              fTrackTrajectoryBWDID;

    // vertex branches
    std::vector<ULong64_t>          fVertexIDNumber;
    std::vector<Float_t>            fVertexX;
    std::vector<Float_t>            fVertexY;
    std::vector<Float_t>            fVertexZ;
    std::vector<ULong64_t>          fVertexT;
    std::vector<Int_t>              fVertexN;
    std::vector<Int_t>              fVertexQ;

    std::vector<ULong64_t>          fVTAssn_VertIDNumber;     // Being the vertex which this Assn belongs to
    std::vector<ULong64_t>          fVTAssn_TrackIDNumber;
    std::vector<gar::rec::TrackEnd> fVTAssn_TrackEnd;

    // Vee branches
    std::vector<ULong64_t>          fVeeIDNumber;
    std::vector<Float_t>            fVeeX;
    std::vector<Float_t>            fVeeY;
    std::vector<Float_t>            fVeeZ;
    std::vector<ULong64_t>          fVeeT;
    std::vector<Float_t>            fVeePXKpipi;
    std::vector<Float_t>            fVeePYKpipi;
    std::vector<Float_t>            fVeePZKpipi;
    std::vector<Float_t>            fVeeEKpipi;
    std::vector<Float_t>            fVeeMKpipi;
    std::vector<Float_t>            fVeePXLppi;
    std::vector<Float_t>            fVeePYLppi;
    std::vector<Float_t>            fVeePZLppi;
    std::vector<Float_t>            fVeeELppi;
    std::vector<Float_t>            fVeeMLppi;
    std::vector<Float_t>            fVeePXLpip;
    std::vector<Float_t>            fVeePYLpip;
    std::vector<Float_t>            fVeePZLpip;
    std::vector<Float_t>            fVeeELpip;
    std::vector<Float_t>            fVeeMLpip;
    std::vector<ULong64_t>          fVeeTAssn_VeeIDNumber;    // Being the Vee which this Assn belongs to
    std::vector<ULong64_t>          fVeeTAssn_TrackIDNumber;
    std::vector<gar::rec::TrackEnd> fVeeTAssn_TrackEnd;

    // raw calo digits data
    UInt_t                          fDiginHits;
    std::vector<Float_t>            fDigiHitX;
    std::vector<Float_t>            fDigiHitY;
    std::vector<Float_t>            fDigiHitZ;
    std::vector<Float_t>            fDigiHitTime;
    std::vector<UInt_t>             fDigiHitADC;              // UInt_t is unsigned 32 bit integer
    std::vector<ULong64_t>          fDigiHitCellID;

    //Muon system raw hits
    UInt_t                          fDiginHits_MuID;
    std::vector<Float_t>            fDigiHitX_MuID;
    std::vector<Float_t>            fDigiHitY_MuID;
    std::vector<Float_t>            fDigiHitZ_MuID;
    std::vector<Float_t>            fDigiHitTime_MuID;
    std::vector<UInt_t>             fDigiHitADC_MuID;         // UInt_t is unsigned 32 bit integer
    std::vector<ULong64_t>          fDigiHitCellID_MuID;

    // reco calo hit data
    UInt_t                          fReconHits;
    std::vector<ULong64_t>          fReconHitIDNumber;
    std::vector<Float_t>            fRecoHitX;
    std::vector<Float_t>            fRecoHitY;
    std::vector<Float_t>            fRecoHitZ;
    std::vector<Float_t>            fRecoHitTime;
    std::vector<Float_t>            fRecoHitEnergy;                     // Hit energies have the sampling fraction
    std::vector<ULong64_t>          fRecoHitCellID;                     // correction factors in Reco/SiPMHitFinder.fcl
    std::vector<Int_t>              fRecoHitLayer;                      // These are propagated to clusters and the same
    Float_t                         fRecoEnergySum;                     // is done for the NuID system.

    //Muon system reco hits
    UInt_t                          fReconHits_MuID;
    std::vector<ULong64_t>          fReconHitIDNumber_MuID;
    std::vector<Float_t>            fRecoHitX_MuID;
    std::vector<Float_t>            fRecoHitY_MuID;
    std::vector<Float_t>            fRecoHitZ_MuID;
    std::vector<Float_t>            fRecoHitTime_MuID;
    std::vector<Float_t>            fRecoHitEnergy_MuID;
    std::vector<ULong64_t>          fRecoHitCellID_MuID;
    std::vector<Int_t>              fRecoHitLayer_MuID;
    Float_t                         fRecoEnergySum_MuID;

    // calo cluster data
    UInt_t                          fnCluster;
    std::vector<ULong64_t>          fClusterIDNumber;
    std::vector<UInt_t>             fClusterNhits;
    std::vector<Float_t>            fClusterEnergy;
    std::vector<Float_t>            fClusterTime;
    std::vector<Float_t>            fClusterTimeDiffFirstLast;
    std::vector<Float_t>            fClusterX;
    std::vector<Float_t>            fClusterY;
    std::vector<Float_t>            fClusterZ;
    std::vector<Float_t>            fClusterTheta;
    std::vector<Float_t>            fClusterPhi;
    std::vector<Float_t>            fClusterPID;
    // std::vector<Float_t> fClusterShape;
    std::vector<Float_t>            fClusterMainAxisX;
    std::vector<Float_t>            fClusterMainAxisY;
    std::vector<Float_t>            fClusterMainAxisZ;
    std::vector<Int_t>              fClusterMCindex;          // Branch index (NOT the GEANT track ID) of MCParticle
    std::vector<Float_t>            fClusterMCfrac;           // that best matches & fraction of ionization therefrom

    // calo cluster data
    UInt_t                          fnCluster_MuID;
    std::vector<ULong64_t>          fClusterIDNumber_MuID;
    std::vector<UInt_t>             fClusterNhits_MuID;
    std::vector<Float_t>            fClusterEnergy_MuID;
    std::vector<Float_t>            fClusterTime_MuID;
    std::vector<Float_t>            fClusterTimeDiffFirstLast_MuID;
    std::vector<Float_t>            fClusterX_MuID;
    std::vector<Float_t>            fClusterY_MuID;
    std::vector<Float_t>            fClusterZ_MuID;
    std::vector<Float_t>            fClusterTheta_MuID;
    std::vector<Float_t>            fClusterPhi_MuID;
    std::vector<Float_t>            fClusterPID_MuID;
    // std::vector<Float_t> fClusterShape;
    std::vector<Float_t>            fClusterMainAxisX_MuID;
    std::vector<Float_t>            fClusterMainAxisY_MuID;
    std::vector<Float_t>            fClusterMainAxisZ_MuID;
    std::vector<Int_t>              fClusterMCindex_MuID;          // Branch index (NOT the GEANT track ID) of MCParticle
    std::vector<Float_t>            fClusterMCfrac_MuID;           // that best matches & fraction of ionization therefrom

    // ECAL cluster to ECAL hits association info
    std::vector<std::vector<ULong64_t>>          fClusterAssn_RecoHitIDNumber;  

    // MuID cluster to MuID hits association info
    std::vector<std::vector<ULong64_t>>          fClusterMuIDAssn_MuIDHitIDNumber;  

    // ECAL cluster to track association info
    std::vector<ULong64_t>          fCALAssn_ClusIDNumber;   // Being the cluster which this Assn belongs to
    std::vector<ULong64_t>          fCALAssn_TrackIDNumber;  // The rec::TrackEnd (see Track.h) that extrapolated to cluster
    std::vector<gar::rec::TrackEnd> fCALAssn_TrackEnd;

    //map of PID TH2 per momentum value
    CLHEP::HepRandomEngine &fEngine;  ///< random engine
    std::unordered_map<int, TH2F*> m_pidinterp;
  };
}



//==============================================================================
//==============================================================================
//==============================================================================
// constructor
gar::anatree::anatree(fhicl::ParameterSet const & p)
  : EDAnalyzer(p),
    fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, p, "Seed")) {
  fGeo     = gar::providerFrom<geo::GeometryGAr>();

  fECALEncoding = fGeo->GetECALCellIDEncoding();
  fFieldDecoder_ECAL = new gar::geo::BitFieldCoder( fECALEncoding );

  if(fGeo->HasMuonDetector()) {
    fMuIDEncoding = fGeo->GetMuIDCellIDEncoding();
    fFieldDecoder_MuID = new gar::geo::BitFieldCoder( fMuIDEncoding );
  }

  bool usegenlabels =
    p.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
  if (!usegenlabels) fGeneratorLabels.clear();

  bool usegeniegenlabels  =
    p.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
  if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

  fPOTtag            = p.get<std::string>("SummaryDataLabel","generator");

  //Sim Hits
  fGeantLabel        = p.get<std::string>("GEANTLabel","geant");
  fInstanceLabelCalo = p.get<std::string>("InstanceLabelCalo","ECAL");
  fInstanceLabelMuID = p.get<std::string>("InstanceLabelMuID","MuID");

  fHitLabel          = p.get<std::string>("HitLabel","hit");
  fTPCClusterLabel   = p.get<std::string>("TPCClusterLabel","tpccluster");
  fTrackLabel        = p.get<std::string>("TrackLabel","track");
  fTrackTrajectoryLabel  = p.get<std::string>("TrackTrajectoryLabel","track");
  fVertexLabel       = p.get<std::string>("VertexLabel","vertex");
  fVeeLabel          = p.get<std::string>("VeeLabel","veefinder1");

  //Calorimetric related ECAL/MuID
  fRawCaloHitLabel   = p.get<std::string>("RawCaloHitLabel","daqsipm");
  fRawMuIDHitLabel   = p.get<std::string>("RawMuIDHitLabel","daqsipmmuid");
  fCaloHitLabel      = p.get<std::string>("CaloHitLabel","sscalohit");
  fMuIDHitLabel      = p.get<std::string>("MuIDHitLabel","sscalohitmuid");

  fClusterLabel      = p.get<std::string>("ClusterLabel","calocluster");
  fClusterMuIDLabel  = p.get<std::string>("MuIDClusterLabel","caloclustermuid");
  fPFLabel           = p.get<std::string>("PFLabel","pandora");
  fECALAssnLabel     = p.get<std::string>("ECALAssnLabel","trkecalassn");

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
void gar::anatree::beginJob() {

  fTPC_X = ItsInTulsa[0] = fGeo->TPCXCent();
  fTPC_Y = ItsInTulsa[1] = fGeo->TPCYCent();
  fTPC_Z = ItsInTulsa[2] = fGeo->TPCZCent();

  xTPC = fGeo->TPCLength() / 2.;
  rTPC = fGeo->TPCRadius();

  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

  fTree->Branch("Run",           &fRun,         "Run/I");
  fTree->Branch("SubRun",        &fSubRun,      "SubRun/I");
  fTree->Branch("Event",         &fEvent,       "Event/I");
  fTree->Branch("TPC_X",         &fTPC_X,       "TPC_X/F");
  fTree->Branch("TPC_Y",         &fTPC_Y,       "TPC_Y/F");
  fTree->Branch("TPC_Z",         &fTPC_Z,       "TPC_Z/F");
  fTree->Branch("POT",           &fTotalPOT,    "POT/I");
  fTree->Branch("NSpills",       &fNSpills,     "NSpills/I");

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

    //GENIE particle list from the event record
    fTree->Branch("nGPart",         &fnGPart);
    fTree->Branch("GPartIntIdx",    &fGPartIntIdx);
    fTree->Branch("GPartIdx",       &fGPartIdx);
    fTree->Branch("GPartName",      &fGPartName);
    fTree->Branch("GPartPdg",       &fGPartPdg);
    fTree->Branch("GPartStatus",    &fGPartStatus);
    fTree->Branch("GPartFirstMom",  &fGPartFirstMom);
    fTree->Branch("GPartLastMom",   &fGPartLastMom);
    fTree->Branch("GPartFirstDaugh",&fGPartFirstDaugh);
    fTree->Branch("GPartLastDaugh", &fGPartLastDaugh);
    fTree->Branch("GPartPx",        &fGPartPx);
    fTree->Branch("GPartPy",        &fGPartPy);
    fTree->Branch("GPartPz",        &fGPartPz);
    fTree->Branch("GPartE",         &fGPartE);
    fTree->Branch("GPartMass",      &fGPartMass);

    fTree->Branch("MCTrkID",     &fMCTrkID);
    fTree->Branch("PDG",         &fMCPDG);
    fTree->Branch("MotherIndex", &fMCMotherIndex);    // Index into these vector branches
    fTree->Branch("MotherTrkID", &fMCMotherTrkID);    // trackid of the mother
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
    fTree->Branch("MCPVertIndex",&fMCPVertIndex);
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
      if (fWriteMCPTrajMomenta) {
        fTree->Branch("TrajMCPPX",          &fTrajMCPPX);
        fTree->Branch("TrajMCPPY",          &fTrajMCPPY);
        fTree->Branch("TrajMCPPZ",          &fTrajMCPPZ);
      }
      fTree->Branch("TrajMCPE",          &fTrajMCPE);
      fTree->Branch("TrajMCPIndex",      &fTrajMCPIndex);
      fTree->Branch("TrajMCPTrackID",    &fTrajMCPTrackID);
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

      if(fGeo->HasMuonDetector() && fWriteMuID) {
        fTree->Branch("SimnHits_MuID",     &fSimnHits_MuID);
        fTree->Branch("SimHitX_MuID",      &fSimHitX_MuID);
        fTree->Branch("SimHitY_MuID",      &fSimHitY_MuID);
        fTree->Branch("SimHitZ_MuID",      &fSimHitZ_MuID);
        fTree->Branch("SimHitTime_MuID",   &fSimHitTime_MuID);
        fTree->Branch("SimHitEnergy_MuID", &fSimHitEnergy_MuID);
        fTree->Branch("SimHitTrkID_MuID",  &fSimHitTrackID_MuID);
        fTree->Branch("SimHitCellID_MuID", &fSimHitCellID_MuID);
        fTree->Branch("SimEnergySum_MuID", &fSimEnergySum_MuID);
      }
    }
  }

  if (fWriteHits) {
    fTree->Branch("HitX",        &fHitX);
    fTree->Branch("HitY",        &fHitY);
    fTree->Branch("HitZ",        &fHitZ);
    fTree->Branch("HitSig",      &fHitSig);
    fTree->Branch("HitRMS",      &fHitRMS);
    fTree->Branch("HitChan",     &fHitChan);
  }

  if (fWriteTPCClusters) {
    fTree->Branch("TPCClusterX",           &fTPCClusterX);
    fTree->Branch("TPCClusterY",           &fTPCClusterY);
    fTree->Branch("TPCClusterZ",           &fTPCClusterZ);
    fTree->Branch("TPCClusterSig",         &fTPCClusterSig);
    fTree->Branch("TPCClusterRMS",         &fTPCClusterRMS);
    fTree->Branch("TPCClusterTrkIDNumber", &fTPCClusterTrkIDNumber);
    fTree->Branch("TPCClusterCovXX",       &fTPCClusterCovXX);
    fTree->Branch("TPCClusterCovXY",       &fTPCClusterCovXY);
    fTree->Branch("TPCClusterCovXZ",       &fTPCClusterCovXZ);
    fTree->Branch("TPCClusterCovYY",       &fTPCClusterCovYY);
    fTree->Branch("TPCClusterCovYZ",       &fTPCClusterCovYZ);
    fTree->Branch("TPCClusterCovZZ",       &fTPCClusterCovZZ);
    fTree->Branch("TPCClusterMCindex",     &fTPCClusterMCindex);
    fTree->Branch("TPCClusterMCfrac",      &fTPCClusterMCfrac);
  }

  // All position, momentum, etc
  if (fWriteTracks) {
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
    fTree->Branch("TrackChi2F",         &fTrackChi2F);
    fTree->Branch("TrackChi2B",         &fTrackChi2B);
    fTree->Branch("NTPCClustersOnTrack",&fNTPCClustersOnTrack);
    fTree->Branch("TrackAvgIonF",       &fTrackAvgIonF);
    fTree->Branch("TrackAvgIonB",       &fTrackAvgIonB);

    fTree->Branch("TrackPIDF",          &fTrackPIDF);
    fTree->Branch("TrackPIDProbF",      &fTrackPIDProbF);
    fTree->Branch("TrackPIDB",          &fTrackPIDB);
    fTree->Branch("TrackPIDProbB",      &fTrackPIDProbB);
    fTree->Branch("TrackMCindex",       &fTrackMCindex);
    fTree->Branch("TrackMCfrac",        &fTrackMCfrac);

    //Track Trajectories
    if(fWriteTrackTrajectories) {
      fTree->Branch("TrackTrajectoryFWDX",    &fTrackTrajectoryFWDX);
      fTree->Branch("TrackTrajectoryFWDY",    &fTrackTrajectoryFWDY);
      fTree->Branch("TrackTrajectoryFWDZ",    &fTrackTrajectoryFWDZ);
      fTree->Branch("TrackTrajectoryFWDID",   &fTrackTrajectoryFWDID);

      fTree->Branch("TrackTrajectoryBWDX",    &fTrackTrajectoryBWDX);
      fTree->Branch("TrackTrajectoryBWDY",    &fTrackTrajectoryBWDY);
      fTree->Branch("TrackTrajectoryBWDZ",    &fTrackTrajectoryBWDZ);
      fTree->Branch("TrackTrajectoryBWDID",     &fTrackTrajectoryBWDID);
    }
  }

  // Reco'd verts & their track-ends
  if (fWriteVertices) {
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

  // Reco'd vees & their track-ends
  if (fWriteVees) {
    if (!fWriteTracks) {
      throw cet::exception("anatree")
        << " fWriteVees, but !fWriteTracks."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    } else {
      fTree->Branch("VeeIDNumber",     &fVeeIDNumber);
      fTree->Branch("VeeX",            &fVeeX);
      fTree->Branch("VeeY",            &fVeeY);
      fTree->Branch("VeeZ",            &fVeeZ);
      fTree->Branch("VeeT",            &fVeeT);
      fTree->Branch("VeePXKpipi",      &fVeePXKpipi);
      fTree->Branch("VeePYKpipi",      &fVeePYKpipi);
      fTree->Branch("VeePZKpipi",      &fVeePZKpipi);
      fTree->Branch("VeeEKpipi",       &fVeeEKpipi);
      fTree->Branch("VeeMKpipi",       &fVeeMKpipi);
      fTree->Branch("VeePXLppi",       &fVeePXLppi);
      fTree->Branch("VeePYLppi",       &fVeePYLppi);
      fTree->Branch("VeePZLppi",       &fVeePZLppi);
      fTree->Branch("VeeELppi",        &fVeeELppi);
      fTree->Branch("VeeMLppi",        &fVeeMLppi);
      fTree->Branch("VeePXLpip",       &fVeePXLpip);
      fTree->Branch("VeePYLpip",       &fVeePYLpip);
      fTree->Branch("VeePZLpip",       &fVeePZLpip);
      fTree->Branch("VeeELpip",        &fVeeELpip);
      fTree->Branch("VeeMLpip",        &fVeeMLpip);
      fTree->Branch("VeeT_VertIDNumber", &fVeeTAssn_VeeIDNumber);
      fTree->Branch("VeeT_TrackIDNumber",&fVeeTAssn_TrackIDNumber);
      fTree->Branch("VeeT_TrackEnd",     &fVeeTAssn_TrackEnd);
    }
  }

  // Write calorimetry digits
  if (fWriteCaloDigits) {
    fTree->Branch("DiginHits",        &fDiginHits);
    fTree->Branch("DigiHitX",         &fDigiHitX);
    fTree->Branch("DigiHitY",         &fDigiHitY);
    fTree->Branch("DigiHitZ",         &fDigiHitZ);
    fTree->Branch("DigiHitTime",      &fDigiHitTime);
    fTree->Branch("DigiHitADC",       &fDigiHitADC);
    fTree->Branch("DigiHitCellID",    &fDigiHitCellID);

    if(fGeo->HasMuonDetector()) {
      fTree->Branch("DiginHits_MuID",        &fDiginHits_MuID);
      fTree->Branch("DigiHitX_MuID",         &fDigiHitX_MuID);
      fTree->Branch("DigiHitY_MuID",         &fDigiHitY_MuID);
      fTree->Branch("DigiHitZ_MuID",         &fDigiHitZ_MuID);
      fTree->Branch("DigiHitTime_MuID",      &fDigiHitTime_MuID);
      fTree->Branch("DigiHitADC_MuID",       &fDigiHitADC_MuID);
      fTree->Branch("DigiHitCellID_MuID",    &fDigiHitCellID_MuID);
    }
  }

  // Write calorimetry hits
  if (fWriteCaloHits) {
    fTree->Branch("ReconHits",        &fReconHits);
    fTree->Branch("ReconHitIDNumber", &fReconHitIDNumber);
    fTree->Branch("RecoHitX",         &fRecoHitX);
    fTree->Branch("RecoHitY",         &fRecoHitY);
    fTree->Branch("RecoHitZ",         &fRecoHitZ);
    fTree->Branch("RecoHitTime",      &fRecoHitTime);
    fTree->Branch("RecoHitEnergy",    &fRecoHitEnergy);
    fTree->Branch("RecoHitCellID",    &fRecoHitCellID);
    fTree->Branch("RecoHitLayer",     &fRecoHitLayer);
    fTree->Branch("RecoEnergySum",    &fRecoEnergySum);

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fTree->Branch("ReconHits_MuID",         &fReconHits_MuID);
      fTree->Branch("ReconHitIDNumber_MuID",    &fReconHitIDNumber_MuID);
      fTree->Branch("RecoHitX_MuID",          &fRecoHitX_MuID);
      fTree->Branch("RecoHitY_MuID",          &fRecoHitY_MuID);
      fTree->Branch("RecoHitZ_MuID",          &fRecoHitZ_MuID);
      fTree->Branch("RecoHitTime_MuID",       &fRecoHitTime_MuID);
      fTree->Branch("RecoHitEnergy_MuID",     &fRecoHitEnergy_MuID);
      fTree->Branch("RecoHitCellID_MuID",     &fRecoHitCellID_MuID);
      fTree->Branch("RecoHitLayer_MuID",      &fRecoHitLayer_MuID);
      fTree->Branch("RecoEnergySum_MuID",     &fRecoEnergySum_MuID);
    }
  }

  // Write calorimetry clusters
  if (fWriteCaloClusters) {
    fTree->Branch("nCluster",                   &fnCluster);
    fTree->Branch("ClusterIDNumber",            &fClusterIDNumber);
    fTree->Branch("ClusterNhits",               &fClusterNhits);
    fTree->Branch("ClusterEnergy",              &fClusterEnergy);
    fTree->Branch("ClusterTime",                &fClusterTime);
    fTree->Branch("ClusterTimeDiffFirstLast",   &fClusterTimeDiffFirstLast);
    fTree->Branch("ClusterX",                   &fClusterX);
    fTree->Branch("ClusterY",                   &fClusterY);
    fTree->Branch("ClusterZ",                   &fClusterZ);
    fTree->Branch("ClusterTheta",               &fClusterTheta);
    fTree->Branch("ClusterPhi",                 &fClusterPhi);
    fTree->Branch("ClusterPID",                 &fClusterPID);
    // fTree->Branch("ClusterShape",            &fClusterShape);
    fTree->Branch("ClusterMainAxisX",           &fClusterMainAxisX);
    fTree->Branch("ClusterMainAxisY",           &fClusterMainAxisY);
    fTree->Branch("ClusterMainAxisZ",           &fClusterMainAxisZ);
    fTree->Branch("ClusterMCindex",             &fClusterMCindex);
    fTree->Branch("ClusterMCfrac",              &fClusterMCfrac);

    fTree->Branch("ClusterAssn_RecoHitIDNumber", &fClusterAssn_RecoHitIDNumber);

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fTree->Branch("nCluster_MuID",                   &fnCluster_MuID);
      fTree->Branch("ClusterIDNumber_MuID",            &fClusterIDNumber_MuID);
      fTree->Branch("ClusterNhits_MuID",               &fClusterNhits_MuID);
      fTree->Branch("ClusterEnergy_MuID",              &fClusterEnergy_MuID);
      fTree->Branch("ClusterTime_MuID",                &fClusterTime_MuID);
      fTree->Branch("ClusterTimeDiffFirstLast_MuID",   &fClusterTimeDiffFirstLast_MuID);
      fTree->Branch("ClusterX_MuID",                   &fClusterX_MuID);
      fTree->Branch("ClusterY_MuID",                   &fClusterY_MuID);
      fTree->Branch("ClusterZ_MuID",                   &fClusterZ_MuID);
      fTree->Branch("ClusterTheta_MuID",               &fClusterTheta_MuID);
      fTree->Branch("ClusterPhi_MuID",                 &fClusterPhi_MuID);
      fTree->Branch("ClusterPID_MuID",                 &fClusterPID_MuID);
      // fTree->Branch("ClusterShape",            &fClusterShape);
      fTree->Branch("ClusterMainAxisX_MuID",           &fClusterMainAxisX_MuID);
      fTree->Branch("ClusterMainAxisY_MuID",           &fClusterMainAxisY_MuID);
      fTree->Branch("ClusterMainAxisZ_MuID",           &fClusterMainAxisZ_MuID);
      fTree->Branch("ClusterMCindex_MuID",             &fClusterMCindex_MuID);
      fTree->Branch("ClusterMCfrac_MuID",              &fClusterMCfrac_MuID);

      fTree->Branch("ClusterMuIDAssn_MuIDHitIDNumber", &fClusterMuIDAssn_MuIDHitIDNumber);
    }
  }

  if (fWriteMatchedTracks) {
    if (!fWriteTracks || !fWriteCaloClusters) {
      throw cet::exception("anatree")
        << " fWriteMatchedTracks, but (!fWriteTracks || !fWriteCaloClusters)."
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    } else {
      fTree->Branch("ECALAssn_ClusIDNumber",  &fCALAssn_ClusIDNumber);
      fTree->Branch("ECALAssn_TrackIDNumber", &fCALAssn_TrackIDNumber);
      fTree->Branch("ECALAssn_TrackEnd",      &fCALAssn_TrackEnd);
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
}  // End of :anatree::beginJob

//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::endRun(art::Run const& run) {
  auto const& ID = run.id();

  auto summaryHandle = run.getHandle<sumdata::POTSummary>(fPOTtag);
  if (!summaryHandle) {
    MF_LOG_DEBUG("anatree") << " No sumdata::POTSummary branch for run " << ID 
                            <<" Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    return;
  }

  sumdata::POTSummary const& RunPOT = *summaryHandle;
  fTotalPOT = RunPOT.TotalPOT();
  fNSpills = RunPOT.TotalSpills();

  MF_LOG_INFO("anatree") << "POT for this file is " << fTotalPOT
                         << " The number of spills is " << fNSpills;
}

//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::analyze(art::Event const & e) {

  ClearVectors();

  fRun    = e.run();
  fSubRun = e.subRun();
  fEvent  = e.id().event();

  // Need a non-constant backtracker instance, for now, in analyze not beginJob
  cheat::BackTrackerCore const* const_bt = gar::providerFrom<cheat::BackTracker>();
  BackTrack = const_cast<cheat::BackTrackerCore*>(const_bt);

  //Fill generator and MC Information
  if (fWriteMCinfo) FillGeneratorMonteCarloInfo(e);

  //Fill raw calorimeter information
  if (fWriteCaloDigits) FillRawInfo(e);

  //Fill reco information
  FillRecoInfo(e);

  //Fill high-leve reco information
  FillHighLevelRecoInfo(e);

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

    fnGPart.clear();
    fGPartIntIdx.clear();
    fGPartIdx.clear();
    fGPartPdg.clear();
    fGPartStatus.clear();
    fGPartName.clear();
    fGPartFirstMom.clear();
    fGPartLastMom.clear();
    fGPartFirstDaugh.clear();
    fGPartLastDaugh.clear();
    fGPartPx.clear();
    fGPartPy.clear();
    fGPartPz.clear();
    fGPartE.clear();
    fGPartMass.clear();

    fMCTrkID.clear();
    fMCPDG.clear();
    fMCMotherIndex.clear();
    fMCMotherTrkID.clear();
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
    fMCPVertIndex.clear();
  }

  if (fWriteMCPTrajectory) {
    fTrajMCPX.clear();
    fTrajMCPY.clear();
    fTrajMCPZ.clear();
    fTrajMCPT.clear();
    if (fWriteMCPTrajMomenta) {
      fTrajMCPPX.clear();
      fTrajMCPPY.clear();
      fTrajMCPPZ.clear();
    }
    fTrajMCPE.clear();
    fTrajMCPIndex.clear();
    fTrajMCPTrackID.clear();
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

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fSimnHits_MuID = 0;
      fSimHitX_MuID.clear();
      fSimHitY_MuID.clear();
      fSimHitZ_MuID.clear();
      fSimHitTime_MuID.clear();
      fSimHitEnergy_MuID.clear();
      fSimHitTrackID_MuID.clear();
      fSimHitCellID_MuID.clear();
      fSimEnergySum_MuID = 0.;
    }
  }

  if (fWriteHits) {
    fHitX.clear();
    fHitY.clear();
    fHitZ.clear();
    fHitSig.clear();
    fHitRMS.clear();
    fHitChan.clear();
  }

  if (fWriteTPCClusters) {
    fTPCClusterX.clear();
    fTPCClusterY.clear();
    fTPCClusterZ.clear();
    fTPCClusterSig.clear();
    fTPCClusterRMS.clear();
    fTPCClusterTrkIDNumber.clear();
    fTPCClusterCovXX.clear();
    fTPCClusterCovXY.clear();
    fTPCClusterCovXZ.clear();
    fTPCClusterCovYY.clear();
    fTPCClusterCovYZ.clear();
    fTPCClusterCovZZ.clear();
    fTPCClusterMCindex.clear();
    fTPCClusterMCfrac.clear();
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
    fTrackChi2F.clear();
    fTrackChi2B.clear();
    fNTPCClustersOnTrack.clear();
    fTrackAvgIonF.clear();
    fTrackAvgIonB.clear();
    fTrackPIDF.clear();
    fTrackPIDProbF.clear();
    fTrackPIDB.clear();
    fTrackPIDProbB.clear();
    fTrackMCindex.clear();
    fTrackMCfrac.clear();

    if(fWriteTrackTrajectories) {
      fTrackTrajectoryFWDX.clear();
      fTrackTrajectoryFWDY.clear();
      fTrackTrajectoryFWDZ.clear();
      fTrackTrajectoryFWDID.clear();

      fTrackTrajectoryBWDX.clear();
      fTrackTrajectoryBWDY.clear();
      fTrackTrajectoryBWDZ.clear();
      fTrackTrajectoryBWDID.clear();
    }
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

  if (fWriteVees) {
    fVeeIDNumber.clear();
    fVeeX.clear();
    fVeeY.clear();
    fVeeZ.clear();
    fVeeT.clear();
    fVeePXKpipi.clear();
    fVeePYKpipi.clear();
    fVeePZKpipi.clear();
    fVeeEKpipi.clear();
    fVeeMKpipi.clear();
    fVeePXLppi.clear();
    fVeePYLppi.clear();
    fVeePZLppi.clear();
    fVeeELppi.clear();
    fVeeMLppi.clear();
    fVeePXLpip.clear();
    fVeePYLpip.clear();
    fVeePZLpip.clear();
    fVeeELpip.clear();
    fVeeMLpip.clear();
    fVeeTAssn_VeeIDNumber.clear();
    fVeeTAssn_TrackIDNumber.clear();
    fVeeTAssn_TrackEnd.clear();
  }

  if (fWriteCaloDigits) {
    fDiginHits = 0;
    fDigiHitX.clear();
    fDigiHitY.clear();
    fDigiHitZ.clear();
    fDigiHitTime.clear();
    fDigiHitADC.clear();
    fDigiHitCellID.clear();

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fDiginHits_MuID = 0;
      fDigiHitX_MuID.clear();
      fDigiHitY_MuID.clear();
      fDigiHitZ_MuID.clear();
      fDigiHitTime_MuID.clear();
      fDigiHitADC_MuID.clear();
      fDigiHitCellID_MuID.clear();
    }
  }

  if (fWriteCaloHits) {
    fReconHits = 0;
    fReconHitIDNumber.clear();
    fRecoHitX.clear();
    fRecoHitY.clear();
    fRecoHitZ.clear();
    fRecoHitTime.clear();
    fRecoHitEnergy.clear();
    fRecoHitCellID.clear();
    fRecoHitLayer.clear();
    fRecoEnergySum = 0.;

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fReconHits_MuID = 0;
      fReconHitIDNumber_MuID.clear();
      fRecoHitX_MuID.clear();
      fRecoHitY_MuID.clear();
      fRecoHitZ_MuID.clear();
      fRecoHitTime_MuID.clear();
      fRecoHitEnergy_MuID.clear();
      fRecoHitCellID_MuID.clear();
      fRecoHitLayer_MuID.clear();
      fRecoEnergySum_MuID = 0.;
    }
  }

  if (fWriteCaloClusters) {
    fnCluster = 0;
    fClusterIDNumber.clear();
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
    // fClusterShape.clear();
    fClusterMainAxisX.clear();
    fClusterMainAxisY.clear();
    fClusterMainAxisZ.clear();
    fClusterMCindex.clear();
    fClusterMCfrac.clear();

    fClusterAssn_RecoHitIDNumber.clear();

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      fnCluster_MuID = 0;
      fClusterIDNumber_MuID.clear();
      fClusterNhits_MuID.clear();
      fClusterEnergy_MuID.clear();
      fClusterTime_MuID.clear();
      fClusterTimeDiffFirstLast_MuID.clear();
      fClusterX_MuID.clear();
      fClusterY_MuID.clear();
      fClusterZ_MuID.clear();
      fClusterTheta_MuID.clear();
      fClusterPhi_MuID.clear();
      fClusterPID_MuID.clear();
      // fClusterShape.clear();
      fClusterMainAxisX_MuID.clear();
      fClusterMainAxisY_MuID.clear();
      fClusterMainAxisZ_MuID.clear();
      fClusterMCindex_MuID.clear();
      fClusterMCfrac_MuID.clear();

      fClusterMuIDAssn_MuIDHitIDNumber.clear();
    }
  }

  if (fWriteMatchedTracks) {
    fCALAssn_ClusIDNumber.clear();
    fCALAssn_TrackIDNumber.clear();
    fCALAssn_TrackEnd.clear();
  }

  return;
} // end :anatree::ClearVectors



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillGeneratorMonteCarloInfo(art::Event const & e) {

  // =============  Get art handles ==========================================
  // Get handles for MCinfo, also good for MCPTrajectory
  std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;
  std::vector< art::Handle< std::vector<simb::GTruth> > > gthandlelist;

  if (fGeneratorLabels.size()<1) {
    mcthandlelist = e.getMany<std::vector<simb::MCTruth> >(); // get them all (even if there are none)
  } else {
    mcthandlelist.resize(fGeneratorLabels.size());
    for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
      // complain if we wanted a specific one but didn't find it
      mcthandlelist.at(i) = e.getHandle<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
      if (!mcthandlelist.at(i)) {
        throw cet::exception("anatree") << " No simb::MCTruth branch."
                                        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }
  }

  if (fGENIEGeneratorLabels.size()<1) {
    gthandlelist = e.getMany< std::vector<simb::GTruth> >();  // get them all (even if there are none)
  } else {
    gthandlelist.resize(fGENIEGeneratorLabels.size());
    for (size_t i=0; i< fGENIEGeneratorLabels.size(); ++i) {
      // complain if we wanted a specific one but didn't find it
      gthandlelist.at(i) = e.getHandle<std::vector<simb::GTruth> >(fGENIEGeneratorLabels.at(i));
      if (!gthandlelist.at(i)) {
        throw cet::exception("anatree") << " No simb::GTruth branch."
                                        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }
  }

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

  // Save the particle list from the GENIE event record
  auto gparthandlelist = e.getMany<std::vector<sdp::GenieParticle>>();

  for (size_t igphl = 0; igphl < gparthandlelist.size(); ++igphl) {
    unsigned int nGPart = 0;
    for ( auto const& gpart : (*gparthandlelist.at(igphl)) ) {
      fGPartIntIdx.push_back(gpart.InteractionIndex());
      fGPartIdx.push_back(gpart.Index());
      fGPartPdg.push_back(gpart.Pdg());
      fGPartStatus.push_back(gpart.Status());
      fGPartName.push_back(gpart.Name());
      fGPartFirstMom.push_back(gpart.FirstMother());
      fGPartLastMom.push_back(gpart.LastMother());
      fGPartFirstDaugh.push_back(gpart.FirstDaughter());
      fGPartLastDaugh.push_back(gpart.LastDaughter());
      fGPartPx.push_back(gpart.Px());
      fGPartPy.push_back(gpart.Py());
      fGPartPz.push_back(gpart.Pz());
      fGPartE.push_back(gpart.E());
      fGPartMass.push_back(gpart.Mass());
      nGPart++;
    }
    fnGPart.push_back(nGPart);
  }

  auto MCPHandle = e.getHandle<std::vector<simb::MCParticle> >(fGeantLabel); 
  if (!MCPHandle) {
    throw cet::exception("anatree") << " No simb::MCParticle branch."
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
    fMCTrkID.push_back(mcp.TrackId());
    fMCPDG.push_back(mcp.PdgCode());

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
        MF_LOG_DEBUG("Anatree_module")
          << " mcp trkid " << mcp.TrackId()
          << " pdg code " << mcp.PdgCode()
          << " could not find mother trk id " << momTrkId
          << " in the TrackIdToIndex map"
          << " creating process is [ " << mcp.Process() << " ]";
      }
    }

    fMCMotherIndex.push_back(momIndex);
    fMCMotherTrkID.push_back(mcp.Mother());//directly trackid not index
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

  // Get vertex for each MCParticle.  In principle, the mct.GetParticle()
  // method should produce all of them; filter out the ones with a non-stable
  // status code and the GENIE pseudo-particles (Pid>= 2000000000) and you
  // should have the contents of (*MCPHandle).  But you don't.  On the 1st
  // event I (Leo Bellantoni) tested, with 40 vertices and 351 primary
  // particles, particle 143, a proton, evidently decayed or something and is
  // not in (*MCPHandle).  So we use the following primitive earth technology.
  size_t nMCParticles = (*MCPHandle).size();
  size_t iMCParticle  = 0;
  fMCPVertIndex.resize(nMCParticles);
  for (; iMCParticle<nMCParticles; ++iMCParticle) {
    // Assign noprimary to start with
    fMCPVertIndex[iMCParticle] = -1;
    // Do the primaries first
    if (fMCMotherIndex[iMCParticle]!=-1) break;
    Float_t trackX = fMCPStartX[iMCParticle];
    Float_t trackY = fMCPStartY[iMCParticle];
    Float_t trackZ = fMCPStartZ[iMCParticle];
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
            fMCPVertIndex[iMCParticle] = vertexIndex;
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
    int momIndex = fMCMotherIndex[iMCParticle];
    int lastMCParticle = iMCParticle;
    while (momIndex != -1) {
      lastMCParticle = momIndex;
      momIndex       = fMCMotherIndex[momIndex];
    }
    fMCPVertIndex[iMCParticle] = fMCPVertIndex[lastMCParticle];
  }

  if (fWriteMCPTrajectory) {
    // It's in the MCParticle table
    Int_t mcpIndex = 0;
    for ( auto const& mcp : (*MCPHandle) ) {
      const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
      const TParticlePDG* definition = databasePDG->GetParticle( mcp.PdgCode() );
      //No charge don't store the trajectory
      // this test fails for alpha particles because they aren't in databasePDG
      // so skip it in this case.
      if (mcp.PdgCode() != 1000020040)
        {
	  if (definition==nullptr || definition->Charge() == 0) continue;
	}
      //TrackID of the mcp to keep track to which mcp this trajectory is
      int trackId = mcp.TrackId();
      for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++) {
        float xTraj = mcp.Trajectory().X(iTraj);
        float yTraj = mcp.Trajectory().Y(iTraj);
        float zTraj = mcp.Trajectory().Z(iTraj);
        float rTraj = std::hypot( yTraj - ItsInTulsa[1], zTraj - ItsInTulsa[2]);

        if (abs(xTraj - ItsInTulsa[0]) > xTPC) continue;
        if (rTraj > rTPC) continue;

        fTrajMCPX.push_back(xTraj);
        fTrajMCPY.push_back(yTraj);
        fTrajMCPZ.push_back(zTraj);
        fTrajMCPT.push_back(mcp.Trajectory().T(iTraj));
        if (fWriteMCPTrajMomenta) {
          fTrajMCPPX.push_back(mcp.Trajectory().Px(iTraj));
          fTrajMCPPY.push_back(mcp.Trajectory().Py(iTraj));
          fTrajMCPPZ.push_back(mcp.Trajectory().Pz(iTraj));
        }
        fTrajMCPE.push_back(mcp.Trajectory().E(iTraj));
        fTrajMCPIndex.push_back(mcpIndex);
        fTrajMCPTrackID.push_back(trackId);
      }
      mcpIndex++;
    }
  }

  // Get handles for MCCaloInfo

  if (fWriteMCCaloInfo) {
    art::InputTag ecalgeanttag(fGeantLabel, fInstanceLabelCalo);
    auto SimHitHandle = e.getHandle<std::vector<gar::sdp::CaloDeposit> >(ecalgeanttag);
    if (!SimHitHandle) {
      throw cet::exception("anatree") << " No gar::sdp::CaloDeposit branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // Save simulation ecal hit info
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

    if(fWriteMuID) {
      art::InputTag muidgeanttag(fGeantLabel, fInstanceLabelMuID);
      art::Handle< std::vector<gar::sdp::CaloDeposit> > MuIDSimHitHandle;
      if (fGeo->HasMuonDetector())
        {
          MuIDSimHitHandle = e.getHandle< std::vector<gar::sdp::CaloDeposit> >(muidgeanttag);
          if (!MuIDSimHitHandle) {
            throw cet::exception("anatree") << " No gar::sdp::CaloDeposit branch."
                                            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
          }
        }

      // Save simulation muon system hit info
      for ( auto const& SimHit : (*MuIDSimHitHandle) ) {
        fSimnHits_MuID++;
        fSimHitX_MuID.push_back(SimHit.X());
        fSimHitY_MuID.push_back(SimHit.Y());
        fSimHitZ_MuID.push_back(SimHit.Z());
        fSimHitTime_MuID.push_back(SimHit.Time());
        fSimHitEnergy_MuID.push_back(SimHit.Energy());
        fSimHitTrackID_MuID.push_back(SimHit.TrackID());
        fSimHitCellID_MuID.push_back(SimHit.CellID());
        fSimEnergySum_MuID += SimHit.Energy();
      }
    }
  }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillRawInfo(art::Event const & e) {

  // save ecal raw digits info
  art::InputTag ecalrawtag(fRawCaloHitLabel, fInstanceLabelCalo);
  auto RawHitHandle = e.getHandle<std::vector<gar::raw::CaloRawDigit> >(ecalrawtag);
  if (!RawHitHandle) {
    throw cet::exception("anatree") << " No :raw::CaloRawDigit branch."
                                    << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
  }
  for ( auto const& DigiHit : (*RawHitHandle) ) {
    fDiginHits++;
    fDigiHitX.push_back(DigiHit.X());
    fDigiHitY.push_back(DigiHit.Y());
    fDigiHitZ.push_back(DigiHit.Z());
    fDigiHitTime.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
    fDigiHitADC.push_back(DigiHit.ADC().first);
    fDigiHitCellID.push_back(DigiHit.CellID());
  }

  // save muon system raw digits info
  if(fWriteMuID) {
    art::InputTag muidrawtag(fRawMuIDHitLabel, fInstanceLabelMuID);
    art::Handle<std::vector<gar::raw::CaloRawDigit> > MuIDRawHitHandle;
    if (fGeo->HasMuonDetector()) {
      MuIDRawHitHandle = e.getHandle<std::vector<gar::raw::CaloRawDigit> >(muidrawtag);
      if (!MuIDRawHitHandle) {
        throw cet::exception("anatree") << " No :raw::CaloRawDigit branch."
                                        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }

    for ( auto const& DigiHit : (*MuIDRawHitHandle) ) {
      fDiginHits_MuID++;
      fDigiHitX_MuID.push_back(DigiHit.X());
      fDigiHitY_MuID.push_back(DigiHit.Y());
      fDigiHitZ_MuID.push_back(DigiHit.Z());
      fDigiHitTime_MuID.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2.0 );
      fDigiHitADC_MuID.push_back(DigiHit.ADC().first);
      fDigiHitCellID_MuID.push_back(DigiHit.CellID());
    }
  }
}



//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillRecoInfo(art::Event const & e) {

  // Get handle for TPC hit data
  if (fWriteHits) {
    auto HitHandle = e.getHandle<std::vector<rec::Hit> >(fHitLabel);
    if (!HitHandle) {
      throw cet::exception("anatree") << " No rec::Hit branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    // save hits in the TPC
    for ( auto const& Hit : (*HitHandle) ) {
      fHitX.push_back(Hit.Position()[0]);
      fHitY.push_back(Hit.Position()[1]);
      fHitZ.push_back(Hit.Position()[2]);
      fHitSig.push_back(Hit.Signal());
      fHitRMS.push_back(Hit.RMS());
      fHitChan.push_back(Hit.Channel());
    }
  }

  if (fWriteCaloHits) {
    art::InputTag ecalrecotag(fCaloHitLabel, fInstanceLabelCalo);
    auto RecoHitHandle = e.getHandle<std::vector<rec::CaloHit> >(ecalrecotag);
    if (!RecoHitHandle) {
      throw cet::exception("anatree") << " No rec::CaloHit branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    // save reco'd Calorimetry hits
    for ( auto const& Hit : (*RecoHitHandle) ) {
      fReconHits++;
      fReconHitIDNumber.push_back(Hit.getIDNumber());
      fRecoHitX.push_back(Hit.Position()[0]);
      fRecoHitY.push_back(Hit.Position()[1]);
      fRecoHitZ.push_back(Hit.Position()[2]);
      fRecoHitTime.push_back(Hit.Time().first);
      fRecoHitEnergy.push_back(Hit.Energy());
      fRecoHitCellID.push_back(Hit.CellID());
      fRecoHitLayer.push_back(fFieldDecoder_ECAL->get(Hit.CellID(), "layer"));
      fRecoEnergySum += Hit.Energy();
    }

    if(fWriteMuID) {
      art::InputTag muirecotag(fMuIDHitLabel, fInstanceLabelMuID);
      art::Handle<std::vector<rec::CaloHit> > MuIDRecoHitHandle;
      if (fGeo->HasMuonDetector()) {
        MuIDRecoHitHandle = e.getHandle<std::vector<rec::CaloHit> >(muirecotag);
        if (!MuIDRecoHitHandle) {
          throw cet::exception("anatree") << " No rec::CaloHit branch."
                                          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
      }
      for ( auto const& Hit : (*MuIDRecoHitHandle) ) {
        fReconHits_MuID++;
        fReconHitIDNumber_MuID.push_back(Hit.getIDNumber());
        fRecoHitX_MuID.push_back(Hit.Position()[0]);
        fRecoHitY_MuID.push_back(Hit.Position()[1]);
        fRecoHitZ_MuID.push_back(Hit.Position()[2]);
        fRecoHitTime_MuID.push_back(Hit.Time().first);
        fRecoHitEnergy_MuID.push_back(Hit.Energy());
        fRecoHitCellID_MuID.push_back(Hit.CellID());
        fRecoHitLayer_MuID.push_back(fFieldDecoder_MuID->get(Hit.CellID(), "layer"));
        fRecoEnergySum_MuID += Hit.Energy();
      }
    }
  }
}

//==============================================================================
//==============================================================================
//==============================================================================
void gar::anatree::FillHighLevelRecoInfo(art::Event const & e) {

  // Get handle for TPCClusters
  art::Handle< std::vector<rec::TPCCluster> > TPCClusterHandle;
  art::FindManyP<rec::Hit>* findManyHits = NULL;
  if (fWriteTPCClusters) {
    TPCClusterHandle = e.getHandle< std::vector<rec::TPCCluster> >(fTPCClusterLabel);
    if (!TPCClusterHandle) {
      throw cet::exception("anatree") << " No rec::TPCCluster branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    findManyHits = new art::FindManyP<rec::Hit>(TPCClusterHandle,e,fTPCClusterLabel);
  }

  // Get handles for Tracks and their ionizations; also Assn's to TPCClusters, TrackIoniz
  // null handles if we switch off reading in the data products
  art::Handle< std::vector<rec::Track> > TrackHandle;
  art::Handle< std::vector<rec::TrackIoniz> > TrackIonHandle;
  art::Handle< std::vector<rec::TrackTrajectory> > TrackTrajHandle;
  art::FindManyP<rec::TPCCluster>* findManyTPCClusters = NULL;
  art::FindOneP<rec::TrackIoniz>*  findIonization = NULL;
  if(fWriteTracks) {
    TrackHandle = e.getHandle< std::vector<rec::Track> >(fTrackLabel);
    if (!TrackHandle) {
      throw cet::exception("anatree") << " No rec::Track branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }
    TrackIonHandle = e.getHandle< std::vector<rec::TrackIoniz> >(fTrackLabel);
    if (!TrackIonHandle) {
      throw cet::exception("anatree") << " No rec::TrackIoniz branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    findManyTPCClusters = new art::FindManyP<rec::TPCCluster>(TrackHandle,e,fTrackLabel);
    findIonization      = new art::FindOneP<rec::TrackIoniz>(TrackHandle,e,fTrackLabel);

    if(fWriteTrackTrajectories) {
      TrackTrajHandle = e.getHandle< std::vector<rec::TrackTrajectory> >(fTrackTrajectoryLabel);
      if (!TrackTrajHandle) {
        throw cet::exception("anatree") << " No rec::TrackTrajectory branch."
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
      throw cet::exception("anatree") << " No rec::Vertex branch."
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
      throw cet::exception("anatree") << " No rec::Vee branch."
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
  art::FindManyP<gar::rec::CaloHit>* findManyClusterRecoHit = NULL;
  art::FindManyP<gar::rec::CaloHit>* findManyClusterMuIDHit = NULL;

  if (fWriteCaloClusters) {
    RecoClusterHandle = e.getHandle< std::vector<rec::Cluster> >(ecalclustertag);
    if (!RecoClusterHandle) {
      throw cet::exception("anatree") << " No rec::Cluster branch."
                                      << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

    if(fWriteMuID) {
      if (fGeo->HasMuonDetector()){
        RecoClusterMuIDHandle = e.getHandle< std::vector<rec::Cluster> >(muidclustertag);
        if (!RecoClusterMuIDHandle) {
          throw cet::exception("anatree") << " No rec::Cluster (MuID) branch."
                                          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
      }
    }

    findManyClusterRecoHit = new art::FindManyP<gar::rec::CaloHit>(RecoClusterHandle,e,ecalclustertag);

    if (fGeo->HasMuonDetector() && fWriteMuID)
      findManyClusterMuIDHit = new art::FindManyP<gar::rec::CaloHit>(RecoClusterMuIDHandle,e,muidclustertag);
        

    if (fWriteTracks)
      findManyCALTrackEnd = new art::FindManyP<rec::Track, rec::TrackEnd>(RecoClusterHandle,e,fECALAssnLabel);

  }

  // save clusters in the TPC. For some reason, can't get FindOneP<rec::Track> or
  // FindManyP<rec::Track> to work; seems the underlying Assn isn't found.  Have
  // to FindManyP<TPCCluster> instead and  iterate if (fWriteTracks).  :(
  if (fWriteTPCClusters) {
    size_t iTPCCluster = 0;
    for ( auto const& TPCCluster : (*TPCClusterHandle) ) {
      fTPCClusterX.push_back(TPCCluster.Position()[0]);
      fTPCClusterY.push_back(TPCCluster.Position()[1]);
      fTPCClusterZ.push_back(TPCCluster.Position()[2]);
      fTPCClusterSig.push_back(TPCCluster.Signal());
      fTPCClusterRMS.push_back(TPCCluster.RMS());
      const float* cov;
      cov = TPCCluster.CovMatPacked();
      fTPCClusterCovXX.push_back(cov[0]);
      fTPCClusterCovXY.push_back(cov[1]);
      fTPCClusterCovXZ.push_back(cov[2]);
      fTPCClusterCovYY.push_back(cov[3]);
      fTPCClusterCovYZ.push_back(cov[4]);
      fTPCClusterCovZZ.push_back(cov[5]);

      // To get MC matching info from TPCClusters, 1st get the associated Hits
      int indexToPush = -1;        float valueToPush = 0;
      if ( findManyHits->isValid() ) {
        std::map<int,float> sumEforTrkID;
        float eTotCluster = 0;
        auto const& hitsInTPCCluster = findManyHits->at(iTPCCluster);
        for (size_t iHits = 0; iHits<hitsInTPCCluster.size(); ++iHits) {
          std::vector<cheat::HitIDE> IDEs = BackTrack->HitToHitIDEs( hitsInTPCCluster[iHits] );
          for (size_t iIDE = 0; iIDE<IDEs.size(); ++iIDE) {
            int  trackID    = IDEs[iIDE].trackID;
            float thisEdepE = IDEs[iIDE].energyTot;
            if ( sumEforTrkID.find(trackID) == sumEforTrkID.end() ) {
              sumEforTrkID[trackID] = 0;
            }
            sumEforTrkID[trackID] += thisEdepE;
            eTotCluster           += thisEdepE;
          }
        }
        if (sumEforTrkID.size()!=0) {
          // Sort the map by value.  Start by declaring the type of the sorting predicate
          typedef std::function<bool(std::pair<int,float>, std::pair<int,float>)> Comparator;
          // Declare a set that will store the pairs using above comparison logic
          std::set<std::pair<int,float>, Comparator> setOfTrkIDs(
                                                                 sumEforTrkID.begin(), sumEforTrkID.end(),
                                                                 [](std::pair<int,float> a ,std::pair<int,float> b) {
                                                                   return a.second > b.second;
                                                                 }
                                                                 );
          auto iReturnSet = setOfTrkIDs.begin();
          indexToPush     = TrackIdToIndex[iReturnSet->first];
          valueToPush     = iReturnSet->second/eTotCluster;
        }
      }
      fTPCClusterMCindex.push_back(indexToPush);
      fTPCClusterMCfrac.push_back(valueToPush);

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
      iTPCCluster++;
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
      fTrackChi2F.push_back(track.ChisqForward());
      fTrackChi2B.push_back(track.ChisqBackward());
      fNTPCClustersOnTrack.push_back(track.NHits());

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
        fTrackPIDF.push_back( pidF.at(ipid).first );
        fTrackPIDProbF.push_back( pidF.at(ipid).second );
      }
      for(size_t ipid = 0; ipid < pidB.size(); ipid++) {
        fTrackPIDB.push_back( pidB.at(ipid).first );
        fTrackPIDProbB.push_back( pidB.at(ipid).second );
      }

      // Matching MCParticle info
      std::vector<std::pair<simb::MCParticle*,float>> trakt;
      trakt = BackTrack->TrackToMCParticles( const_cast<rec::Track*>(&track) );
      int eileen = -1;
      if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
        eileen = TrackIdToIndex[trakt[0].first->TrackId()];
      }
      fTrackMCindex.push_back(eileen);
      if (eileen > -1) {
        fTrackMCfrac.push_back(trakt[0].second);
      } else {
        fTrackMCfrac.push_back(0.0);
      }

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

  //TrackTrajectories
  if(fWriteTrackTrajectories) {
    size_t iTrackTraj = 0;
    for ( auto const& tracktraj : (*TrackTrajHandle) ) {
            
      std::vector<TVector3> temp = tracktraj.getFWDTrajectory();
      for(size_t i = 0; i < temp.size(); i++) {
        fTrackTrajectoryFWDX.push_back(temp.at(i).X());
        fTrackTrajectoryFWDY.push_back(temp.at(i).Y());
        fTrackTrajectoryFWDZ.push_back(temp.at(i).Z());
        fTrackTrajectoryFWDID.push_back(iTrackTraj);
      }

      temp = tracktraj.getBAKTrajectory();
      for(size_t i = 0; i < temp.size(); i++) {
        fTrackTrajectoryBWDX.push_back(temp.at(i).X());
        fTrackTrajectoryBWDY.push_back(temp.at(i).Y());
        fTrackTrajectoryBWDZ.push_back(temp.at(i).Z());
        fTrackTrajectoryBWDID.push_back(iTrackTraj);
      }
      iTrackTraj++;
    }
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

  // save Vee and Track-Vee association info
  if (fWriteVees) {
    size_t iVee = 0;
    for ( auto const& vee : (*VeeHandle) ) {
      fVeeIDNumber.push_back(vee.getIDNumber());
      fVeeX.push_back(vee.Position()[0]);
      fVeeY.push_back(vee.Position()[1]);
      fVeeZ.push_back(vee.Position()[2]);
      fVeeT.push_back(vee.Time());
      fVeePXKpipi.push_back(vee.FourMomentum(0).X());
      fVeePYKpipi.push_back(vee.FourMomentum(0).Y());
      fVeePZKpipi.push_back(vee.FourMomentum(0).Z());
      fVeeEKpipi.push_back(vee.FourMomentum(0).E());
      fVeeMKpipi.push_back(vee.FourMomentum(0).M());
      fVeePXLppi.push_back(vee.FourMomentum(1).X());
      fVeePYLppi.push_back(vee.FourMomentum(1).Y());
      fVeePZLppi.push_back(vee.FourMomentum(1).Z());
      fVeeELppi.push_back(vee.FourMomentum(1).E());
      fVeeMLppi.push_back(vee.FourMomentum(1).M());
      fVeePXLpip.push_back(vee.FourMomentum(2).X());
      fVeePYLpip.push_back(vee.FourMomentum(2).Y());
      fVeePZLpip.push_back(vee.FourMomentum(2).Z());
      fVeeELpip.push_back(vee.FourMomentum(2).E());
      fVeeMLpip.push_back(vee.FourMomentum(2).M());

      int nVeeTracks = 0;
      if ( findManyVeeTrackEnd->isValid() ) {
        nVeeTracks = findManyVeeTrackEnd->at(iVee).size();
      }

      for (int iVeeTrack=0; iVeeTrack<nVeeTracks; ++iVeeTrack) {
        fVeeTAssn_VeeIDNumber.push_back(vee.getIDNumber());

        // Get this vertexed track.
        rec::Track track = *(findManyVeeTrackEnd->at(iVee).at(iVeeTrack));
        fVeeTAssn_TrackIDNumber.push_back(track.getIDNumber());

        rec::TrackEnd fee = *(findManyVeeTrackEnd->data(iVee).at(iVeeTrack));
        // TrackEnd is defined in Track.h; 1 means use Beg values, 0 means use End
        fVeeTAssn_TrackEnd.push_back(fee);
      }
      ++iVee;
    } // end loop over VeeHandle
  }

  // save Cluster info
  if (fWriteCaloClusters) {
    size_t iCluster = 0;
    for ( auto const& cluster : (*RecoClusterHandle) ) {
      fnCluster++;
      fClusterIDNumber.push_back(cluster.getIDNumber());
      fClusterNhits.push_back(cluster.CalorimeterHits().size());
      fClusterEnergy.push_back(cluster.Energy());
      fClusterTime.push_back(cluster.Time());
      fClusterTimeDiffFirstLast.push_back(cluster.TimeDiffFirstLast());
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

      //Get the associated reco hit
      std::vector<ULong64_t> fVecHitIDs = {};
      if (findManyClusterRecoHit->isValid()) {
        int nClusterHit = findManyClusterRecoHit->at(iCluster).size();
        for (int iClusterHit=0; iClusterHit<nClusterHit; ++iClusterHit) {
          rec::CaloHit hit  = *(findManyClusterRecoHit->at(iCluster).at(iClusterHit));
          fVecHitIDs.push_back(hit.getIDNumber());
        }
      }

      fClusterAssn_RecoHitIDNumber.push_back(fVecHitIDs);

      // Matching MCParticle info
      std::vector<std::pair<simb::MCParticle*,float>> trakt;
      trakt = BackTrack->ClusterToMCParticles( const_cast<rec::Cluster*>(&cluster) );
      int eileen = -1;
      if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
        eileen = TrackIdToIndex[trakt[0].first->TrackId()];
      }
      fClusterMCindex.push_back(eileen);
      if (eileen > -1) {
        fClusterMCfrac.push_back(trakt[0].second);
      } else {
        fClusterMCfrac.push_back(0.0);
      }
      iCluster++;
    }

    if(fGeo->HasMuonDetector() && fWriteMuID) {
      size_t iCluster_local = 0;
      for ( auto const& cluster : (*RecoClusterMuIDHandle) ) {
        fnCluster_MuID++;
        fClusterIDNumber_MuID.push_back(cluster.getIDNumber());
        fClusterNhits_MuID.push_back(cluster.CalorimeterHits().size());
        fClusterEnergy_MuID.push_back(cluster.Energy());
        fClusterTime_MuID.push_back(cluster.Time());
        fClusterTimeDiffFirstLast_MuID.push_back(cluster.TimeDiffFirstLast());
        fClusterX_MuID.push_back(cluster.Position()[0]);
        fClusterY_MuID.push_back(cluster.Position()[1]);
        fClusterZ_MuID.push_back(cluster.Position()[2]);
        fClusterTheta_MuID.push_back(cluster.ITheta());
        fClusterPhi_MuID.push_back(cluster.IPhi());
        fClusterPID_MuID.push_back(cluster.ParticleID());
        // fClusterShape.push_back(cluster.Shape());
        fClusterMainAxisX_MuID.push_back(cluster.EigenVectors()[0]);
        fClusterMainAxisY_MuID.push_back(cluster.EigenVectors()[1]);
        fClusterMainAxisZ_MuID.push_back(cluster.EigenVectors()[2]);

        //Get the associated reco hit
        std::vector<ULong64_t> fVecHitIDs = {};
        if (findManyClusterMuIDHit->isValid()) {
          int nClusterHit = findManyClusterMuIDHit->at(iCluster_local).size();
          for (int iClusterHit=0; iClusterHit<nClusterHit; ++iClusterHit) {
            rec::CaloHit hit  = *(findManyClusterMuIDHit->at(iCluster_local).at(iClusterHit));
            fVecHitIDs.push_back(hit.getIDNumber());
          }
        }

        fClusterMuIDAssn_MuIDHitIDNumber.push_back(fVecHitIDs);

        /*** TO DO ***/
        // Matching MCParticle info
        // std::vector<std::pair<simb::MCParticle*,float>> trakt;
        // trakt = BackTrack->ClusterToMCParticles( const_cast<rec::Cluster*>(&cluster) );
        // int eileen = -1;
        // if (trakt.size()>0 && TrackIdToIndex.size()!=0) {
        //     eileen = TrackIdToIndex[trakt[0].first->TrackId()];
        // }
        // fClusterMCindex_MuID.push_back(eileen);
        // if (eileen > -1) {
        //     fClusterMCfrac_MuID.push_back(trakt[0].second);
        // } else {
        //     fClusterMCfrac_MuID.push_back(0.0);
        // }
        fClusterMCindex_MuID.push_back(0);
        fClusterMCfrac_MuID.push_back(0.0);
        iCluster_local++;
      }
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

  return;
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



//==============================================================================
//==============================================================================
//==============================================================================
std::vector< std::pair<int, float> > gar::anatree::processPIDInfo( float p ) {

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

DEFINE_ART_MODULE(gar::anatree)
