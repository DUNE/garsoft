////////////////////////////////////////////////////////////////////////
// Class:       CaloAnaTree
// Plugin Type: analyzer (art v2_11_02)
// File:        CaloAnaTree_module.cc
//
// Generated at Mon Oct 8 2018 by Eldwan Brianne
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"

#include "TTree.h"

#include <vector>

namespace gar
{

    class CaloAnaTree : public art::EDAnalyzer {

    public:
        explicit CaloAnaTree(fhicl::ParameterSet const & p);
        // The compiler-generated destructor is fine for non-base
        // classes without bare pointers or other resource use.

        // Plugins should not be copied or assigned.
        CaloAnaTree(CaloAnaTree const &) = delete;
        CaloAnaTree(CaloAnaTree &&) = delete;
        CaloAnaTree & operator = (CaloAnaTree const &) = delete;
        CaloAnaTree & operator = (CaloAnaTree &&) = delete;

        virtual void beginJob() override;

        // Required functions.
        void analyze(art::Event const & e) override;

        void ClearVectors();

        void FillVectors(art::Event const & e);

    private:

        // Input data labels

        std::string fGeneratorLabel;
        std::string fGeantLabel;
        std::string fRawCaloHitLabel;
        std::string fCaloHitLabel;
        std::string fClusterLabel;

        // the analysis tree

        TTree *fTree;

        // global event info

        int fEvent;        ///< number of the event being processed
        int fRun;          ///< number of the run being processed
        int fSubRun;       ///< number of the sub-run being processed

        // MCTruth data

        std::vector<int> fNeutrinoType;
        std::vector<int> fCCNC;

        // MCParticle data

        std::vector<int>   fMCPPDGID;
        std::vector<float> fMCPStartX;
        std::vector<float> fMCPStartY;
        std::vector<float> fMCPStartZ;
        std::vector<float> fMCPEndX;
        std::vector<float> fMCPEndY;
        std::vector<float> fMCPEndZ;
        std::vector<float> fMCPPX;
        std::vector<float> fMCPPY;
        std::vector<float> fMCPPZ;
        std::vector<std::string> fMCPProc;
        std::vector<std::string> fMCPEndProc;
        std::vector<bool> fMCPPrimaryConv;

        // sim calo hit data

        std::vector<float> fSimHitX;
        std::vector<float> fSimHitY;
        std::vector<float> fSimHitZ;
        std::vector<float> fSimHitTime;
        std::vector<float> fSimHitEnergy;
        std::vector<int> fSimHitTrackID;
        std::vector<long long int> fSimHitCellID;
        float fSimEnergySum;
        unsigned int fSimnHits;

        // // raw calo hit data

        std::vector<float> fDigiHitX;
        std::vector<float> fDigiHitY;
        std::vector<float> fDigiHitZ;
        std::vector<float> fDigiHitTime;
        std::vector<unsigned int> fDigiHitADC;
        std::vector<long long int> fDigiHitCellID;
        unsigned int fDiginHits;

        // reco calo hit data

        std::vector<float> fRecoHitX;
        std::vector<float> fRecoHitY;
        std::vector<float> fRecoHitZ;
        std::vector<float> fRecoHitTime;
        std::vector<float> fRecoHitEnergy;
        std::vector<long long int> fRecoHitCellID;
        float fRecoEnergySum;
        unsigned int fReconHits;

        // calo cluster data
        unsigned int fnCluster;
        std::vector<unsigned int> fClusterNhits;
        std::vector<float> fClusterEnergy;
        std::vector<float> fClusterCoGX;
        std::vector<float> fClusterCoGY;
        std::vector<float> fClusterCoGZ;
        std::vector<float> fClusterDirectionX;
        std::vector<float> fClusterDirectionY;
        std::vector<float> fClusterDirectionZ;
        std::vector<float> fClusterMainAxisX;
        std::vector<float> fClusterMainAxisY;
        std::vector<float> fClusterMainAxisZ;
        std::vector<float> fClusterSecondAxisX;
        std::vector<float> fClusterSecondAxisY;
        std::vector<float> fClusterSecondAxisZ;
        std::vector<float> fClusterThirdAxisX;
        std::vector<float> fClusterThirdAxisY;
        std::vector<float> fClusterThirdAxisZ;

    };

}

// constructor

gar::CaloAnaTree::CaloAnaTree(fhicl::ParameterSet const & p)
:
    EDAnalyzer(p)  // ,
    // More initializers here.

    {
        fGeneratorLabel = p.get<std::string>("GeneratorLabel","generator");
        fGeantLabel     = p.get<std::string>("GEANTLabel","geant");
        fRawCaloHitLabel       = p.get<std::string>("RawCaloHitLabel","daqecal");
        fCaloHitLabel     = p.get<std::string>("CaloHitLabel","calohit");
        fClusterLabel     = p.get<std::string>("ClusterLabel","cluster");

        consumes<std::vector<simb::MCParticle> >(fGeantLabel);
        consumes<std::vector<simb::MCTruth> >(fGeneratorLabel);
        //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
        consumes<std::vector<gar::sdp::CaloDeposit> >(fGeantLabel);
        consumes<std::vector<gar::raw::CaloRawDigit> >(fRawCaloHitLabel);
        consumes<std::vector<gar::rec::CaloHit> >(fCaloHitLabel);
        consumes<std::vector<gar::rec::Cluster> >(fClusterLabel);
    }

    void gar::CaloAnaTree::beginJob()
    {
        art::ServiceHandle<art::TFileService> tfs;

        //Make the TTree
        fTree = tfs->make<TTree>("ECALAnaTree","ECALAnaTree");

        fTree->Branch("Event",       &fEvent,          "Event/I");
        fTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
        fTree->Branch("Run",         &fRun,            "Run/I");

        fTree->Branch("NType", &fNeutrinoType);
        fTree->Branch("CCNC",  &fCCNC);

        fTree->Branch("PDG", &fMCPPDGID);
        fTree->Branch("MCPStartX", &fMCPStartX);
        fTree->Branch("MCPStartY", &fMCPStartY);
        fTree->Branch("MCPStartZ", &fMCPStartZ);
        fTree->Branch("MCPEndX", &fMCPEndX);
        fTree->Branch("MCPEndY", &fMCPEndY);
        fTree->Branch("MCPEndZ", &fMCPEndZ);
        fTree->Branch("MCPPX", &fMCPPX);
        fTree->Branch("MCPPY", &fMCPPY);
        fTree->Branch("MCPPZ", &fMCPPZ);
        fTree->Branch("MCPProc",  &fMCPProc);
        fTree->Branch("MCPEndProc", &fMCPEndProc);

        fTree->Branch("SimnHits", &fSimnHits);
        fTree->Branch("SimHitX", &fSimHitX);
        fTree->Branch("SimHitY", &fSimHitY);
        fTree->Branch("SimHitZ", &fSimHitZ);
        fTree->Branch("SimHitTime", &fSimHitTime);
        fTree->Branch("SimHitEnergy", &fSimHitEnergy);
        fTree->Branch("SimHitTrkID", &fSimHitTrackID);
        fTree->Branch("SimHitCellID", &fSimHitCellID);
        fTree->Branch("SimEnergySum", &fSimEnergySum);

        fTree->Branch("DiginHits", &fDiginHits);
        fTree->Branch("DigiHitX", &fDigiHitX);
        fTree->Branch("DigiHitY", &fDigiHitY);
        fTree->Branch("DigiHitZ", &fDigiHitZ);
        fTree->Branch("DigiHitTime", &fDigiHitTime);
        fTree->Branch("DigiHitEnergy", &fDigiHitADC);
        fTree->Branch("DigiHitCellID", &fDigiHitCellID);

        fTree->Branch("ReconHits", &fReconHits);
        fTree->Branch("RecoHitX", &fRecoHitX);
        fTree->Branch("RecoHitY", &fRecoHitY);
        fTree->Branch("RecoHitZ", &fRecoHitZ);
        fTree->Branch("RecoHitTime", &fRecoHitTime);
        fTree->Branch("RecoHitEnergy", &fRecoHitEnergy);
        fTree->Branch("RecoHitCellID", &fRecoHitCellID);
        fTree->Branch("RecoEnergySum", &fRecoEnergySum);

        fTree->Branch("nCluster", &fnCluster);
        fTree->Branch("ClusterNhits", &fClusterNhits);
        fTree->Branch("ClusterEnergy", &fClusterEnergy);
        fTree->Branch("ClusterCoGX", &fClusterCoGX);
        fTree->Branch("ClusterCoGY", &fClusterCoGY);
        fTree->Branch("ClusterCoGZ", &fClusterCoGZ);
        fTree->Branch("ClusterDirectionX", &fClusterDirectionX);
        fTree->Branch("ClusterDirectionY", &fClusterDirectionY);
        fTree->Branch("ClusterDirectionZ", &fClusterDirectionZ);
        fTree->Branch("ClusterMainAxisX", &fClusterMainAxisX);
        fTree->Branch("ClusterMainAxisY", &fClusterMainAxisY);
        fTree->Branch("ClusterMainAxisZ", &fClusterMainAxisZ);
        fTree->Branch("ClusterSecondAxisX", &fClusterSecondAxisX);
        fTree->Branch("ClusterSecondAxisY", &fClusterSecondAxisY);
        fTree->Branch("ClusterSecondAxisZ", &fClusterSecondAxisZ);
        fTree->Branch("ClusterThirdAxisX", &fClusterThirdAxisX);
        fTree->Branch("ClusterThirdAxisY", &fClusterThirdAxisY);
        fTree->Branch("ClusterThirdAxisZ", &fClusterThirdAxisZ);
    }


    void gar::CaloAnaTree::analyze(art::Event const & e)
    {
        // clear out all our vectors
        this->ClearVectors();

        // Fill the vectors
        this->FillVectors(e);

        fTree->Fill();
    }

    void gar::CaloAnaTree::ClearVectors()
    {
        fNeutrinoType.clear();
        fCCNC.clear();

        fMCPPDGID.clear();
        fMCPStartX.clear();
        fMCPStartY.clear();
        fMCPStartZ.clear();
        fMCPEndX.clear();
        fMCPEndY.clear();
        fMCPEndZ.clear();
        fMCPPX.clear();
        fMCPPY.clear();
        fMCPPZ.clear();
        fMCPProc.clear();
        fMCPEndProc.clear();

        // sim calo hit data

        fSimnHits = 0;
        fSimHitX.clear();
        fSimHitY.clear();
        fSimHitZ.clear();
        fSimHitTime.clear();
        fSimHitEnergy.clear();
        fSimHitTrackID.clear();
        fSimHitCellID.clear();
        fSimEnergySum = 0.;

        // raw calo hit data

        fDiginHits = 0;
        fDigiHitX.clear();
        fDigiHitY.clear();
        fDigiHitZ.clear();
        fDigiHitTime.clear();
        fDigiHitADC.clear();
        fDigiHitCellID.clear();

        // reco calo hit data

        fReconHits = 0;
        fRecoHitX.clear();
        fRecoHitY.clear();
        fRecoHitZ.clear();
        fRecoHitTime.clear();
        fRecoHitEnergy.clear();
        fRecoHitCellID.clear();
        fRecoEnergySum = 0.;

        // calo cluster data

        fnCluster = 0;
        fClusterNhits.clear();
        fClusterEnergy.clear();
        fClusterCoGX.clear();
        fClusterCoGY.clear();
        fClusterCoGZ.clear();
        fClusterDirectionX.clear();
        fClusterDirectionY.clear();
        fClusterDirectionZ.clear();
        fClusterMainAxisX.clear();
        fClusterMainAxisY.clear();
        fClusterMainAxisZ.clear();
        fClusterSecondAxisX.clear();
        fClusterSecondAxisY.clear();
        fClusterSecondAxisZ.clear();
        fClusterThirdAxisX.clear();
        fClusterThirdAxisY.clear();
        fClusterThirdAxisZ.clear();
    }

    void gar::CaloAnaTree::FillVectors(art::Event const & e)
    {

        fEvent  = e.id().event();
        fRun    = e.run();
        fSubRun = e.subRun();

        art::Handle< std::vector<simb::MCTruth> > MCTHandle;
        if (!e.getByLabel(fGeneratorLabel, MCTHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No simb::MCTruth branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        art::Handle< std::vector<simb::MCParticle> > MCPHandle;
        if (!e.getByLabel(fGeantLabel, MCPHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No simb::MCParticle branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        art::Handle< std::vector<gar::sdp::CaloDeposit> > SimHitHandle;
        if (!e.getByLabel(fGeantLabel, SimHitHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No gar::sdp::CaloDeposit branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        art::Handle< std::vector<gar::raw::CaloRawDigit> > RawHitHandle;
        if (!e.getByLabel(fRawCaloHitLabel, RawHitHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No gar::raw::CaloRawDigit branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        art::Handle< std::vector<gar::rec::CaloHit> > RecoHitHandle;
        if (!e.getByLabel(fCaloHitLabel, RecoHitHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No gar::rec::CaloHit branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        art::Handle< std::vector<gar::rec::Cluster> > RecoClusterHandle;
        if (!e.getByLabel(fClusterLabel, RecoClusterHandle))
        {
            throw cet::exception("CaloAnaTree")
            << " No gar::rec::Cluster branch - "
            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        // save MCTruth info

        for ( auto const& mct : (*MCTHandle) )
        {
            fNeutrinoType.push_back(mct.GetNeutrino().Nu().PdgCode());
            fCCNC.push_back(mct.GetNeutrino().CCNC());
        }

        // save MCParticle info

        for ( auto const& mcp : (*MCPHandle) )
        {
            fMCPPDGID.push_back(mcp.PdgCode());
            const TLorentzVector& pos = mcp.Position(0);
            const TLorentzVector& mom = mcp.Momentum(0);
            fMCPStartX.push_back(pos.X());
            fMCPStartY.push_back(pos.Y());
            fMCPStartZ.push_back(pos.Z());
            fMCPEndX.push_back(mcp.EndX());
            fMCPEndY.push_back(mcp.EndY());
            fMCPEndZ.push_back(mcp.EndZ());
            fMCPPX.push_back(mom.Px());
            fMCPPY.push_back(mom.Py());
            fMCPPZ.push_back(mom.Pz());
            fMCPProc.push_back(mcp.Process());
            fMCPEndProc.push_back(mcp.EndProcess());
        }

        //Save Sim Hit info

        for ( auto const& SimHit : (*SimHitHandle) )
        {
            fSimHitX.push_back(SimHit.X());
            fSimHitY.push_back(SimHit.Y());
            fSimHitZ.push_back(SimHit.Z());
            fSimHitTime.push_back(SimHit.Time());
            fSimHitEnergy.push_back(SimHit.Energy());
            fSimHitTrackID.push_back(SimHit.TrackID());
            fSimHitCellID.push_back(SimHit.CellID());

            fSimEnergySum += SimHit.Energy();
            fSimnHits++;
        }

        //Save Digit Hit info

        for ( auto const& DigiHit : (*RawHitHandle) )
        {
            fDigiHitX.push_back(DigiHit.X());
            fDigiHitY.push_back(DigiHit.Y());
            fDigiHitZ.push_back(DigiHit.Z());
            fDigiHitTime.push_back( (DigiHit.Time().first + DigiHit.Time().second) / 2. );
            fDigiHitADC.push_back(DigiHit.ADC());
            fDigiHitCellID.push_back(DigiHit.CellID());

            fDiginHits++;
        }

        //Save Reco Hit info

        for ( auto const& Hit : (*RecoHitHandle) )
        {
            fRecoHitX.push_back(Hit.Position()[0]);
            fRecoHitY.push_back(Hit.Position()[1]);
            fRecoHitZ.push_back(Hit.Position()[2]);
            fRecoHitTime.push_back(Hit.Time());
            fRecoHitEnergy.push_back(Hit.Energy());
            fRecoHitCellID.push_back(Hit.CellID());

            fRecoEnergySum += Hit.Energy();
            fReconHits++;
        }

        // save Cluster info
        for ( auto const& Cluster : (*RecoClusterHandle) )
        {
            fnCluster++;

            fClusterNhits.push_back(Cluster.NCaloHits());
            fClusterEnergy.push_back(Cluster.Energy());
            fClusterCoGX.push_back(Cluster.CenterOfGravity()[0]);
            fClusterCoGY.push_back(Cluster.CenterOfGravity()[1]);
            fClusterCoGZ.push_back(Cluster.CenterOfGravity()[2]);

            fClusterDirectionX.push_back(Cluster.Direction().x());
            fClusterDirectionY.push_back(Cluster.Direction().y());
            fClusterDirectionZ.push_back(Cluster.Direction().z());

            fClusterMainAxisX.push_back(Cluster.EigenVectors()[0].x());
            fClusterMainAxisY.push_back(Cluster.EigenVectors()[0].y());
            fClusterMainAxisZ.push_back(Cluster.EigenVectors()[0].z());

            fClusterSecondAxisX.push_back(Cluster.EigenVectors()[1].x());
            fClusterSecondAxisY.push_back(Cluster.EigenVectors()[1].y());
            fClusterSecondAxisZ.push_back(Cluster.EigenVectors()[1].z());

            fClusterThirdAxisX.push_back(Cluster.EigenVectors()[2].x());
            fClusterThirdAxisY.push_back(Cluster.EigenVectors()[2].y());
            fClusterThirdAxisZ.push_back(Cluster.EigenVectors()[2].z());
        }
    }

    DEFINE_ART_MODULE(gar::CaloAnaTree)
