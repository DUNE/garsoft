////////////////////////////////////////////////////////////////////////
// Class:       anatree
// Plugin Type: analyzer (art v2_11_02)
// File:        anatree_module.cc
//
// Generated at Mon Aug 27 16:41:13 2018 by Thomas Junk using cetskelgen
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
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Vertex.h"

#include "TTree.h"

#include <vector>

namespace gar
{

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

    // Input data labels

    std::string fGeneratorLabel;
    std::string fGeantLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fVertexLabel;

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
    std::vector<float> fMCPPX;
    std::vector<float> fMCPPY;
    std::vector<float> fMCPPZ;

    // hit data

    std::vector<float> fHitX;
    std::vector<float> fHitY;
    std::vector<float> fHitZ;
    std::vector<float> fHitSignal;
    std::vector<float> fHitRMS;

    // track branches

    std::vector<float> fTrackStartX;
    std::vector<float> fTrackStartY;
    std::vector<float> fTrackStartZ;
    std::vector<float> fTrackStartPX;
    std::vector<float> fTrackStartPY;
    std::vector<float> fTrackStartPZ;

    std::vector<float> fTrackEndX;
    std::vector<float> fTrackEndY;
    std::vector<float> fTrackEndZ;
    std::vector<float> fTrackEndPX;
    std::vector<float> fTrackEndPY;
    std::vector<float> fTrackEndPZ;

    // vertex branches

    std::vector<float> fVertexX;
    std::vector<float> fVertexY;
    std::vector<float> fVertexZ;
    std::vector<int> fVertexNTracks;    
 
  };

}

// constructor

gar::anatree::anatree(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
  // More initializers here.

{   
  fGeneratorLabel = p.get<std::string>("GeneratorLabel","generator");
  fGeantLabel     = p.get<std::string>("GEANTLabel","geant");
  fHitLabel       = p.get<std::string>("HitLabel","hit");
  fTrackLabel     = p.get<std::string>("TrackLabel","track");
  fVertexLabel    = p.get<std::string>("VertexLabel","vertex");

  consumes<std::vector<simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<simb::MCTruth> >(fGeneratorLabel);
  //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<gar::rec::Hit> >(fHitLabel);
  consumes<std::vector<gar::rec::Track> >(fTrackLabel);
  consumes<std::vector<gar::rec::Vertex> >(fVertexLabel);
  consumes<art::Assns<gar::rec::Track, gar::rec::Vertex> >(fVertexLabel);
}

void gar::anatree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

  fTree->Branch("Event",       &fEvent,          "Event/I");
  fTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
  fTree->Branch("Run",         &fRun,            "Run/I");

  fTree->Branch("NType", &fNeutrinoType);
  fTree->Branch("CCNC",  &fCCNC);

  fTree->Branch("PDG", &fMCPPDGID);
  fTree->Branch("MCPStartX", &fMCPStartX);
  fTree->Branch("MCPStartY", &fMCPStartY);
  fTree->Branch("MCPStartZ", &fMCPStartZ);
  fTree->Branch("MCPPX", &fMCPPX);
  fTree->Branch("MCPPY", &fMCPPY);
  fTree->Branch("MCPPZ", &fMCPPZ);

  fTree->Branch("HitX", &fHitX);
  fTree->Branch("HitY", &fHitY);
  fTree->Branch("HitZ", &fHitZ);
  fTree->Branch("HitSig", &fHitSignal);
  fTree->Branch("HitRMS", &fHitRMS);

  fTree->Branch("TrackStartX",  &fTrackStartX);
  fTree->Branch("TrackStartY",  &fTrackStartY);
  fTree->Branch("TrackStartZ",  &fTrackStartZ);
  fTree->Branch("TrackStartPX", &fTrackStartPX);
  fTree->Branch("TrackStartPY", &fTrackStartPY);
  fTree->Branch("TrackStartPZ", &fTrackStartPZ);

  fTree->Branch("TrackEndX",  &fTrackEndX);
  fTree->Branch("TrackEndY",  &fTrackEndY);
  fTree->Branch("TrackEndZ",  &fTrackEndZ);
  fTree->Branch("TrackEndPX", &fTrackEndPX);
  fTree->Branch("TrackEndPY", &fTrackEndPY);
  fTree->Branch("TrackEndPZ", &fTrackEndPZ);


  fTree->Branch("VertX", &fVertexX);
  fTree->Branch("VertY", &fVertexY);
  fTree->Branch("VertZ", &fVertexZ);
  fTree->Branch("VertN", &fVertexNTracks);    
}


void gar::anatree::analyze(art::Event const & e)
{

  // clear out all our vectors

  fNeutrinoType.clear();
  fCCNC.clear();
  fMCPPDGID.clear();
  fMCPStartX.clear();
  fMCPStartY.clear();
  fMCPStartZ.clear();
  fMCPPX.clear();
  fMCPPY.clear();
  fMCPPZ.clear();
  fHitX.clear();
  fHitY.clear();
  fHitZ.clear();
  fHitSignal.clear();
  fHitRMS.clear();
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
  fVertexX.clear();
  fVertexY.clear();
  fVertexZ.clear();
  fVertexNTracks.clear();

  fEvent  = e.id().event(); 
  fRun    = e.run();
  fSubRun = e.subRun();

  art::Handle< std::vector<simb::MCTruth> > MCTHandle;
  if (!e.getByLabel(fGeneratorLabel, MCTHandle)) 
    {
      throw cet::exception("anatree") 
	<< " No simb::MCTruth branch - "
	<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<simb::MCParticle> > MCPHandle;
  if (!e.getByLabel(fGeantLabel, MCPHandle)) 
    {
      throw cet::exception("anatree") 
	<< " No simb::MCParticle branch - "
	<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<gar::rec::Hit> > HitHandle;
  if (!e.getByLabel(fHitLabel, HitHandle)) 
    {
      throw cet::exception("anatree") 
	<< " No gar::rec::Hit branch - "
	<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<gar::rec::Track> > TrackHandle;
  if (!e.getByLabel(fTrackLabel, TrackHandle)) 
    {
      throw cet::exception("anatree") 
	<< " No gar::rec::Track branch - "
	<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<gar::rec::Vertex> > VertexHandle;
  if (!e.getByLabel(fVertexLabel, VertexHandle)) 
    {
      throw cet::exception("anatree") 
	<< " No gar::rec::Vertex branch - "
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
      fMCPPX.push_back(mom.Px());
      fMCPPY.push_back(mom.Py());
      fMCPPZ.push_back(mom.Pz());
    }

  for ( auto const& hit : (*HitHandle) )
    {
      fHitX.push_back(hit.Position()[0]);
      fHitY.push_back(hit.Position()[1]);
      fHitZ.push_back(hit.Position()[2]);
      fHitSignal.push_back(hit.Signal());
      fHitRMS.push_back(hit.RMS());
    }

  for ( auto const& track : (*TrackHandle) )
    {
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
    }

  for ( auto const& vertex : (*VertexHandle) )
    {
      fVertexX.push_back(vertex.Position()[0]);
      fVertexY.push_back(vertex.Position()[1]);
      fVertexZ.push_back(vertex.Position()[2]);

      // count up the tracks belonging to this vertex
      int ntracks = 0;
      const art::FindManyP<gar::rec::Track> findManyTrack(VertexHandle,e,fTrackLabel);
      if ( findManyTrack.isValid() )
	{
	  ntracks = findManyTrack.size();
	}
      fVertexNTracks.push_back(ntracks);
    }

  fTree->Fill();
}


DEFINE_ART_MODULE(gar::anatree)
