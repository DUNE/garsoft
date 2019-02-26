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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/SimChannel.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "Ana/anautil.h"

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

    // Get T from Q, pion in MCTruth
    double computeT( simb::MCTruth theMCTruth );

    // Input data labels

    std::string fGeneratorLabel;
    std::string fGeantLabel;
    std::string fHitLabel;
    std::string fTrackLabel;
    std::string fVertexLabel;

    // Optionally keep/drop parts of the analysis tree

    bool fWriteMCinfo;        ///< Info from MCTruth, GTruth into tree.  Default=true
    bool fWriteHits;          ///< Write info about hits into tree.  Default=true
    bool fWriteCohInfo;       ///< Write variables for coherent pi analysis.  Default=true
    bool fWriteHitsInTracks;  ///< Write all hits that are used in a track with approprate association Default=false

    // the analysis tree

    TTree *fTree;

    // global event info

    Int_t fEvent;        ///< number of the event being processed
    Int_t fRun;          ///< number of the run being processed
    Int_t fSubRun;       ///< number of the sub-run being processed

    // MCTruth data
    // GENIE kinematics computed ignoring Fermi momentum and the off-shellness of the bound nucleon.
    // Get the right kinematics here, except T has to be computed
    // Use Rtypes.h here, as these data get used by root

    std::vector<Int_t>     fNeutrinoType;
    std::vector<Int_t>     fCCNC;
    std::vector<Int_t>     fMode;
    std::vector<Int_t>     fInteractionType;
    std::vector<Float_t>   fQ2;
    std::vector<Float_t>   fW;
    std::vector<Float_t>   fX;
    std::vector<Float_t>   fY;
    std::vector<Float_t>   fTheta;
    std::vector<Float_t>   fT;
    std::vector<Float_t>   fMCVertexX;
    std::vector<Float_t>   fMCVertexY;
    std::vector<Float_t>   fMCVertexZ;
    std::vector<Float_t>   fMCnuPx;
    std::vector<Float_t>   fMCnuPy;
    std::vector<Float_t>   fMCnuPz;

    // GTruth data
    std::vector<Int_t>     fGint;
    std::vector<Int_t>     fTgtPDG;
    std::vector<Float_t>   fWeight;
    std::vector<Float_t>   fgT;

    // MCParticle data
    std::vector<Int_t>   fMCPPDGID;
    std::vector<Int_t>   fMCPDGMom;
    std::vector<Float_t> fMCPStartX;
    std::vector<Float_t> fMCPStartY;
    std::vector<Float_t> fMCPStartZ;
    std::vector<Float_t> fMCPPX;
    std::vector<Float_t> fMCPPY;
    std::vector<Float_t> fMCPPZ;

    // true trajectory hits data
    std::vector<Float_t> fTrajHitX;
    std::vector<Float_t> fTrajHitY;
    std::vector<Float_t> fTrajHitZ;
    std::vector<Int_t>   fTrajHitTrajIndex;

    // hit data
    std::vector<Float_t> fHitX;
    std::vector<Float_t> fHitY;
    std::vector<Float_t> fHitZ;
    std::vector<Float_t> fHitSignal;
    std::vector<Float_t> fHitRMS;
         
    // hits belonging to tracks data
    std::vector<Float_t> fTrkHitX;
    std::vector<Float_t> fTrkHitY;
    std::vector<Float_t> fTrkHitZ;
    std::vector<Float_t> fTrkHitSignal;
    std::vector<Float_t> fTrkHitRMS;
    std::vector<Int_t>   fTrkHitTrkIndex;

    // track data
    std::vector<Float_t> fTrackStartX;
    std::vector<Float_t> fTrackStartY;
    std::vector<Float_t> fTrackStartZ;
    std::vector<Float_t> fTrackStartPX;
    std::vector<Float_t> fTrackStartPY;
    std::vector<Float_t> fTrackStartPZ;

    std::vector<Float_t> fTrackEndX;
    std::vector<Float_t> fTrackEndY;
    std::vector<Float_t> fTrackEndZ;
    std::vector<Float_t> fTrackEndPX;
    std::vector<Float_t> fTrackEndPY;
    std::vector<Float_t> fTrackEndPZ;

    std::vector<Float_t> fTrackAvgIon;

    // vertex branches
    std::vector<Float_t> fVertexX;
    std::vector<Float_t> fVertexY;
    std::vector<Float_t> fVertexZ;
    std::vector<Int_t>   fVertexN;
    std::vector<ULong64_t> fVertexT;

    // vertex-track Assn branches
    std::vector<Int_t>   fVTAssn_Vertex;
    std::vector<Float_t> fVTAssn_TrackPx;
    std::vector<Float_t> fVTAssn_TrackPy;
    std::vector<Float_t> fVTAssn_TrackPz;
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

  fWriteMCinfo    = p.get<bool>("WriteMCinfo",true);
  fWriteHits      = p.get<bool>("WriteHits",true);
  fWriteCohInfo   = p.get<bool>("WriteCohInfo",true);
  fWriteHitsInTracks   = p.get<bool>("WriteHitsInTracks",false);

  consumes<std::vector<simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<simb::MCTruth> >(fGeneratorLabel);
  consumes<std::vector<simb::GTruth> >(fGeneratorLabel);
  //consumes<art::Assns<simb::MCTruth, simb::MCParticle> >(fGeantLabel);
  consumes<std::vector<gar::rec::Hit> >(fHitLabel);
  consumes<std::vector<gar::rec::Track> >(fTrackLabel);
  consumes<std::vector<gar::rec::Vertex> >(fVertexLabel);
  consumes<art::Assns<gar::rec::Track, gar::rec::Vertex> >(fVertexLabel);
  consumes<std::vector<gar::rec::TrackIoniz>>(fTrackLabel);
  consumes<art::Assns<rec::TrackIoniz, rec::Track>>(fTrackLabel);
}

void gar::anatree::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("GArAnaTree","GArAnaTree");

  fTree->Branch("Run",         &fRun,            "Run/I");
  fTree->Branch("SubRun",      &fSubRun,         "SubRun/I");
  fTree->Branch("Event",       &fEvent,          "Event/I");

  if (fWriteMCinfo)
    {
      fTree->Branch("NType",       &fNeutrinoType);
      fTree->Branch("CCNC",        &fCCNC);
      fTree->Branch("Mode",        &fMode);
      fTree->Branch("InterT",      &fInteractionType);
      fTree->Branch("MC_Q2",       &fQ2);
      fTree->Branch("MC_W",        &fW);
      fTree->Branch("MC_X",        &fX);
      fTree->Branch("MC_Y",        &fY);
      fTree->Branch("MC_Theta",    &fTheta);
      if (fWriteCohInfo) fTree->Branch("MC_T",        &fT);

      fTree->Branch("Gint",        &fGint);
      fTree->Branch("TgtPDG",      &fTgtPDG);
      fTree->Branch("Weight",      &fWeight);
      fTree->Branch("GT_T",        &fgT);

      fTree->Branch("MCVertX",     &fMCVertexX);
      fTree->Branch("MCVertY",     &fMCVertexY);
      fTree->Branch("MCVertZ",     &fMCVertexZ);
      fTree->Branch("MCNuPx",      &fMCnuPx);
      fTree->Branch("MCNuPy",      &fMCnuPy);
      fTree->Branch("MCNuPz",      &fMCnuPz);

      fTree->Branch("PDG",         &fMCPPDGID);
      fTree->Branch("PDGMom",      &fMCPDGMom);
      fTree->Branch("MCPStartX",   &fMCPStartX);
      fTree->Branch("MCPStartY",   &fMCPStartY);
      fTree->Branch("MCPStartZ",   &fMCPStartZ);
      fTree->Branch("MCPPX",       &fMCPPX);
      fTree->Branch("MCPPY",       &fMCPPY);
      fTree->Branch("MCPPZ",       &fMCPPZ);

      if(fWriteHitsInTracks)        // Write hits in MC tracks
        {
          fTree->Branch("TrajHitX", &fTrajHitX);
          fTree->Branch("TrajHitY", &fTrajHitY);
          fTree->Branch("TrajHitZ", &fTrajHitZ);
          fTree->Branch("TrajHitTrajIndex", &fTrajHitTrajIndex);
        }
    }

  if (fWriteHits)
    {
      fTree->Branch("HitX",        &fHitX);
      fTree->Branch("HitY",        &fHitY);
      fTree->Branch("HitZ",        &fHitZ);
      fTree->Branch("HitSig",      &fHitSignal);
      fTree->Branch("HitRMS",      &fHitRMS);
    }

  if(fWriteHitsInTracks)       // Write hits in reco tracks
    {
      fTree->Branch("TrkHitX", &fTrkHitX);
      fTree->Branch("TrkHitY", &fTrkHitY);
      fTree->Branch("TrkHitZ", &fTrkHitZ);
      fTree->Branch("TrkHitTrkIndex", &fTrkHitTrkIndex);
    }

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

  fTree->Branch("TrackAvgIon",     &fTrackAvgIon);

  fTree->Branch("VertX",           &fVertexX);
  fTree->Branch("VertY",           &fVertexY);
  fTree->Branch("VertZ",           &fVertexZ);
  fTree->Branch("VertN",           &fVertexN);
  fTree->Branch("VertT",           &fVertexT);

  fTree->Branch("VT_Vertex",       &fVTAssn_Vertex);
  fTree->Branch("VT_TrackPx",      &fVTAssn_TrackPx);
  fTree->Branch("VT_TrackPy",      &fVTAssn_TrackPy);
  fTree->Branch("VT_TrackPz",      &fVTAssn_TrackPz);
}


void gar::anatree::analyze(art::Event const & e)
{

  // clear out all our vectors
  if (fWriteMCinfo)
    {
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
      fMCPPDGID.clear();
      fMCPDGMom.clear(); 
      fMCPStartX.clear();
      fMCPStartY.clear();
      fMCPStartZ.clear();
      fMCPPX.clear();
      fMCPPY.clear();
      fMCPPZ.clear();
      if(fWriteHitsInTracks)
        {
          fTrajHitX.clear();
          fTrajHitY.clear();
          fTrajHitZ.clear();
          fTrajHitTrajIndex.clear();
        }
    }
  if (fWriteHits)
    {
      fHitX.clear();
      fHitY.clear();
      fHitZ.clear();
      fHitSignal.clear();
      fHitRMS.clear();
    }
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
  fTrackAvgIon.clear();
  if(fWriteHitsInTracks)
    {
      fTrkHitX.clear();
      fTrkHitY.clear();
      fTrkHitZ.clear();
      fTrkHitTrkIndex.clear();
    }
  fVertexX.clear();
  fVertexY.clear();
  fVertexZ.clear();
  fVertexN.clear();
  fVertexT.clear();
  fVTAssn_Vertex.clear();
  fVTAssn_TrackPx.clear();
  fVTAssn_TrackPy.clear();
  fVTAssn_TrackPz.clear();



  fRun    = e.run();
  fSubRun = e.subRun();
  fEvent  = e.id().event();

  // Get a grip!
  art::Handle< std::vector<simb::MCTruth> > MCTHandle;
  art::Handle< std::vector<simb::GTruth> > GTHandle;
  art::Handle< std::vector<simb::MCParticle> > MCPHandle;
  if (fWriteMCinfo)
    {
    if (!e.getByLabel(fGeneratorLabel, MCTHandle)) 
      {
        throw cet::exception("anatree") 
          << " No simb::MCTruth branch - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    if (!e.getByLabel(fGeneratorLabel, GTHandle)) 
      {
        throw cet::exception("anatree") 
          << " No simb::GTruth branch - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    if (!e.getByLabel(fGeantLabel, MCPHandle)) 
      {
        throw cet::exception("anatree") 
          << " No simb::MCParticle branch - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }

  art::Handle< std::vector<gar::rec::Hit> > HitHandle;
  if (fWriteHits)
    {
    if (!e.getByLabel(fHitLabel, HitHandle)) 
      {
        throw cet::exception("anatree") 
          << " No gar::rec::Hit branch - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }
    }

  art::Handle< std::vector<gar::rec::Track> > TrackHandle;
  if (!e.getByLabel(fTrackLabel, TrackHandle)) 
    {
      throw cet::exception("anatree") 
        << " No gar::rec::Track branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<gar::rec::TrackIoniz> > TrackIonHandle;
  if (!e.getByLabel(fTrackLabel, TrackIonHandle)) 
    {
      throw cet::exception("anatree") 
        << " No gar::rec::TrackIoniz branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }

  art::Handle< std::vector<gar::rec::Vertex> > VertexHandle;
  if (!e.getByLabel(fVertexLabel, VertexHandle)) 
    {
      throw cet::exception("anatree") 
        << " No gar::rec::Vertex branch - "
        << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
    }


  // Pull on the handle now!

  if (fWriteMCinfo)
    {
    // save MCTruth info
    for ( auto const& mct : (*MCTHandle) )
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
        if (fWriteCohInfo)
          {
            double getT = gar::computeT(mct);
            fT.push_back( static_cast<Float_t>(getT) );
          }
        fMCVertexX.push_back(nuw.Nu().EndX());
        fMCVertexY.push_back(nuw.Nu().EndY());
        fMCVertexZ.push_back(nuw.Nu().EndZ());
        fMCnuPx.push_back(nuw.Nu().Px());
        fMCnuPy.push_back(nuw.Nu().Py());
        fMCnuPz.push_back(nuw.Nu().Pz());
      }
    if (fNeutrinoType.size() != 1)
      {
        throw cet::exception("anatree")
          << " MCTruth size != 1 - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    // save GTruth info
    for ( auto const& gt : (*GTHandle) )
      {
        fGint.push_back(gt.fGint);
        fTgtPDG.push_back(gt.ftgtPDG);
        fWeight.push_back(gt.fweight);
        fgT.push_back(gt.fgT);
      }
    if (fGint.size() != 1)
      {
        throw cet::exception("anatree")
          << " GTruth size != 1 - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

    // save MCParticle info
    Int_t mcpIndex = 0;
    for ( auto const& mcp : (*MCPHandle) )
      {
        fMCPPDGID.push_back(mcp.PdgCode());
        int momNumber = mcp.Mother();    int momPDG = 0;
        // shorter loop to find the mother.  MCParticle.fmother appears to be the TrackID of
        // of the mother particle, minus 1.  StatusCode==1 at this point, no point checking that.
        if (momNumber>0)
          {
            for ( auto const& mcp_short : (*MCPHandle) )
              {
                if (mcp_short.TrackId()-1 == momNumber)
                  {
                    momPDG = mcp_short.PdgCode();
                    break;
                  }
              }
          }
        fMCPDGMom.push_back(momPDG);
        const TLorentzVector& position = mcp.Position(0);
        const TLorentzVector& momentum = mcp.Momentum(0);
        fMCPStartX.push_back(position.X());
        fMCPStartY.push_back(position.Y());
        fMCPStartZ.push_back(position.Z());
        fMCPPX.push_back(momentum.Px());
        fMCPPY.push_back(momentum.Py());
        fMCPPZ.push_back(momentum.Pz());

        if(fWriteHitsInTracks)
          {
            for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++)
              {
                fTrajHitX.push_back(mcp.Trajectory().X(iTraj));
                fTrajHitY.push_back(mcp.Trajectory().Y(iTraj));
                fTrajHitZ.push_back(mcp.Trajectory().Z(iTraj));
                fTrajHitTrajIndex.push_back(mcpIndex);
              }
            mcpIndex++;
          }
      }
    }

  // save Hit info
  if (fWriteHits)
    {
    for ( auto const& hit : (*HitHandle) )
      {
        fHitX.push_back(hit.Position()[0]);
        fHitY.push_back(hit.Position()[1]);
        fHitZ.push_back(hit.Position()[2]);
        fHitSignal.push_back(hit.Signal());
        fHitRMS.push_back(hit.RMS());
      }
    }

  // save Track info
  const art::FindManyP<gar::rec::Hit>       findManyHits(TrackHandle,e,fTrackLabel);
  const art::FindManyP<gar::rec::TrackIoniz> findIonizations(TrackHandle,e,fTrackLabel);
  size_t iTrack = 0;
  for ( auto const& track : (*TrackHandle) )
    { // track is a gar::rec::Track, not a gar::rec::TrackPar
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

      if (fWriteHitsInTracks)
        {
          int nTrackedHits = 0;
          if (findManyHits.isValid())
            {   // hits is a vector of gar::rec::Hit
              auto const& hits = findManyHits.at(iTrack);
              nTrackedHits = hits.size();
            }
          for (int iTrackedHit=0; iTrackedHit<nTrackedHits; iTrackedHit++)
            {
              auto const& hit = *(findManyHits.at(iTrack).at(iTrackedHit));
              fTrkHitX.push_back(hit.Position()[0]);
              fTrkHitY.push_back(hit.Position()[1]);
              fTrkHitZ.push_back(hit.Position()[2]);
              fTrkHitTrkIndex.push_back((Int_t)iTrack);
            }
        }

      if (findIonizations.isValid())
        {
          // No calibration for now
          rec::TrackIoniz ionization = *(findIonizations.at(iTrack).at(0));
          gar::processIonizationInfo(ionization);
          fTrackAvgIon.push_back( gar::AverageIonization(ionization) );
        }
      iTrack++;
    }

  // save Vertex and Track-Vertex association info
  const art::FindManyP<gar::rec::Track> findManyTrack(VertexHandle,e,fVertexLabel);
  size_t iVertex = 0;
  for ( auto const& vertex : (*VertexHandle) )
    {
      fVertexX.push_back(vertex.Position()[0]);
      fVertexY.push_back(vertex.Position()[1]);
      fVertexZ.push_back(vertex.Position()[2]);
      fVertexT.push_back(vertex.Time());

      int nVertexedTracks = 0;
      if ( findManyTrack.isValid() )
        { // tracks is a vector of gar::rec::Track
          auto const& tracks = findManyTrack.at(iVertex);
          nVertexedTracks = tracks.size();
        }
      fVertexN.push_back(nVertexedTracks);

      for (int iVertexedTrack=0; iVertexedTrack<nVertexedTracks; ++iVertexedTrack)
        {
          fVTAssn_Vertex.push_back(iVertex);  // 1 push of iVertex for each track in vertex
          // Get this vertexed gar::rec::Track
          auto const& track = *(findManyTrack.at(iVertex).at(iVertexedTrack));

          // Determine which end of track is in the vertex.  Someday, upgrade the Assns
          // created in vertexfinder1_module.cc to include this info
          float beginsAt[3], endsAt[3];
          for (int i=0; i<3; ++i)
            {
              beginsAt[i] = track.Vertex()[i] - vertex.Position()[i];
              endsAt[i]   = track.End()[i]    - vertex.Position()[i];
            }
          float d2Begin = (beginsAt[0]*beginsAt[0] +beginsAt[1]*beginsAt[1] +beginsAt[2]*beginsAt[2]);
          float d2End   = (endsAt[0]*endsAt[0]     +endsAt[1]*endsAt[1]     +endsAt[2]*endsAt[2]);
          bool usebeg = d2Begin < d2End;
          if (usebeg)
            {
              fVTAssn_TrackPx.push_back( track.Momentum_beg()*track.VtxDir()[0] );
              fVTAssn_TrackPy.push_back( track.Momentum_beg()*track.VtxDir()[1] );
              fVTAssn_TrackPz.push_back( track.Momentum_beg()*track.VtxDir()[2] );
            }
          else
            {
              fVTAssn_TrackPx.push_back( track.Momentum_end()*track.EndDir()[0] );
              fVTAssn_TrackPy.push_back( track.Momentum_end()*track.EndDir()[1] );
              fVTAssn_TrackPz.push_back( track.Momentum_end()*track.EndDir()[2] );
            }
        }
      ++iVertex;
     }

  fTree->Fill();
}



DEFINE_ART_MODULE(gar::anatree)
