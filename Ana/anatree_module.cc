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

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "nusimdata/SimulationBase/GTruth.h"
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

			// the analysis tree

			TTree *fTree;

			// global event info

			int fEvent;        ///< number of the event being processed
			int fRun;          ///< number of the run being processed
			int fSubRun;       ///< number of the sub-run being processed

			// MCTruth data
			// GENIE kinematics computed ignoring Fermi momentum and the off-shellness of the bound nucleon.
			// Get the right kinematics here, except T has to be computed

			std::vector<int>   fNeutrinoType;
			std::vector<int>   fCCNC;
			std::vector<int>   fMode;
			std::vector<int>   fInteractionType;
			std::vector<float> fQ2;
			std::vector<float> fW;
			std::vector<float> fX;
			std::vector<float> fY;
			std::vector<float> fTheta;
			std::vector<float> fT;

			// GTruth data

			std::vector<int>   fGint;
			std::vector<int>   fTgtPDG;
			std::vector<float> fWeight;
			std::vector<float> fgT;

			// MCParticle data

			std::vector<int>   fMCPPDGID;
			std::vector<int>   fMCPDGMom;
			std::vector<float> fMCPStartX;
			std::vector<float> fMCPStartY;
			std::vector<float> fMCPStartZ;
			std::vector<float> fMCPPX;
			std::vector<float> fMCPPY;
			std::vector<float> fMCPPZ;
			std::vector< std::vector<float> > fTrajHitX;
			std::vector< std::vector<float> > fTrajHitY;
			std::vector< std::vector<float> > fTrajHitZ;


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
			std::vector< std::vector<float> > fTrackHitX;
			std::vector< std::vector<float> > fTrackHitY;
			std::vector< std::vector<float> > fTrackHitZ;

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

	fWriteMCinfo    = p.get<bool>("WriteMCinfo",true);
	fWriteHits      = p.get<bool>("WriteHits",true);
	fWriteCohInfo   = p.get<bool>("WriteCohInfo",true);

	consumes<std::vector<simb::MCParticle> >(fGeantLabel);
	consumes<std::vector<simb::MCTruth> >(fGeneratorLabel);
	consumes<std::vector<simb::GTruth> >(fGeneratorLabel);
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

		fTree->Branch("PDG",       &fMCPPDGID);
		fTree->Branch("PDGMom",    &fMCPDGMom);
		fTree->Branch("MCPStartX", &fMCPStartX);
		fTree->Branch("MCPStartY", &fMCPStartY);
		fTree->Branch("MCPStartZ", &fMCPStartZ);
		fTree->Branch("MCPPX", &fMCPPX);
		fTree->Branch("MCPPY", &fMCPPY);
		fTree->Branch("MCPPZ", &fMCPPZ);
		fTree->Branch("TrajHitX", &fTrajHitX);
		fTree->Branch("TrajHitY", &fTrajHitY);
		fTree->Branch("TrajHitZ", &fTrajHitZ);

	}

	if (fWriteHits)
	{
		fTree->Branch("HitX", &fHitX);
		fTree->Branch("HitY", &fHitY);
		fTree->Branch("HitZ", &fHitZ);
		fTree->Branch("HitSig", &fHitSignal);
		fTree->Branch("HitRMS", &fHitRMS);
	}

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
	fTree->Branch("TrackHitX", &fTrackHitX);
	fTree->Branch("TrackHitY", &fTrackHitY);
	fTree->Branch("TrackHitZ", &fTrackHitZ);

	fTree->Branch("VertX", &fVertexX);
	fTree->Branch("VertY", &fVertexY);
	fTree->Branch("VertZ", &fVertexZ);
	fTree->Branch("VertN", &fVertexNTracks);    
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
		fTrajHitX.clear();
		fTrajHitY.clear();
		fTrajHitZ.clear();

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
	fTrackHitX.clear();
	fTrackHitY.clear();
	fTrackHitZ.clear();
	fVertexX.clear();
	fVertexY.clear();
	fVertexZ.clear();
	fVertexNTracks.clear();

	fEvent  = e.id().event(); 
	fRun    = e.run();
	fSubRun = e.subRun();

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
			fNeutrinoType.push_back(mct.GetNeutrino().Nu().PdgCode());
			fCCNC.push_back(mct.GetNeutrino().CCNC());
			fMode.push_back(mct.GetNeutrino().Mode());
			fInteractionType.push_back(mct.GetNeutrino().InteractionType());
			fQ2.push_back(mct.GetNeutrino().QSqr());
			fW.push_back(mct.GetNeutrino().W());
			fX.push_back(mct.GetNeutrino().X());
			fY.push_back(mct.GetNeutrino().Y());
			fTheta.push_back(mct.GetNeutrino().Theta());
			if (fWriteCohInfo)
			{
				double getT = computeT(mct);
				fT.push_back( static_cast<float>(getT) );
			}
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

		for ( auto const& mcp : (*MCPHandle) )
		{
			fMCPPDGID.push_back(mcp.PdgCode());
			int momNumber = mcp.Mother();	int momPDG = 0;
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
			const TLorentzVector& pos = mcp.Position(0);
			const TLorentzVector& mom = mcp.Momentum(0);
			fMCPStartX.push_back(pos.X());
			fMCPStartY.push_back(pos.Y());
			fMCPStartZ.push_back(pos.Z());
			fMCPPX.push_back(mom.Px());
			fMCPPY.push_back(mom.Py());
			fMCPPZ.push_back(mom.Pz());
			//TODO: save traj hits here
			std::vector<float> tempTrjHitX;
			std::vector<float> tempTrjHitY;
			std::vector<float> tempTrjHitZ;
			for(uint iTraj=0; iTraj < mcp.Trajectory().size(); iTraj++){
				tempTrjHitX.push_back(mcp.Trajectory().X(iTraj));
				tempTrjHitY.push_back(mcp.Trajectory().Y(iTraj));
				tempTrjHitZ.push_back(mcp.Trajectory().Z(iTraj));
			}
			fTrajHitX.push_back(tempTrjHitX);
			fTrajHitY.push_back(tempTrjHitY);
			fTrajHitZ.push_back(tempTrjHitZ);
		}
	}

	if (fWriteHits)
	{
		for ( auto const& hit : (*HitHandle) )
		{
			fHitX.push_back(hit.Position()[0]);
			fHitY.push_back(hit.Position()[1]);
			fHitZ.push_back(hit.Position()[2]);
			fHitSignal.push_back(hit.Signal());
			fHitRMS.push_back(hit.RMS());

		} //for hits
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
		std::vector<float> tempTrkHitX;
		std::vector<float> tempTrkHitY;
		std::vector<float> tempTrkHitZ;
		for(size_t ii=0; (ii<track.NHits() && ii<200); ii++){
			tempTrkHitX.push_back(track.HitXs()[ii]);
			tempTrkHitY.push_back(track.HitYs()[ii]);
			tempTrkHitZ.push_back(track.HitZs()[ii]);
		}
		fTrackHitX.push_back(tempTrkHitX);
		fTrackHitY.push_back(tempTrkHitY);
		fTrackHitZ.push_back(tempTrkHitZ);
	}

	for ( auto const& vertex : (*VertexHandle) )
	{
		fVertexX.push_back(vertex.Position()[0]);
		fVertexY.push_back(vertex.Position()[1]);
		fVertexZ.push_back(vertex.Position()[2]);

		// count up the tracks belonging to this vertex
		int ntracks = 0;
		// Leo: a bug, I guess:
		// const art::FindManyP<gar::rec::Track> findManyTrack(VertexHandle,e,fTrackLabel);
		const art::FindManyP<gar::rec::Track> findManyTrack(VertexHandle,e,fVertexLabel);
		if ( findManyTrack.isValid() )
		{
			ntracks = findManyTrack.size();
		}
		fVertexNTracks.push_back(ntracks);
	}

	fTree->Fill();
}


// Coherent pion analysis specific code

double gar::anatree::computeT( simb::MCTruth theMCTruth )
{
	int nPart = theMCTruth.NParticles();
	enum { nu, mu, pi};
	double E[3], Px[3], Py[3], Pz[3];
	E[nu] = E[mu] = E[pi] = -1e42;

	for (int i=0; i<3;++i)
	{
		Px[i] = 0; 
		Py[i] = 0;
		Pz[i] = 0;
		E[i] = 0;
	}  

	// Find t from the MCParticles via the
	for (int iPart=0; iPart<nPart; iPart++)
	{
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
	}

	// Compute t; reuse nu 4-vector to get first q, then t.
	E[nu] -= E[mu];   Px[nu] -= Px[mu];   Py[nu] -= Py[mu];   Pz[nu] -= Pz[mu];
	E[nu] -= E[pi];   Px[nu] -= Px[pi];   Py[nu] -= Py[pi];   Pz[nu] -= Pz[pi];
	double t = E[nu]*E[nu] -Px[nu]*Px[nu] -Py[nu]*Py[nu] -Pz[nu]*Pz[nu];
	return t*t;
}



DEFINE_ART_MODULE(gar::anatree)
