#include "MCP_Skimmer.h"

#include "TDirectory.h"
#include "TVector3.h"

#include <iostream>
#include <algorithm>
#include <functional>
#include <map>

MCP_Skimmer::MCP_Skimmer()
: _intfile(nullptr), _skimfile(nullptr), _skimtree(nullptr), _inttree(nullptr), _debug(false), _util(new Utils()), _infile(""), _outfile("")
{

}

MCP_Skimmer::MCP_Skimmer(std::string infilename, std::string outfilename)
: _skimfile(nullptr), _skimtree(nullptr), _debug(false), _util(new Utils()), _infile(infilename), _outfile(outfilename)
{

}

MCP_Skimmer::~MCP_Skimmer()
{

}

void MCP_Skimmer::ClearVectors()
{
    _MC_Q2.clear();
    _MC_W.clear();
    _MC_Y.clear();
    _MC_X.clear();
    _MC_Theta.clear();
    _MC_T.clear();
    _MCVertX.clear();
    _MCVertY.clear();
    _MCVertZ.clear();
    _MCNuPx.clear();
    _MCNuPy.clear();
    _MCNuPz.clear();
    _NType.clear();
    _CCNC.clear();
    _Mode.clear();
    _Gint.clear();
    _TgtPDG.clear();
    _GT_T.clear();
    _InterT.clear();
    _Weight.clear();

    _PDG.clear();
    _Mother.clear();
    _PDGMother.clear();
    _MCPTrkID.clear();
    _MCPTime.clear();
    _MCPStartX.clear();
    _MCPStartY.clear();
    _MCPStartZ.clear();
    _MCPStartPX.clear();
    _MCPStartPY.clear();
    _MCPStartPZ.clear();
    _MCPProc.clear();
    _MCPEndProc.clear();
    _MCPEndX.clear();
    _MCPEndY.clear();
    _MCPEndZ.clear();
    _TrajMCPX.clear();
    _TrajMCPY.clear();
    _TrajMCPZ.clear();
    _TrajMCPTrajIndex.clear();
}

void MCP_Skimmer::SkimMCParticle()
{
    int           	 Event = 0;
    int           	 SubRun = 0;
    int           	 Run = 0;

    std::vector<float> *MC_Q2 = 0;
    std::vector<float> *MC_W = 0;
    std::vector<float> *MC_Y = 0;
    std::vector<float> *MC_X = 0;
    std::vector<float> *MC_Theta = 0;
    std::vector<float> *MC_T = 0;
    std::vector<float>   *MCVertX = 0;
    std::vector<float>   *MCVertY = 0;
    std::vector<float>   *MCVertZ = 0;
    std::vector<float>   *MCNuPx = 0;
    std::vector<float>   *MCNuPy = 0;
    std::vector<float>   *MCNuPz = 0;
    std::vector<int>     *NType = 0;
    std::vector<int>     *CCNC = 0;
    std::vector<int>     *Mode = 0;
    std::vector<int>     *Gint=0;
    std::vector<int>     *TgtPDG=0;
    std::vector<int>     *GT_T=0;
    std::vector<int>     *InterT=0;
    std::vector<float>   *Weight=0;

    std::vector<int>     *PDG = 0;
    std::vector<int>     *Mother=0;
    std::vector<int>     *PDGMother=0;
    std::vector<int>     *MCPTrkID = 0;
    std::vector<float>   *MCPTime = 0;
    std::vector<float>   *MCPStartX = 0;
    std::vector<float>   *MCPStartY = 0;
    std::vector<float>   *MCPStartZ = 0;
    std::vector<float>   *MCPStartPX = 0;
    std::vector<float>   *MCPStartPY = 0;
    std::vector<float>   *MCPStartPZ = 0;
    std::vector<std::string>   *MCPProc = 0;
    std::vector<std::string>   *MCPEndProc = 0;
    std::vector<float>   *MCPEndX = 0;
    std::vector<float>   *MCPEndY = 0;
    std::vector<float>   *MCPEndZ = 0;
    std::vector<float>   *TrajMCPX = 0;
    std::vector<float>   *TrajMCPY = 0;
    std::vector<float>   *TrajMCPZ = 0;
    std::vector<int>     *TrajMCPTrajIndex = 0;

    _inttree->SetBranchStatus("*", 0);
    _inttree->SetBranchStatus("Event", 1);
    _inttree->SetBranchStatus("SubRun", 1);
    _inttree->SetBranchStatus("Run", 1);

    _inttree->SetBranchStatus("NType", 1);
    _inttree->SetBranchStatus("CCNC", 1);
    _inttree->SetBranchStatus("MC_Q2", 1);
    _inttree->SetBranchStatus("MC_W", 1);
    _inttree->SetBranchStatus("MC_Y", 1);
    _inttree->SetBranchStatus("MC_X", 1);
    _inttree->SetBranchStatus("MC_Theta", 1);
    _inttree->SetBranchStatus("MC_T", 1);
    _inttree->SetBranchStatus("Mode", 1);
    _inttree->SetBranchStatus("Gint", 1);
    _inttree->SetBranchStatus("TgtPDG", 1);
    _inttree->SetBranchStatus("GT_T", 1);
    _inttree->SetBranchStatus("MCVertX", 1);
    _inttree->SetBranchStatus("MCVertY", 1);
    _inttree->SetBranchStatus("MCVertZ", 1);
    _inttree->SetBranchStatus("MCNuPx", 1);
    _inttree->SetBranchStatus("MCNuPy", 1);
    _inttree->SetBranchStatus("MCNuPz", 1);
    _inttree->SetBranchStatus("InterT", 1);
    _inttree->SetBranchStatus("Weight", 1);

    //MC info
    _inttree->SetBranchStatus("PDG", 1);
    _inttree->SetBranchStatus("MCPTrkID", 1);
    _inttree->SetBranchStatus("MCPTime", 1);
    _inttree->SetBranchStatus("MCPStartX", 1);
    _inttree->SetBranchStatus("MCPStartY", 1);
    _inttree->SetBranchStatus("MCPStartZ", 1);
    _inttree->SetBranchStatus("MCPEndX", 1);
    _inttree->SetBranchStatus("MCPEndY", 1);
    _inttree->SetBranchStatus("MCPEndZ", 1);
    _inttree->SetBranchStatus("Mother", 1);
    _inttree->SetBranchStatus("PDGMother", 1);
    _inttree->SetBranchStatus("MCPStartPX", 1);
    _inttree->SetBranchStatus("MCPStartPY", 1);
    _inttree->SetBranchStatus("MCPStartPZ", 1);
    _inttree->SetBranchStatus("MCPProc", 1);
    _inttree->SetBranchStatus("MCPEndProc", 1);
    _inttree->SetBranchStatus("TrajMCPX", 1);
    _inttree->SetBranchStatus("TrajMCPY", 1);
    _inttree->SetBranchStatus("TrajMCPZ", 1);
    _inttree->SetBranchStatus("TrajMCPTrajIndex", 1);

    //-------------------

    _inttree->SetBranchAddress("Event", &Event);
    _inttree->SetBranchAddress("SubRun", &SubRun);
    _inttree->SetBranchAddress("Run", &Run);

    //Generator info
    _inttree->SetBranchAddress("NType", &NType);
    _inttree->SetBranchAddress("CCNC", &CCNC);
    _inttree->SetBranchAddress("MC_Q2", &MC_Q2);
    _inttree->SetBranchAddress("MC_W", &MC_W);
    _inttree->SetBranchAddress("MC_Y", &MC_Y);
    _inttree->SetBranchAddress("MC_X", &MC_X);
    _inttree->SetBranchAddress("MC_Theta", &MC_Theta);
    _inttree->SetBranchAddress("MC_T", &MC_T);
    _inttree->SetBranchAddress("Mode", &Mode);
    _inttree->SetBranchAddress("Gint", &Gint);
    _inttree->SetBranchAddress("TgtPDG", &TgtPDG);
    _inttree->SetBranchAddress("GT_T", &GT_T);
    _inttree->SetBranchAddress("MCVertX", &MCVertX);
    _inttree->SetBranchAddress("MCVertY", &MCVertY);
    _inttree->SetBranchAddress("MCVertZ", &MCVertZ);
    _inttree->SetBranchAddress("MCNuPx", &MCNuPx);
    _inttree->SetBranchAddress("MCNuPy", &MCNuPy);
    _inttree->SetBranchAddress("MCNuPz", &MCNuPz);
    _inttree->SetBranchAddress("InterT", &InterT);
    _inttree->SetBranchAddress("Weight", &Weight);

    //MC info
    _inttree->SetBranchAddress("PDG", &PDG);
    _inttree->SetBranchAddress("MCPTrkID", &MCPTrkID);
    _inttree->SetBranchAddress("MCPTime", &MCPTime);
    _inttree->SetBranchAddress("MCPStartX", &MCPStartX);
    _inttree->SetBranchAddress("MCPStartY", &MCPStartY);
    _inttree->SetBranchAddress("MCPStartZ", &MCPStartZ);
    _inttree->SetBranchAddress("MCPEndX", &MCPEndX);
    _inttree->SetBranchAddress("MCPEndY", &MCPEndY);
    _inttree->SetBranchAddress("MCPEndZ", &MCPEndZ);
    _inttree->SetBranchAddress("Mother", &Mother);
    _inttree->SetBranchAddress("PDGMother", &PDGMother);
    _inttree->SetBranchAddress("MCPStartPX", &MCPStartPX);
    _inttree->SetBranchAddress("MCPStartPY", &MCPStartPY);
    _inttree->SetBranchAddress("MCPStartPZ", &MCPStartPZ);
    _inttree->SetBranchAddress("MCPProc", &MCPProc);
    _inttree->SetBranchAddress("MCPEndProc", &MCPEndProc);
    _inttree->SetBranchAddress("TrajMCPX", &TrajMCPX);
    _inttree->SetBranchAddress("TrajMCPY", &TrajMCPY);
    _inttree->SetBranchAddress("TrajMCPZ", &TrajMCPZ);
    _inttree->SetBranchAddress("TrajMCPTrajIndex", &TrajMCPTrajIndex);

    // Main event loop
    for( int entry = 0; entry < _inttree->GetEntries(); entry++ )
    {
        _inttree->GetEntry(entry);

        this->ClearVectors();

        if(Event%100 == 0)
        std::cout << "Treating Evt " << Event << std::endl;

        _Event = Event;
        _Run = Run;
        _SubRun = SubRun;

        for(size_t i = 0; i < NType->size(); i++)
        {
            _CCNC.push_back(CCNC->at(i));
            _NType.push_back(NType->at(i));
            _MC_Q2.push_back(MC_Q2->at(i));
            _MC_W.push_back(MC_W->at(i));
            _MC_Y.push_back(MC_Y->at(i));
            _MC_X.push_back(MC_X->at(i));
            _MC_Theta.push_back(MC_Theta->at(i));
            _Mode.push_back(Mode->at(i));
            _MC_T.push_back(MC_T->at(i));
            _InterT.push_back(InterT->at(i));
            _MCVertX.push_back(MCVertX->at(i));
            _MCVertY.push_back(MCVertY->at(i));
            _MCVertZ.push_back(MCVertZ->at(i));
            _MCNuPx.push_back(MCNuPx->at(i));
            _MCNuPy.push_back(MCNuPy->at(i));
            _MCNuPz.push_back(MCNuPz->at(i));
        }

        for(size_t i = 0; i < Gint->size(); i++)
        {
            _Gint.push_back(Gint->at(i));
            _TgtPDG.push_back(TgtPDG->at(i));
            _GT_T.push_back(GT_T->at(i));
            _Weight.push_back(Weight->at(i));
        }

        std::map<int, int> mcp_kid_mother_map;

        for(size_t i = 0; i < MCPStartPX->size(); i++ )
        {
            //check if daughter has mother until we find the original mother
            //get current kid trk id
            int this_kid_trkid = MCPTrkID->at(i);
            int this_kid_mother = Mother->at(i);
            //Mother->at(i) returns the mother trkid

            //Has no mother -> primaries
            if(this_kid_mother == 0) {
                // std::cout << "Kid " << this_kid_trkid << " Mother " << this_kid_mother << std::endl;
                continue;
            }

            //loop over the mcp again and go back in history need to check Mother->at(j) until it is 0 then we found the original mcp
            for(size_t j = 0; j < MCPStartPX->size(); j++ )
            {
                int mother_trk_id = MCPTrkID->at(j);
                int mother = Mother->at(j);

                //found the mother
                if(mother_trk_id == this_kid_mother)
                {
                    //need to check if it has mother if yes loop until mother is 0
                    if(mother == 0)
                    {
                        //found the original particle
                        // std::cout << "Kid " << this_kid_trkid << " Mother " << this_kid_mother << std::endl;
                        //fill a map linking kid track id to mother track id to know which ones to skip and replace their mother by the original one
                        mcp_kid_mother_map[this_kid_trkid] = this_kid_mother;
                        break;
                    }
                    else
                    {
                        this_kid_mother = mother;
                        j = 0;
                    }
                }
            }
        }

        //Vector of indexes to keep
        std::vector<int> IndexesToKeep;

        //Check which indexes to keep
        for(size_t i = 0; i < MCPStartPX->size(); i++)
        {
            TVector3 spoint(MCPStartX->at(i), MCPStartY->at(i), MCPStartZ->at(i));//start point
            TVector3 epoint(MCPEndX->at(i), MCPEndY->at(i), MCPEndZ->at(i));//end point
            std::string process = MCPProc->at(i);
            std::string endprocess = MCPEndProc->at(i);

            //found a particle that has already an ancestor
            if( mcp_kid_mother_map.find(MCPTrkID->at(i)) !=  mcp_kid_mother_map.end() )
            {
                // However keep the daughters for specific ones
                // - D0s and V0s (kinks and decays) such as gamma, pi0, K0, muon, pion and kaon
                // - backscatters from calo to tracker important for background
                // - breamstrahlung photons / delta rays

                //keep D0s and V0s from decays, conversions (a lot of compton or Ioni)
                if( std::find(daughtersToKeep.begin(), daughtersToKeep.end(), PDGMother->at(i)) != daughtersToKeep.end() && _util->hasOriginInTracker(spoint) && (process == "Decay" || process == "conv") ) {
                    if(_debug){
                        std::cout << "Keeping D0 or V0" << std::endl;
                        std::cout << "Index " << i << " TrkID " << MCPTrkID->at(i) << " pdg " << PDG->at(i) << " mother pdg " << PDGMother->at(i);
                        std::cout << " process " << process << " endprocess " << endprocess << std::endl;
                        std::cout << " Start point X: " << spoint.X() << " Y: " << spoint.Y() << " Z: " << spoint.Z() << std::endl;
                        std::cout << " End point X: " << epoint.X() << " Y: " << epoint.Y() << " Z: " << epoint.Z() << std::endl;
                        std::cout << " Origin in tracker " << _util->hasOriginInTracker(spoint) << std::endl;
                        std::cout << " Decayed in Calo " << _util->hasDecayedInCalo(epoint) << std::endl;
                    }
                    IndexesToKeep.push_back(i);
                }
                //check for backscatter
                else if( _util->isBackscatter(spoint, epoint) ) {
                    if(_debug){
                        std::cout << "Keeping Backscatter" << std::endl;
                        std::cout << "Index " << i << " TrkID " << MCPTrkID->at(i) << " pdg " << PDG->at(i) << " mother pdg " << PDGMother->at(i);
                        std::cout << " process " << process << " endprocess " << endprocess << std::endl;
                        std::cout << " Start point X: " << spoint.X() << " Y: " << spoint.Y() << " Z: " << spoint.Z() << std::endl;
                        std::cout << " End point X: " << epoint.X() << " Y: " << epoint.Y() << " Z: " << epoint.Z() << std::endl;
                        std::cout << " Origin in tracker " << _util->hasOriginInTracker(spoint) << std::endl;
                        std::cout << " Decayed in Calo " << _util->hasDecayedInCalo(epoint) << std::endl;
                    }
                    IndexesToKeep.push_back(i);
                }
                //if Bremsstrahlung photon
                else if ( _util->isBremsstrahlung(spoint, PDG->at(i), PDGMother->at(i)) && process == "eBrem" ) {
                    if(_debug){
                        std::cout << "Keeping Bremsstrahlung" << std::endl;
                        std::cout << "Index " << i << " TrkID " << MCPTrkID->at(i) << " pdg " << PDG->at(i) << " mother pdg " << PDGMother->at(i);
                        std::cout << " process " << process << " endprocess " << endprocess << std::endl;
                        std::cout << " Start point X: " << spoint.X() << " Y: " << spoint.Y() << " Z: " << spoint.Z() << std::endl;
                        std::cout << " End point X: " << epoint.X() << " Y: " << epoint.Y() << " Z: " << epoint.Z() << std::endl;
                        std::cout << " Origin in tracker " << _util->hasOriginInTracker(spoint) << std::endl;
                        std::cout << " Decayed in Calo " << _util->hasDecayedInCalo(epoint) << std::endl;
                    }
                    IndexesToKeep.push_back(i);
                }
                else {
                    continue;
                }
            }

            //keep only primaries whatever they are from (inside TPC or outside)
            if(MCPProc->at(i) == "primary"){
                if(_debug){
                    std::cout << "Keeping Primary particle" << std::endl;
                    std::cout << "Index " << i << " TrkID " << MCPTrkID->at(i) << " pdg " << PDG->at(i) << " mother pdg " << PDGMother->at(i);
                    std::cout << " process " << process << " endprocess " << endprocess << std::endl;
                    std::cout << " Start point X: " << spoint.X() << " Y: " << spoint.Y() << " Z: " << spoint.Z() << std::endl;
                    std::cout << " End point X: " << epoint.X() << " Y: " << epoint.Y() << " Z: " << epoint.Z() << std::endl;
                    std::cout << " Origin in tracker " << _util->hasOriginInTracker(spoint) << std::endl;
                    std::cout << " Decayed in Calo " << _util->hasDecayedInCalo(epoint) << std::endl;
                }
                IndexesToKeep.push_back(i);
            }
        }


        for(size_t i = 0; i < IndexesToKeep.size(); i++)
        {
            _PDG.push_back(PDG->at(i));
            _MCPStartPX.push_back(MCPStartPX->at(i));
            _MCPStartPY.push_back(MCPStartPY->at(i));
            _MCPStartPZ.push_back(MCPStartPZ->at(i));
            _MCPTrkID.push_back(MCPTrkID->at(i));
            _MCPTime.push_back(MCPTime->at(i));
            _MCPStartX.push_back(MCPStartX->at(i));
            _MCPStartY.push_back(MCPStartY->at(i));
            _MCPStartZ.push_back(MCPStartZ->at(i));
            _MCPEndX.push_back(MCPEndX->at(i));
            _MCPEndY.push_back(MCPEndY->at(i));
            _MCPEndZ.push_back(MCPEndZ->at(i));
            _MCPProc.push_back(MCPProc->at(i));
            _MCPEndProc.push_back(MCPEndProc->at(i));
            _Mother.push_back(Mother->at(i));
            _PDGMother.push_back(PDGMother->at(i));

            for(size_t itraj = 0; itraj < TrajMCPX->size(); itraj++){
                //keep only traj from kept mcp
	      if(TrajMCPTrajIndex->at(itraj) == (int) i){
                    _TrajMCPX.push_back(TrajMCPX->at(itraj));
                    _TrajMCPY.push_back(TrajMCPY->at(itraj));
                    _TrajMCPZ.push_back(TrajMCPZ->at(itraj));
                    _TrajMCPTrajIndex.push_back(TrajMCPTrajIndex->at(itraj));
                }
            }
        }

        this->FillTTree();
    }//end tree loop
}

bool MCP_Skimmer::BookTFile()
{
    _intfile = new TFile( _infile.c_str() );
    if(nullptr == _intfile)
    return false;

    _inttree = (TTree*) _intfile->Get( "anatree/GArAnaTree" );
    if(nullptr == _inttree)
    return false;

    _skimfile = new TFile(_outfile.c_str(), "RECREATE");
    if(nullptr == _skimfile){
        return false;
    }
    else
    {

        TDirectory *anatree = _skimfile->mkdir("anatree");
        anatree->cd();
        _skimtree = new TTree( "GArAnaTree", "GArAnaTree" );

        _skimtree->Branch("Event", &_Event);
        _skimtree->Branch("SubRun", &_SubRun);
        _skimtree->Branch("Run", &_Run);

        //Generator info
        _skimtree->Branch("NType", &_NType);
        _skimtree->Branch("CCNC", &_CCNC);
        _skimtree->Branch("MC_Q2", &_MC_Q2);
        _skimtree->Branch("MC_W", &_MC_W);
        _skimtree->Branch("MC_Y", &_MC_Y);
        _skimtree->Branch("MC_X", &_MC_X);
        _skimtree->Branch("MC_Theta", &_MC_Theta);
        _skimtree->Branch("MC_T", &_MC_T);
        _skimtree->Branch("Mode", &_Mode);
        _skimtree->Branch("Gint", &_Gint);
        _skimtree->Branch("TgtPDG", &_TgtPDG);
        _skimtree->Branch("GT_T", &_GT_T);
        _skimtree->Branch("MCVertX", &_MCVertX);
        _skimtree->Branch("MCVertY", &_MCVertY);
        _skimtree->Branch("MCVertZ", &_MCVertZ);
        _skimtree->Branch("MCNuPx", &_MCNuPx);
        _skimtree->Branch("MCNuPy", &_MCNuPy);
        _skimtree->Branch("MCNuPz", &_MCNuPz);
        _skimtree->Branch("InterT", &_InterT);
        _skimtree->Branch("Weight", &_Weight);

        //MC info
        _skimtree->Branch("PDG", &_PDG);
        _skimtree->Branch("MCPTrkID", &_MCPTrkID);
        _skimtree->Branch("MCPTime", &_MCPTime);
        _skimtree->Branch("MCPStartX", &_MCPStartX);
        _skimtree->Branch("MCPStartY", &_MCPStartY);
        _skimtree->Branch("MCPStartZ", &_MCPStartZ);
        _skimtree->Branch("MCPEndX", &_MCPEndX);
        _skimtree->Branch("MCPEndY", &_MCPEndY);
        _skimtree->Branch("MCPEndZ", &_MCPEndZ);
        _skimtree->Branch("Mother", &_Mother);
        _skimtree->Branch("PDGMother", &_PDGMother);
        _skimtree->Branch("MCPStartPX", &_MCPStartPX);
        _skimtree->Branch("MCPStartPY", &_MCPStartPY);
        _skimtree->Branch("MCPStartPZ", &_MCPStartPZ);
        _skimtree->Branch("MCPProc", &_MCPProc);
        _skimtree->Branch("MCPEndProc", &_MCPEndProc);
        _skimtree->Branch("TrajMCPX", &_TrajMCPX);
        _skimtree->Branch("TrajMCPY", &_TrajMCPY);
        _skimtree->Branch("TrajMCPZ", &_TrajMCPZ);
        _skimtree->Branch("TrajMCPTrajIndex", &_TrajMCPTrajIndex);

        return true;
    }
}

void MCP_Skimmer::FillTTree()
{
    if(_skimtree != nullptr) {
        _skimfile->cd("anatree");
        _skimtree->Fill();
    }

    return;
}

void MCP_Skimmer::CloseTFile()
{
    if(_skimfile != nullptr) {
        _skimfile->Close();
        return;
    }
    else{ return; }
}

void MCP_Skimmer::WriteTTree()
{
    if(_skimtree != nullptr) {
        _skimfile->cd("anatree");
        _skimtree->Write();
        _skimfile->Close();
    }

    return;
}
