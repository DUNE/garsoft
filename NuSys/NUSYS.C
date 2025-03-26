#define CAF_cxx
#ifdef CAF_cxx

#include "NUSYS.h"

CAF::CAF(std::string filename, bool isGas)
{
  cafFile = new TFile(filename.c_str(), "RECREATE");
  cafMVA = new TTree("caf", "caf");
  cafPOT = new TTree("meta", "meta");
  genie = new TTree("genieEvt", "genieEvt");

  // initialize the GENIE record
  mcrec = NULL;
  cafMVA->Branch("ievt", &ievt, "ievt/I");
  cafMVA->Branch("run", &run, "run/I");
  cafMVA->Branch("subrun", &subrun, "subrun/I");
  cafMVA->Branch("event", &event, "event/I");

  genie->Branch("genie_record", &mcrec);

  cafPOT->Branch("pot", &pot, "pot/D");
  cafPOT->Branch("run", &meta_run, "run/I");
  cafPOT->Branch("subrun", &meta_subrun, "subrun/I");
  cafPOT->Branch("version", &version, "version/I");
}

CAF::~CAF() {}

void CAF::fill()
{
  cafMVA->Fill();
  genie->Fill();
}

void CAF::Print()
{
  printf("Event %d:\n", event);
  //printf( "   Truth: Ev = %3.3f Elep = %3.3f Q2 = %3.3f W = %3.3f x = %3.3f y = %3.3f lepton %d mode %d\n", Ev, LepE, Q2, W, X, Y, LepPDG, mode );
  //printf( "    Reco: Ev = %3.3f Elep = %3.3f q %d mu/e/nc %d%d%d cont/trk/ecal/exit %d%d%d%d had veto %2.1f\n\n", Ev_reco, Elep_reco, reco_q, reco_numu, reco_nue, reco_nc, muon_contained, muon_tracker, muon_ecal, muon_exit, Ehad_veto );
}

void CAF::fillPOT()
{
  printf("Filling metadata\n");
  cafPOT->Fill();
}

void CAF::write()
{
  cafFile->cd();
  cafMVA->Write();
  cafPOT->Write();
  genie->Write();
  cafFile->Close();
}

void CAF::addRWbranch(int parId, std::string name, std::string wgt_var, std::vector<double>& vars)
{
  cafMVA->Branch(
    Form("%s_nshifts", name.c_str()), &nwgt[parId], Form("%s_nshifts/I", name.c_str()));
  cafMVA->Branch(Form("%s_cv%s", name.c_str(), wgt_var.c_str()),
                 &cvwgt[parId],
                 Form("%s_cv%s/D", name.c_str(), wgt_var.c_str()));
  cafMVA->Branch(Form("%s_%s", wgt_var.c_str(), name.c_str()),
                 wgt[parId],
                 Form("%s_%s[%s_nshifts]/D", wgt_var.c_str(), name.c_str(), name.c_str()));
}

#endif
