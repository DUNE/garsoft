#include "EVGCore/EventRecord.h"
#include "NUSYS.C"
#include "Ntuple/NtpMCEventRecord.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TVector3.h"
#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>

TRandom3* rando;

struct params {
  double OA_xcoord;
  bool fhc, grid;
  int seed, run, subrun, first, n, nfiles;
};

// main loop function
void loop(CAF& caf, params& par, TTree* tree, std::string ghepdir, std::string fhicl_filename)
{

  // Get GHEP file for genie::EventRecord from other file
  int ievt, ifileNo;
  tree->SetBranchAddress("ievt", &ievt);
  tree->SetBranchAddress("ifileNo", &ifileNo);
  int current_file = -1;
  TFile* ghep_file = NULL;
  TTree* gtree = NULL;
  std::string mode = (par.fhc ? "neutrino" : "antineutrino");

  // DUNE reweight getter
  nusyst::response_helper rh(fhicl_filename);
  // Get list of variations, and make tree branch for each one
  std::vector<unsigned int> parIds = rh.GetParameters();
  for (unsigned int i = 0; i < parIds.size(); ++i) {
    systtools::SystParamHeader head = rh.GetHeader(parIds[i]);
    printf("Adding reweight branch %u for %s with %lu shifts\n",
           parIds[i],
           head.prettyName.c_str(),
           head.paramVariations.size());
    bool is_wgt = head.isWeightSystematicVariation;
    std::string wgt_var = (is_wgt ? "wgt" : "var");
    caf.addRWbranch(parIds[i], head.prettyName, wgt_var, head.paramVariations);
    caf.iswgt[parIds[i]] = is_wgt;
  }

  // Main event loop
  int N = tree->GetEntries();
  if (par.n > 0 && par.n < N) N = par.n + par.first;
  for (int ii = par.first; ii < N; ++ii) {

    tree->GetEntry(ii);
    if (ii % 100 == 0) printf("Event %d of %d...\n", ii, N);

    // set defaults
    for (int j = 0; j < 100; ++j) {
      caf.nwgt[j] = 7;
      caf.cvwgt[j] = (caf.iswgt[j] ? 1. : 0.);
      for (unsigned int k = 0; k < 100; ++k) {
        caf.wgt[j][k] = (caf.iswgt[j] ? 1. : 0.);
      }
    }

    // make sure ghep file matches the current one, otherwise update to the current ghep file
    if (ifileNo != current_file) {
      // close the previous file
      if (ghep_file) ghep_file->Close();

      if (par.grid)
        ghep_file = new TFile(Form("genie.%d.root", ifileNo));
      else
        ghep_file =
          new TFile(Form("%s/GAr.%s.%d.ghep.root", ghepdir.c_str(), mode.c_str(), ifileNo));

      gtree = (TTree*)ghep_file->Get("gtree");

      gtree->SetBranchAddress("gmcrec", &caf.mcrec);
      current_file = ifileNo;
    }

    // get GENIE event record

    gtree->GetEntry(ievt);
    genie::EventRecord* event = caf.mcrec->event;
    genie::Interaction* in = event->Summary();

    // Get truth stuff out of GENIE ghep record

    TLorentzVector lepP4;
    TLorentzVector nuP4nuc = *(in->InitState().GetProbeP4(genie::kRfHitNucRest));
    TLorentzVector nuP4 = *(in->InitState().GetProbeP4(genie::kRfLab));

    // Add DUNErw weights to tree

    systtools::event_unit_response_w_cv_t resp = rh.GetEventVariationAndCVResponse(*event);
    for (systtools::event_unit_response_w_cv_t::iterator it = resp.begin(); it != resp.end();
         ++it) {
      caf.nwgt[(*it).pid] = (*it).responses.size();
      caf.cvwgt[(*it).pid] = (*it).CV_response;
      for (unsigned int i = 0; i < (*it).responses.size(); ++i) {
        caf.wgt[(*it).pid][i] = (*it).responses[i];
      }
    }

    caf.fill();
  }

  // set POT
  caf.meta_run = par.run;
  caf.meta_subrun = par.subrun;
}

int main(int argc, char const* argv[])
{

  if ((argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1]))) {
    std::cout << "Help yourself by looking at the source code to see what the options are."
              << std::endl;
    return 0;
  }

  // get command line options
  std::string ghepdir;
  std::string outfile;
  std::string edepfile;
  std::string fhicl_filename;

  // Make parameter object and set defaults
  params par;
  par.OA_xcoord = 0.; // on-axis by default
  par.fhc = true;
  par.grid = false;
  par.seed = 7; // a very random number
  par.run = 1;  // CAFAna doesn't like run number 0
  par.subrun = 0;
  par.n = -1;
  par.nfiles = 1;
  par.first = 0;

  int i = 0;
  while (i < argc) {
    if (argv[i] == std::string("--edepfile")) {
      edepfile = argv[i + 1];
      i += 2;
    }
    else if (argv[i] == std::string("--ghepdir")) {
      ghepdir = argv[i + 1];
      i += 2;
    }
    else if (argv[i] == std::string("--outfile")) {
      outfile = argv[i + 1];
      i += 2;
    }
    else if (argv[i] == std::string("--fhicl")) {
      fhicl_filename = argv[i + 1];
      i += 2;
    }
    else if (argv[i] == std::string("--seed")) {
      par.seed = atoi(argv[i + 1]);
      par.run = par.seed;
      i += 2;
    }
    else if (argv[i] == std::string("--nevents")) {
      par.n = atoi(argv[i + 1]);
      i += 2;
    }
    else if (argv[i] == std::string("--nfiles")) {
      par.nfiles = atoi(argv[i + 1]);
      i += 2;
    }
    else if (argv[i] == std::string("--first")) {
      par.first = atoi(argv[i + 1]);
      i += 2;
    }
    else if (argv[i] == std::string("--oa")) {
      par.OA_xcoord = atof(argv[i + 1]);
      i += 2;
    }
    else if (argv[i] == std::string("--grid")) {
      par.grid = true;
      i += 1;
    }
    else if (argv[i] == std::string("--rhc")) {
      par.fhc = false;
      i += 1;
    }
    else
      i += 1; // look for next thing
  }

  printf("Making CAF from edep-sim tree dump: %s\n", edepfile.c_str());
  printf("Searching for GENIE ghep files here: %s\n", ghepdir.c_str());
  if (par.fhc)
    printf("Running neutrino mode (FHC)\n");
  else
    printf("Running antineutrino mode (RHC)\n");
  if (par.grid)
    printf("Running grid mode\n");
  else
    printf("Running test mode\n");
  printf("Output CAF file: %s\n", outfile.c_str());

  rando = new TRandom3(par.seed);
  CAF caf(outfile);

  TFile* tf = new TFile(edepfile.c_str());
  TTree* tree = (TTree*)tf->Get("tree");

  loop(caf, par, tree, ghepdir, fhicl_filename);

  caf.version = 4;
  printf("Run %d POT %g\n", caf.meta_run, caf.pot);
  caf.fillPOT();
  caf.write();

  caf.cafFile->Close();

  printf("-30-\n");
}
