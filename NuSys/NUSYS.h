#ifndef CAF_h
#define CAF_h

#include "Ntuple/NtpMCEventRecord.h"
#include "TFile.h"
#include "TTree.h"

class CAF {

public:
  CAF(std::string filename, bool isGas = false);
  ~CAF();
  void fill();
  void fillPOT();
  void write();
  void addRWbranch(int parId, std::string name, std::string wgt_var, std::vector<double>& vars);
  void Print();
  void setToBS();

  // Make ntuple variables public so they can be set from other file

  // configuration variables
  int isFD, isFHC;
  // event accounting
  int run, subrun, event, ievt;

  // reweights -- make sure big enough to hold all the variations for each knob, and all the knobs
  // the names, and what they actually mean, are determined automatically from the fhicl input file
  int nwgt[100];
  double cvwgt[100];
  double wgt[100][100];
  bool iswgt[100];

  // store the GENIE record as a branch
  genie::NtpMCEventRecord* mcrec;

  // meta
  double pot;
  int meta_run, meta_subrun;
  int version;

  TFile* cafFile;
  TTree* cafMVA;
  TTree* cafPOT;
  TTree* genie;
};

#endif
