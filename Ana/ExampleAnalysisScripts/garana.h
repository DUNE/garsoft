//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 22 17:38:47 2018 by ROOT version 6.12/06
// from TTree GArAnaTree/GArAnaTree
// found on file: anatree.root
//////////////////////////////////////////////////////////

#ifndef garana_h
#define garana_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class garana {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   Int_t           SubRun;
   Int_t           Run;
   vector<int>     *NType;
   vector<int>     *CCNC;
   vector<int>     *PDG;
   vector<float>   *MCPStartX;
   vector<float>   *MCPStartY;
   vector<float>   *MCPStartZ;
   vector<float>   *MCPPX;
   vector<float>   *MCPPY;
   vector<float>   *MCPPZ;
   vector<float>   *HitX;
   vector<float>   *HitY;
   vector<float>   *HitZ;
   vector<float>   *HitSig;
   vector<float>   *HitRMS;
   vector<float>   *TrackStartX;
   vector<float>   *TrackStartY;
   vector<float>   *TrackStartZ;
   vector<float>   *TrackStartPX;
   vector<float>   *TrackStartPY;
   vector<float>   *TrackStartPZ;
   vector<float>   *TrackEndX;
   vector<float>   *TrackEndY;
   vector<float>   *TrackEndZ;
   vector<float>   *TrackEndPX;
   vector<float>   *TrackEndPY;
   vector<float>   *TrackEndPZ;
   vector<float>   *VertX;
   vector<float>   *VertY;
   vector<float>   *VertZ;
   vector<int>     *VertN;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_SubRun;   //!
   TBranch        *b_Run;   //!
   TBranch        *b_NType;   //!
   TBranch        *b_CCNC;   //!
   TBranch        *b_PDG;   //!
   TBranch        *b_MCPStartX;   //!
   TBranch        *b_MCPStartY;   //!
   TBranch        *b_MCPStartZ;   //!
   TBranch        *b_MCPPX;   //!
   TBranch        *b_MCPPY;   //!
   TBranch        *b_MCPPZ;   //!
   TBranch        *b_HitX;   //!
   TBranch        *b_HitY;   //!
   TBranch        *b_HitZ;   //!
   TBranch        *b_HitSig;   //!
   TBranch        *b_HitRMS;   //!
   TBranch        *b_TrackStartX;   //!
   TBranch        *b_TrackStartY;   //!
   TBranch        *b_TrackStartZ;   //!
   TBranch        *b_TrackStartPX;   //!
   TBranch        *b_TrackStartPY;   //!
   TBranch        *b_TrackStartPZ;   //!
   TBranch        *b_TrackEndX;   //!
   TBranch        *b_TrackEndY;   //!
   TBranch        *b_TrackEndZ;   //!
   TBranch        *b_TrackEndPX;   //!
   TBranch        *b_TrackEndPY;   //!
   TBranch        *b_TrackEndPZ;   //!
   TBranch        *b_VertX;   //!
   TBranch        *b_VertY;   //!
   TBranch        *b_VertZ;   //!
   TBranch        *b_VertN;   //!

   garana(TTree *tree=0);
   virtual ~garana();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef garana_cxx
garana::garana(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("anatree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("anatree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("anatree.root:/anatree");
      dir->GetObject("GArAnaTree",tree);

   }
   Init(tree);
}

garana::~garana()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t garana::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t garana::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void garana::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   NType = 0;
   CCNC = 0;
   PDG = 0;
   MCPStartX = 0;
   MCPStartY = 0;
   MCPStartZ = 0;
   MCPPX = 0;
   MCPPY = 0;
   MCPPZ = 0;
   HitX = 0;
   HitY = 0;
   HitZ = 0;
   HitSig = 0;
   HitRMS = 0;
   TrackStartX = 0;
   TrackStartY = 0;
   TrackStartZ = 0;
   TrackStartPX = 0;
   TrackStartPY = 0;
   TrackStartPZ = 0;
   TrackEndX = 0;
   TrackEndY = 0;
   TrackEndZ = 0;
   TrackEndPX = 0;
   TrackEndPY = 0;
   TrackEndPZ = 0;
   VertX = 0;
   VertY = 0;
   VertZ = 0;
   VertN = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("SubRun", &SubRun, &b_SubRun);
   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("NType", &NType, &b_NType);
   fChain->SetBranchAddress("CCNC", &CCNC, &b_CCNC);
   fChain->SetBranchAddress("PDG", &PDG, &b_PDG);
   fChain->SetBranchAddress("MCPStartX", &MCPStartX, &b_MCPStartX);
   fChain->SetBranchAddress("MCPStartY", &MCPStartY, &b_MCPStartY);
   fChain->SetBranchAddress("MCPStartZ", &MCPStartZ, &b_MCPStartZ);
   fChain->SetBranchAddress("MCPPX", &MCPPX, &b_MCPPX);
   fChain->SetBranchAddress("MCPPY", &MCPPY, &b_MCPPY);
   fChain->SetBranchAddress("MCPPZ", &MCPPZ, &b_MCPPZ);
   fChain->SetBranchAddress("HitX", &HitX, &b_HitX);
   fChain->SetBranchAddress("HitY", &HitY, &b_HitY);
   fChain->SetBranchAddress("HitZ", &HitZ, &b_HitZ);
   fChain->SetBranchAddress("HitSig", &HitSig, &b_HitSig);
   fChain->SetBranchAddress("HitRMS", &HitRMS, &b_HitRMS);
   fChain->SetBranchAddress("TrackStartX", &TrackStartX, &b_TrackStartX);
   fChain->SetBranchAddress("TrackStartY", &TrackStartY, &b_TrackStartY);
   fChain->SetBranchAddress("TrackStartZ", &TrackStartZ, &b_TrackStartZ);
   fChain->SetBranchAddress("TrackStartPX", &TrackStartPX, &b_TrackStartPX);
   fChain->SetBranchAddress("TrackStartPY", &TrackStartPY, &b_TrackStartPY);
   fChain->SetBranchAddress("TrackStartPZ", &TrackStartPZ, &b_TrackStartPZ);
   fChain->SetBranchAddress("TrackEndX", &TrackEndX, &b_TrackEndX);
   fChain->SetBranchAddress("TrackEndY", &TrackEndY, &b_TrackEndY);
   fChain->SetBranchAddress("TrackEndZ", &TrackEndZ, &b_TrackEndZ);
   fChain->SetBranchAddress("TrackEndPX", &TrackEndPX, &b_TrackEndPX);
   fChain->SetBranchAddress("TrackEndPY", &TrackEndPY, &b_TrackEndPY);
   fChain->SetBranchAddress("TrackEndPZ", &TrackEndPZ, &b_TrackEndPZ);
   fChain->SetBranchAddress("VertX", &VertX, &b_VertX);
   fChain->SetBranchAddress("VertY", &VertY, &b_VertY);
   fChain->SetBranchAddress("VertZ", &VertZ, &b_VertZ);
   fChain->SetBranchAddress("VertN", &VertN, &b_VertN);
   Notify();
}

Bool_t garana::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void garana::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t garana::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef garana_cxx
