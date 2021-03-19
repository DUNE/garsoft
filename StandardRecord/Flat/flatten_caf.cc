#include "StandardRecord/StandardRecord.h"

#include "StandardRecord/Flat/FlatRecord.h"

//#include "CAFAna/Core/Progress.h"

#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

int main(int argc, char** argv)
{
  if(argc != 3){
    std::cout << "Usage: flatten_caf input.events.root output.flat.root"
              << std::endl;
    return 1;
  }

  const std::string inname = argv[1];
  const std::string outname = argv[2];

  TFile* fin = TFile::Open(inname.c_str());

  TTree* tr = (TTree*)fin->Get("anatree/GArAnaTree");
  if(!tr){
    std::cout << "Couldn't find tree 'anatree/GArAnaTree' in " << inname << std::endl;
    return 1;
  }

  caf::StandardRecord* event = 0;
  tr->SetBranchAddress("rec", &event);

  // LZ4 is the fastest format to decompress. I get 3x faster loading with
  // this compared to the default, and the files are only slightly larger.
  TFile fout(outname.c_str(), "RECREATE", "",
             ROOT::CompressionSettings(ROOT::kLZ4, 1));

  TTree* trout = new TTree("recTree", "recTree");
  // On NOvA, had trouble with memory usage (because several trees are open at
  // once?). Setting the maximum buffer usage (per tree) to 3MB (10x less than
  // default) fixed it. But it doesn't seem necessary for now on SBN.
  //  trout->SetAutoFlush(-3*1000*1000);

  flat::Flat<caf::StandardRecord> rec(trout, "rec", "", 0);//policy);

  //  ana::Progress prog("Converting '"+inname+"' to '"+outname+"'");
  for(int i = 0; i < tr->GetEntries(); ++i){
    //    prog.SetProgress(double(i)/tr->GetEntries());

    tr->GetEntry(i);

    rec.Clear();
    rec.Fill(*event);
    trout->Fill();
  }
  //  prog.Done();

  trout->Write();

  TH1* hPOT = (TH1*)fin->Get("TotalPOT");
  TH1* hEvts = (TH1*)fin->Get("TotalEvents");
  fout.cd();
  if(hPOT) hPOT->Write("TotalPOT");
  if(hEvts) hEvts->Write("TotalEvents");

  TTree* metain = (TTree*)fin->Get("metadata/metatree");
  if(metain){
    TDirectory* metadir = fout.mkdir("metadata");
    metadir->cd();
    std::string key, value;
    std::string *pkey = &key, *pvalue = &value;
    TTree* metaout = new TTree("metatree", "metatree");
    metaout->Branch("key", &key);
    metaout->Branch("value", &value);
    metain->SetBranchAddress("key", &pkey);
    metain->SetBranchAddress("value", &pvalue);
    for(int i = 0; i < metain->GetEntries(); ++i){
      metain->GetEntry(i);

      if(key == "data_tier") value = "\"flat_caf\"";
      if(key == "file_format") value = "\"flat_caf\"";

      metaout->Fill();
    }
    metaout->Write();
  }

  return 0;
}