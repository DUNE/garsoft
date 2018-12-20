//
//  geovis.C
//  garsoft-mrb
//
//  Original from Brian Rebel
//  
//
#include "TGeoManager.h"
#include "TFile.h"
#include "TSystem.h"

typedef struct _drawopt {
  const char* volume;
  int         color;
} drawopt;

//void geoVis(TString volName="volWorld",TString filebase="hpgtpc_v2", bool checkoverlaps=true, bool writerootfile=true){

void geoVis(TString volName="volWorld",TString filebase="ND_Strawman_Concept_v01", bool checkoverlaps=false, bool writerootfile=false){

//  void geoVis(TString volName="volWorld",TString filebase="ND_Concept_all", bool checkoverlaps=false, bool writerootfile=false){
  
  gSystem->Load("libGeom");
  gSystem->Load("libGdml");
  
  TString gdmlfile=filebase;
  gdmlfile += ".gdml";
  TGeoManager::Import(gdmlfile);
  
  drawopt opt[] = {
    {"volWorld",           0},
    {"volDetEnclosure",    kWhite},
    {"volCryostat",        kOrange},
    {"volTPCWidthFace",    kRed},
    {"volTPCLengthFace",   kCyan+5},
    {"volTPCBottomFace",   kOrange},
    {"volTubBottom",       kOrange+7},
    {"voltheX",            kOrange+7},
    {"volArgon_solid_L",   kRed},
    {"volArgon_cap_L",     kOrange},
    {"volArgon_cap_front", kOrange},
    {"volTPC",             kOrange},
    {"TPCGas_vol",         kOrange},
    {"volDetEnclosure",    kBlue},
    {"volTPCActive",       kGreen},
    {0, 0}
  };
  
  /*for (int i=0;; ++i) {
    if (opt[i].volume==0) break;
    gGeoManager->FindVolumeFast(opt[i].volume)->SetLineColor(opt[i].color);
    }*/
  
  // TList* mat = gGeoManager->GetListOfMaterials();
  // TIter next(mat);
  // TObject *obj;
  // while (obj = next()) {
  //  obj->Print();
  // }

  if (checkoverlaps)
    {
      gGeoManager->CheckOverlaps(0.01);
      gGeoManager->PrintOverlaps();
    }
  gGeoManager->SetMaxVisNodes(70000);

  TObjArray *vollist = gGeoManager->GetListOfVolumes();
  TIter next(vollist);
  TObject *obj = next();
  while (obj) {
    obj->Print();
    obj = next();
   }

  //gGeoManager->GetTopVolume()->Draw();
  gGeoManager->FindVolumeFast(volName)->Draw("ogl");
  //gGeoManager->FindVolumeFast(volName)->Draw("");

  if (writerootfile)
    {
      TString rootfile=filebase;
      rootfile += ".root";
      TFile *tf = new TFile(rootfile, "RECREATE");
      gGeoManager->Write();
      tf->Close();
    }
}
