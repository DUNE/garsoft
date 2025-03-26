///
/// \file    CalorPad.cxx
/// \brief   Drawing pad showing calorimetric particle ID information
/// \author  msoderbe@syr.edu
///

#include "EventDisplay/EVD/CalorPad.h"
#include "EventDisplay/EVD/AnalysisBaseDrawer.h"
#include "EventDisplay/EVD/AnalysisDrawingOptions.h"
#include "EventDisplay/EVD/EvdLayoutOptions.h"
#include "EventDisplay/EVD/RecoDrawingOptions.h"
#include "EventDisplay/EVD/Style.h"
#include "nuevdb/EventDisplayBase/EventHolder.h"
#include "nuevdb/EventDisplayBase/View2D.h"
#include <iostream>

#include "TFile.h"
#include "TFrame.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TPad.h"

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "cetlib/search_path.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
///
/// Create a pad to show calorimety/PID info. for reconstructed tracks.
/// @param name : Name of the pad
/// @param title : Title of the pad
/// @param x1 : Location of left  edge of pad (0-1)
/// @param x2 : Location of right edge of pad (0-1)
/// @param y1 : Location of bottom edge of pad (0-1)
/// @param y2 : Location of top    edge of pad (0-1)
///

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn, cet::exception const& e)
  {
    mf::LogWarning("CalorPad") << "CalorPad::" << fcn << " failed with message:\n" << e;
  }
}

//------------------------------------------------------------------------------
gar::evd::CalorPad::CalorPad(const char* name,
                             const char* title,
                             double x1,
                             double y1,
                             double x2,
                             double y2,
                             int curvetype)
  : DrawingPad(name, title, x1, y1, x2, y2), fcurvetype(curvetype)
{

  // Set up pad.
  this->Pad()->cd();
  this->Pad()->SetBit(kCannotPick);
  this->Pad()->SetBit(TPad::kCannotMove);
  this->Pad()->SetFillColor(kWhite);
  this->Pad()->SetLeftMargin(0.10);
  this->Pad()->SetRightMargin(0.025);
  this->Pad()->SetTopMargin(0.025);
  this->Pad()->SetBottomMargin(0.10);
  this->Pad()->Draw();

  dedx_range_pro = 0;
  dedx_range_ka = 0;
  dedx_range_pi = 0;
  dedx_range_mu = 0;
  ke_range_pro = 0;
  ke_range_ka = 0;
  ke_range_pi = 0;
  ke_range_mu = 0;

  fView = new evdb::View2D();
}

//......................................................................
// Destructor.
gar::evd::CalorPad::~CalorPad()
{
  if (dedx_range_pro) {
    delete dedx_range_pro;
    dedx_range_pro = 0;
  }
  if (dedx_range_ka) {
    delete dedx_range_ka;
    dedx_range_ka = 0;
  }
  if (dedx_range_pi) {
    delete dedx_range_pi;
    dedx_range_pi = 0;
  }
  if (dedx_range_mu) {
    delete dedx_range_mu;
    dedx_range_mu = 0;
  }
  if (ke_range_pro) {
    delete ke_range_pro;
    ke_range_pro = 0;
  }
  if (ke_range_ka) {
    delete ke_range_ka;
    ke_range_ka = 0;
  }
  if (ke_range_pi) {
    delete ke_range_pi;
    ke_range_pi = 0;
  }
  if (ke_range_mu) {
    delete ke_range_mu;
    ke_range_mu = 0;
  }
  if (fView) {
    delete fView;
    fView = 0;
  }
}

//......................................................................
// Draw selected objects.

void gar::evd::CalorPad::Draw(const char* /*opt*/)
{

  this->Pad()->cd();

  //Remove all previous objects from Pad's primitive list
  this->Pad()->Clear();

  //Remove all previous TPolyMarkers, TLatexs, etc... from list of such objects
  fView->Clear();

  //Draw coordinate axis and also GEANT based dE/dx vs. Range, or KE vs. Range, curves.
  DrawRefCurves();

  // grab the event from the singleton
  const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
  art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;

  // Insert graphic objects into fView collection.
  if (evt) {
    try {
      if (fcurvetype == 1)
        AnalysisBaseDraw()->DrawDeDx(*evt, fView);
      else if (fcurvetype == 0)
        AnalysisBaseDraw()->DrawKineticEnergy(*evt, fView);
    }
    catch (cet::exception& e) {
      if (fcurvetype == 1)
        writeErrMsg("Draw->DrawDeDx", e);
      else if (fcurvetype == 0)
        writeErrMsg("Draw->DrawKineticEnergy", e);
      else if (fcurvetype == 2)
        writeErrMsg("Draw->CalorShower", e);
    }
  }
  //      try{
  //	AnalysisBaseDraw()->CalorShower(*evt, fView);
  //      }
  //      catch (cet::exception e){
  //        writeErrMsg("Draw->CalorShower",e);
  //      }

  // Draw objects on pad.
  fView->Draw();
  fPad->Modified();
  fPad->Update();
}

//......................................................................
// Draw truth curves

void gar::evd::CalorPad::DrawRefCurves()
{

  if (dedx_range_pro) {
    delete dedx_range_pro;
    dedx_range_pro = 0;
  }
  if (dedx_range_ka) {
    delete dedx_range_ka;
    dedx_range_ka = 0;
  }
  if (dedx_range_pi) {
    delete dedx_range_pi;
    dedx_range_pi = 0;
  }
  if (dedx_range_mu) {
    delete dedx_range_mu;
    dedx_range_mu = 0;
  }
  if (ke_range_pro) {
    delete ke_range_pro;
    ke_range_pro = 0;
  }
  if (ke_range_ka) {
    delete ke_range_ka;
    ke_range_ka = 0;
  }
  if (ke_range_pi) {
    delete ke_range_pi;
    ke_range_pi = 0;
  }
  if (ke_range_mu) {
    delete ke_range_mu;
    ke_range_mu = 0;
  }

  double ymax;
  if (fcurvetype == 1)
    ymax = 50.0;
  else
    ymax = 200.0;
  TH1F* h = this->Pad()->DrawFrame(0.0, 0.0, 25.0, ymax);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->CenterTitle();

  if (fcurvetype == 1) {
    h->GetXaxis()->SetTitle("Residual Range (cm)");
    h->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  }
  else {
    h->GetXaxis()->SetTitle("Total Range (cm)");
    h->GetYaxis()->SetTitle("T (MeV)");
  }

  art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;

  cet::search_path sp("FW_SEARCH_PATH");
  if (!sp.find_file(anaOpt->fCalorTemplateFileName + ".root", fROOTfile))
    throw cet::exception("Chi2ParticleID")
      << "cannot find the root template file: \n"
      << anaOpt->fCalorTemplateFileName << "\n bail ungracefully.\n";

  TFile* file = TFile::Open(fROOTfile.c_str());
  if (fcurvetype == 1) {
    dedx_range_pro = (TGraph*)file->Get("dedx_range_pro");
    dedx_range_ka = (TGraph*)file->Get("dedx_range_ka");
    dedx_range_pi = (TGraph*)file->Get("dedx_range_pi");
    dedx_range_mu = (TGraph*)file->Get("dedx_range_mu");

    dedx_range_pro->SetMarkerStyle(7);
    dedx_range_ka->SetMarkerStyle(7);
    dedx_range_pi->SetMarkerStyle(7);
    dedx_range_mu->SetMarkerStyle(7);

    dedx_range_pro->SetMarkerColor(kBlack);
    dedx_range_ka->SetMarkerColor(kGray + 2);
    dedx_range_pi->SetMarkerColor(kGray + 1);
    dedx_range_mu->SetMarkerColor(kGray);

    dedx_range_mu->Draw("P,same");
    dedx_range_pi->Draw("P,same");
    dedx_range_ka->Draw("P,same");
    dedx_range_pro->Draw("P,same");
  }
  else {
    ke_range_pro = (TGraph*)file->Get("kinen_range_pro");
    ke_range_ka = (TGraph*)file->Get("kinen_range_ka");
    ke_range_pi = (TGraph*)file->Get("kinen_range_pi");
    ke_range_mu = (TGraph*)file->Get("kinen_range_mu");

    ke_range_pro->SetMarkerStyle(7);
    ke_range_ka->SetMarkerStyle(7);
    ke_range_pi->SetMarkerStyle(7);
    ke_range_mu->SetMarkerStyle(7);

    ke_range_pro->SetMarkerColor(kBlack);
    ke_range_ka->SetMarkerColor(kGray + 2);
    ke_range_pi->SetMarkerColor(kGray + 1);
    ke_range_mu->SetMarkerColor(kGray);

    ke_range_mu->Draw("P,same");
    ke_range_pi->Draw("P,same");
    ke_range_ka->Draw("P,same");
    ke_range_pro->Draw("P,same");
  }
  file->Close();
}

////////////////////////////////////////////////////////////////////////
