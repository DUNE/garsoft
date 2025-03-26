///
/// \file    Display3DView.cxx
/// \brief   The "main" event display view that most people will want to use
/// \author  messier@indiana.edu
/// \version $Id: Display3DView.cxx,v 1.2 2011/01/19 16:40:59 p-novaart Exp $
///
#include "TCanvas.h"
#include "TVirtualViewer3D.h"

#include "EventDisplay/EVD/ColorDrawingOptions.h"
#include "EventDisplay/EVD/Display3DPad.h"
#include "EventDisplay/EVD/Display3DView.h"
#include "EventDisplay/EVD/EvdLayoutOptions.h"
#include "EventDisplay/EVD/RawDataDrawer.h"
#include "EventDisplay/EVD/RawDrawingOptions.h"
#include "EventDisplay/EVD/RecoBaseDrawer.h"
#include "EventDisplay/EVD/RecoDrawingOptions.h"
#include "EventDisplay/EVD/SimulationDrawingOptions.h"

#include "art/Framework/Principal/Event.h"

namespace gar {
  namespace evd {

    //......................................................................
    Display3DView::Display3DView(TGMainFrame* mf)
      : evdb::Canvas(mf), fMCOn(nullptr), fRawDraw(nullptr), fRecoDraw(nullptr)
    {
      art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutopt;
      art::ServiceHandle<evd::SimulationDrawingOptions> sdo;

      fHeaderPad = new HeaderPad("fHeaderPad", "Header", 0.00, 0.00, 0.15, 0.13, "");
      evdb::Canvas::fCanvas->cd();

      fMC = new MCBriefPad("fMCPad", "MC Info.", 0.15, 0.13, 1.00, 0.17, "");

      evdb::Canvas::fCanvas->cd();
      fHeaderPad->Draw();

      evdb::Canvas::fCanvas->cd();
      fMC->Draw();

      if (evdlayoutopt->fEnableMCTruthCheckBox) {
        fMCOn = new TGCheckButton(fFrame, "MC Truth", 5);
        fMCOn->Connect("Clicked()", "gar::evd::Display3DView", this, "SetMCInfo()");
        if (sdo->fShowMCTruthText == 1) fMCOn->SetState(kButtonDown);
      }

      // radio buttons to toggle drawing raw vs calibrated information
      fRecoDraw = new TGRadioButton(fFrame, "Reconstructed", 3);
      fRawDraw = new TGRadioButton(fFrame, "Raw", 4);
      fRawDraw->Connect("Clicked()", "gar::evd::Display3DView", this, "SetRawReco()");
      fRecoDraw->Connect("Clicked()", "gar::evd::Display3DView", this, "SetRawReco()");

      if (evdlayoutopt->fEnableMCTruthCheckBox) {
        fFrame->AddFrame(fMCOn, new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 5, 1));
      }
      fFrame->AddFrame(fRecoDraw, new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 5, 1));
      fFrame->AddFrame(fRawDraw, new TGLayoutHints(kLHintsBottom | kLHintsRight, 0, 0, 5, 1));
      evdb::Canvas::fCanvas->cd();
      fDisplay3DPad = new Display3DPad("fDisplay3DPad", "3D Display", 0.0, 0.0, 1.0, 1.0, "");

      this->Connect("CloseWindow()", "gar::evd::Display3DView", this, "CloseWindow()");

      fDisplay3DPad->Draw();

      evdb::Canvas::fCanvas->Update();
    }

    //......................................................................
    Display3DView::~Display3DView()
    {
      if (fHeaderPad) {
        delete fHeaderPad;
        fHeaderPad = nullptr;
      }
      if (fMC) {
        delete fMC;
        fMC = nullptr;
      }
    }

    //......................................................................
    void Display3DView::CloseWindow()
    {
      delete this;
    }

    //......................................................................
    void Display3DView::Draw(const char* /*opt*/)
    {

      fDisplay3DPad->Draw();
      evdb::Canvas::fCanvas->Update();

      TVirtualViewer3D* viewer = fDisplay3DPad->Pad()->GetViewer3D("pad");
      viewer->PreferLocalFrame();
      viewer->ResetCameras();
      viewer->PadPaint(fDisplay3DPad->Pad());
    }

    //......................................................................
    void Display3DView::SetRawReco()
    {
      art::ServiceHandle<evd::RawDrawingOptions> rawopt;

      TGButton* b = (TGButton*)gTQSender;
      int id = b->WidgetId();

      // id values are set in lines 41 - 42
      if (id == 4) {
        rawopt->fDrawRawOrReco = 0;
        fRawDraw->SetState(kButtonDown);
        fRecoDraw->SetState(kButtonUp);
      }
      else if (id == 3) {
        rawopt->fDrawRawOrReco = 1;
        fRawDraw->SetState(kButtonUp);
        fRecoDraw->SetState(kButtonDown);
      }

      TVirtualPad* ori = gPad;

      evdb::Canvas::fCanvas->cd();
      evdb::Canvas::fCanvas->Modified();
      evdb::Canvas::fCanvas->Update();

      ori->cd();

      return;
    }

    //......................................................................
    void Display3DView::SetMCInfo()
    {
      art::ServiceHandle<evd::SimulationDrawingOptions> simopt;

      //TGButton *b = (TGButton *)gTQSender;
      //int id = b->WidgetId();
      // set button states TODO

      TVirtualPad* ori = gPad;

      evdb::Canvas::fCanvas->cd();
      evdb::Canvas::fCanvas->Modified();
      evdb::Canvas::fCanvas->Update();

      ori->cd();

      return;
    }
  }
} // end namespace
////////////////////////////////////////////////////////////////////////
