///
/// \file    Display3DPad.cxx
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
/// \version $Id: Display3DPad.cxx,v 1.5 2011/04/12 22:07:16 bckhouse Exp $
///
#include <iostream>
#include "TPad.h"
#include "TView3D.h"
#include "TGLViewer.h"
#include "TPolyLine3D.h"

#include "EventDisplay/Display3DPad.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "Geometry/Geometry.h"
#include "EventDisplay/HeaderDrawer.h"
#include "EventDisplay/GeometryDrawer.h"
#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplay/SimulationDrawer.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplay/EvdLayoutOptions.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace gar {
namespace evd{

  ///
  /// Create a pad to show a 3D rendering of the detector and events
  /// @param nm : Name of the pad
  /// @param ti : Title of the pad
  /// @param x1 : Location of left  edge of pad (0-1)
  /// @param x2 : Location of right edge of pad (0-1)
  /// @param y1 : Location of bottom edge of pad (0-1)
  /// @param y2 : Location of top    edge of pad (0-1)
  /// @param opt: Options. Currently just a place holder
  ///
  Display3DPad::Display3DPad(const char* nm,
                             const char* ti,
                             double x1,
                             double y1,
                             double x2,
                             double y2,
                             const char* /*opt*/)
  : DrawingPad(nm, ti, x1, y1, x2, y2)
  {
    this->Pad()->SetFillColor(kBlack);
    this->Pad()->Draw();
    this->Pad()->cd();
    fView = new evdb::View3D();
  }

  //......................................................................

  Display3DPad::~Display3DPad() 
  {
    if (fView) { delete fView; fView = 0; }
  }

  //......................................................................

  void Display3DPad::Draw() 
  {
    fView->Clear();

    art::ServiceHandle<geo::Geometry> geo;

    // grab the event from the singleton
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();

    if(evt){
      this->HeaderDraw()    ->Header       (fView);
      this->GeometryDraw()  ->DetOutline3D (fView);
      this->SimulationDraw()->MCTruth3D    (*evt, fView);
      this->RecoBaseDraw()  ->Track3D      (*evt, fView);
      this->RecoBaseDraw()  ->Hit3D        (*evt, fView);
      this->RecoBaseDraw()  ->TPCCluster3D (*evt, fView);
      this->RecoBaseDraw()  ->CaloCluster3D(*evt, fView);
      this->RecoBaseDraw()  ->Vertex3D     (*evt, fView);
      this->RecoBaseDraw()  ->VecHit3D     (*evt, fView);
      this->RawDataDraw()   ->RawDigit3D   (*evt, fView);
    }
    
    this->Pad()->Clear();
    this->Pad()->cd();
    if (fPad->GetView()==0) {
      //int irep=0;
      double rmin[]={geo->TPCXCent() - 0.7*geo->DetHalfWidth(), geo->TPCYCent() - 0.7*geo->DetHalfHeight(), geo->TPCZCent() - 0.7*geo->DetLength()};
      double rmax[]={geo->TPCXCent() + 0.7*geo->DetHalfWidth(), geo->TPCYCent() + 0.7*geo->DetHalfHeight(), geo->TPCZCent() + 0.7*geo->DetLength()};
      TView3D* v = new TView3D(1,rmin,rmax);
      v->SetPerspective();
      //v->SetView(-90.0,75.0,0,irep);
      //v->ZoomView(0, 3);
      v->ResizePad();
      //v->ToggleZoom(0);
      //v->SetView(geo->TPCXCent(),geo->TPCYCent(),geo->TPCZCent(),irep);
      fPad->SetView(v); // ROOT takes ownership of object *v
    }
    fView->Draw();
    fPad->Update();
  }


  ////////////////////////////////////////////////////////////////////////
}
}//namespace
