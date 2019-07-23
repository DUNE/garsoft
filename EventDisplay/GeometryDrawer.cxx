/// \file    GeometryDrawer.cxx
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
/// \version $Id: GeometryDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TMath.h"

#include "Geometry/Geometry.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "EventDisplay/GeometryDrawer.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace gar {
  namespace evd{

    //......................................................................
    GeometryDrawer::GeometryDrawer()
    {
    }

    //......................................................................
    GeometryDrawer::~GeometryDrawer()
    {
    }

    //......................................................................
    void GeometryDrawer::DetOutline3D(evdb::View3D*        view)
    {
      art::ServiceHandle<geo::Geometry> geo;

      // these are relative to the center of the TPC.  We will add the detector offsets
      // on each point

      double xlo =  -geo->TPCLength()/2;
      double xhi =  geo->TPCLength()/2;
      double ylo = -geo->TPCHalfHeight();
      double yhi =  geo->TPCHalfHeight();
      //double zlo =  -geo->TPCHalfWidth();
      double zhi =  geo->TPCHalfWidth();
      double r=yhi;

      double fracinner = 788.0/2580.;

      double ang = TMath::Pi()*2.0/18.0;

      int c = kGray;
      int s = 1;
      int w = 1;
      TPolyLine3D& spos = view->AddPolyLine3D(19, c, w, s);
      TPolyLine3D& sposi = view->AddPolyLine3D(19, c, w, s);
      for (int i=0;i<19;++i)
      {
        spos.SetPoint(i,xhi+geo->TPCXCent(),r*TMath::Sin(i*ang)+geo->TPCYCent(),r*TMath::Cos(i*ang)+geo->TPCZCent());
        sposi.SetPoint(i,xhi+geo->TPCXCent(),fracinner*r*TMath::Sin(i*ang)+geo->TPCYCent(),fracinner*r*TMath::Cos(i*ang)+geo->TPCZCent());
      }
      TPolyLine3D& sneg = view->AddPolyLine3D(19, c, w, s);
      TPolyLine3D& snegi = view->AddPolyLine3D(19, c, w, s);
      for (int i=0;i<19;++i)
      {
        sneg.SetPoint(i,xlo+geo->TPCXCent(),r*TMath::Sin(i*ang)+geo->TPCYCent(),r*TMath::Cos(i*ang)+geo->TPCZCent());
        snegi.SetPoint(i,xlo+geo->TPCXCent(),fracinner*r*TMath::Sin(i*ang)+geo->TPCYCent(),fracinner*r*TMath::Cos(i*ang)+geo->TPCZCent());
      }

      c = kGray+2;
      s = 1;
      w = 1;
      for (int i=0;i<18;++i)
      {
        TPolyLine3D& gridt = view->AddPolyLine3D(4, c, w, s);
        gridt.SetPoint(0,xlo+geo->TPCXCent(),fracinner*r*TMath::Sin(i*ang)+geo->TPCYCent(),fracinner*r*TMath::Cos(i*ang)+geo->TPCZCent());
        gridt.SetPoint(1,xlo+geo->TPCXCent(),r*TMath::Sin(i*ang)+geo->TPCYCent(),r*TMath::Cos(i*ang)+geo->TPCZCent());
        gridt.SetPoint(2,xhi+geo->TPCXCent(),r*TMath::Sin(i*ang)+geo->TPCYCent(),r*TMath::Cos(i*ang)+geo->TPCZCent());
        gridt.SetPoint(3,xhi+geo->TPCXCent(),fracinner*r*TMath::Sin(i*ang)+geo->TPCYCent(),fracinner*r*TMath::Cos(i*ang)+geo->TPCZCent());
      }

      // Indicate coordinate system
      double x0 = -0.20;     // Center location of the key
      double y0 =  1.10*ylo; // Center location of the key
      double z0 = -0.10*zhi; // Center location of the key
      double sz =  0.20*zhi; // Scale size of the key in z direction

      c = kCyan;

      TPolyLine3D& xaxis = view->AddPolyLine3D(2, c, w, s);
      TPolyLine3D& yaxis = view->AddPolyLine3D(2, c, w, s);
      TPolyLine3D& zaxis = view->AddPolyLine3D(2, c, s, w);
      xaxis.SetPoint(0, x0+geo->TPCXCent(),    y0+geo->TPCYCent(), z0+geo->TPCZCent());
      xaxis.SetPoint(1, sz+x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0+geo->TPCZCent());

      yaxis.SetPoint(0, x0+geo->TPCXCent(), y0+geo->TPCYCent(),     z0+geo->TPCZCent());
      yaxis.SetPoint(1, x0+geo->TPCXCent(), y0+sz+geo->TPCYCent(),  z0+geo->TPCZCent());

      zaxis.SetPoint(0, x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0+geo->TPCZCent());
      zaxis.SetPoint(1, x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0+sz+geo->TPCZCent());

      TPolyLine3D& xpoint = view->AddPolyLine3D(3, c, w, s);
      TPolyLine3D& ypoint = view->AddPolyLine3D(3, c, w, s);
      TPolyLine3D& zpoint = view->AddPolyLine3D(3, c, w, s);

      xpoint.SetPoint(0, 0.95*sz+x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0-0.05*sz+geo->TPCZCent());
      xpoint.SetPoint(1, 1.00*sz+x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0+geo->TPCZCent());
      xpoint.SetPoint(2, 0.95*sz+x0+geo->TPCXCent(), y0+geo->TPCYCent(), z0+0.05*sz+geo->TPCZCent());

      ypoint.SetPoint(0, x0+geo->TPCXCent(), 0.95*sz+y0+geo->TPCYCent(), z0-0.05*sz+geo->TPCZCent());
      ypoint.SetPoint(1, x0+geo->TPCXCent(), 1.00*sz+y0+geo->TPCYCent(), z0+geo->TPCZCent());
      ypoint.SetPoint(2, x0+geo->TPCXCent(), 0.95*sz+y0+geo->TPCYCent(), z0+0.05*sz+geo->TPCZCent());

      zpoint.SetPoint(0, x0-0.05*sz+geo->TPCXCent(), y0+geo->TPCYCent(), 0.95*sz+z0+geo->TPCZCent());
      zpoint.SetPoint(1, x0+0.00*sz+geo->TPCXCent(), y0+geo->TPCYCent(), 1.00*sz+z0+geo->TPCZCent());
      zpoint.SetPoint(2, x0+0.05*sz+geo->TPCXCent(), y0+geo->TPCYCent(), 0.95*sz+z0+geo->TPCZCent());

      TPolyLine3D& zleg = view->AddPolyLine3D(4, c, w, s);
      zleg.SetPoint(0,  x0-0.05*sz+geo->TPCXCent(), y0+0.05*sz+geo->TPCYCent(), z0+1.05*sz+geo->TPCZCent());
      zleg.SetPoint(1,  x0+0.05*sz+geo->TPCXCent(), y0+0.05*sz+geo->TPCYCent(), z0+1.05*sz+geo->TPCZCent());
      zleg.SetPoint(2,  x0-0.05*sz+geo->TPCXCent(), y0-0.05*sz+geo->TPCYCent(), z0+1.05*sz+geo->TPCZCent());
      zleg.SetPoint(3,  x0+0.05*sz+geo->TPCXCent(), y0-0.05*sz+geo->TPCYCent(), z0+1.05*sz+geo->TPCZCent());

      TPolyLine3D& yleg = view->AddPolyLine3D(5, c, w, s);
      yleg.SetPoint(0,  x0-0.05*sz+geo->TPCXCent(), y0+1.15*sz+geo->TPCYCent(), z0+geo->TPCZCent());
      yleg.SetPoint(1,  x0+0.00*sz+geo->TPCXCent(), y0+1.10*sz+geo->TPCYCent(), z0+geo->TPCZCent());
      yleg.SetPoint(2,  x0+0.00*sz+geo->TPCXCent(), y0+1.05*sz+geo->TPCYCent(), z0+geo->TPCZCent());
      yleg.SetPoint(3,  x0+0.00*sz+geo->TPCXCent(), y0+1.10*sz+geo->TPCYCent(), z0+geo->TPCZCent());
      yleg.SetPoint(4,  x0+0.05*sz+geo->TPCXCent(), y0+1.15*sz+geo->TPCYCent(), z0+geo->TPCZCent());

      TPolyLine3D& xleg = view->AddPolyLine3D(7, c, w, s);
      xleg.SetPoint(0,  x0+1.05*sz+geo->TPCXCent(), y0+0.05*sz+geo->TPCYCent(), z0-0.05*sz+geo->TPCZCent());
      xleg.SetPoint(1,  x0+1.05*sz+geo->TPCXCent(), y0+0.00*sz+geo->TPCYCent(), z0-0.00*sz+geo->TPCZCent());
      xleg.SetPoint(2,  x0+1.05*sz+geo->TPCXCent(), y0+0.05*sz+geo->TPCYCent(), z0+0.05*sz+geo->TPCZCent());
      xleg.SetPoint(3,  x0+1.05*sz+geo->TPCXCent(), y0+0.00*sz+geo->TPCYCent(), z0-0.00*sz+geo->TPCZCent());
      xleg.SetPoint(4,  x0+1.05*sz+geo->TPCXCent(), y0-0.05*sz+geo->TPCYCent(), z0-0.05*sz+geo->TPCZCent());
      xleg.SetPoint(5,  x0+1.05*sz+geo->TPCXCent(), y0+0.00*sz+geo->TPCYCent(), z0-0.00*sz+geo->TPCZCent());
      xleg.SetPoint(6,  x0+1.05*sz+geo->TPCXCent(), y0-0.05*sz+geo->TPCYCent(), z0+0.05*sz+geo->TPCZCent());

    }

  }
}// namespace
////////////////////////////////////////////////////////////////////////
