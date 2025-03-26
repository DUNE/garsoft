/// \file    GeometryDrawer.cxx
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
/// \version $Id: GeometryDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
#include "TBox.h"
#include "TLine.h"
#include "TMath.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TText.h"

#include "EventDisplay/EVD/GeometryDrawer.h"
#include "Geometry/GeometryGAr.h"
#include "nuevdb/EventDisplayBase/View2D.h"
#include "nuevdb/EventDisplayBase/View3D.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

namespace gar {
  namespace evd {

    //......................................................................
    GeometryDrawer::GeometryDrawer() {}

    //......................................................................
    GeometryDrawer::~GeometryDrawer() {}

    //......................................................................
    void GeometryDrawer::DetOutline3D(evdb::View3D* view)
    {
      art::ServiceHandle<geo::GeometryGAr> geo;

      // these are relative to the center of the TPC.  We will add the detector offsets
      // on each point

      double xlo = -geo->TPCLength() / 2;
      double xhi = geo->TPCLength() / 2;
      double ylo = -geo->TPCRadius();
      double yhi = geo->TPCRadius();
      //double zlo =  -geo->TPCRadius();
      double zhi = geo->TPCRadius();
      double r = yhi;

      double fracinner = 788.0 / 2580.;

      double ang = TMath::Pi() * 2.0 / 18.0;

      int c = kGray;
      int s = 1;
      int w = 1;
      TPolyLine3D& spos = view->AddPolyLine3D(19, c, w, s);
      TPolyLine3D& sposi = view->AddPolyLine3D(19, c, w, s);
      for (int i = 0; i < 19; ++i) {
        // nb. coordinates are so Y is up, but ROOT thinks Z is up, so report (z,x,y)
        spos.SetPoint(i,
                      r * TMath::Cos(i * ang) + geo->TPCZCent(),
                      xhi + geo->TPCXCent(),
                      r * TMath::Sin(i * ang) + geo->TPCYCent());
        sposi.SetPoint(i,
                       fracinner * r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xhi + geo->TPCXCent(),
                       fracinner * r * TMath::Sin(i * ang) + geo->TPCYCent());
      }
      TPolyLine3D& sneg = view->AddPolyLine3D(19, c, w, s);
      TPolyLine3D& snegi = view->AddPolyLine3D(19, c, w, s);
      for (int i = 0; i < 19; ++i) {
        // nb. coordinates are so Y is up, but ROOT thinks Z is up, so report (z,x,y)
        sneg.SetPoint(i,
                      r * TMath::Cos(i * ang) + geo->TPCZCent(),
                      xlo + geo->TPCXCent(),
                      r * TMath::Sin(i * ang) + geo->TPCYCent());
        snegi.SetPoint(i,
                       fracinner * r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xlo + geo->TPCXCent(),
                       fracinner * r * TMath::Sin(i * ang) + geo->TPCYCent());
      }

      c = kGray + 2;
      s = 1;
      w = 1;
      for (int i = 0; i < 18; ++i) {
        TPolyLine3D& gridt = view->AddPolyLine3D(4, c, w, s);
        // nb. coordinates are so Y is up, but ROOT thinks Z is up, so report (z,x,y)
        gridt.SetPoint(0,
                       fracinner * r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xlo + geo->TPCXCent(),
                       fracinner * r * TMath::Sin(i * ang) + geo->TPCYCent());
        gridt.SetPoint(1,
                       r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xlo + geo->TPCXCent(),
                       r * TMath::Sin(i * ang) + geo->TPCYCent());
        gridt.SetPoint(2,
                       r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xhi + geo->TPCXCent(),
                       r * TMath::Sin(i * ang) + geo->TPCYCent());
        gridt.SetPoint(3,
                       fracinner * r * TMath::Cos(i * ang) + geo->TPCZCent(),
                       xhi + geo->TPCXCent(),
                       fracinner * r * TMath::Sin(i * ang) + geo->TPCYCent());
      }

      // Indicate coordinate system
      double x0 = -0.20;       // Center location of the key
      double y0 = 1.10 * ylo;  // Center location of the key
      double z0 = -0.10 * zhi; // Center location of the key
      double sz = 0.20 * zhi;  // Scale size of the key in z direction

      c = kCyan;

      TPolyLine3D& xaxis = view->AddPolyLine3D(2, c, w, s);
      TPolyLine3D& yaxis = view->AddPolyLine3D(2, c, w, s);
      TPolyLine3D& zaxis = view->AddPolyLine3D(2, c, s, w);
      // nb. coordinates are so Y is up, but ROOT thinks Z is up, so report (z,x,y)
      xaxis.SetPoint(0, z0 + geo->TPCZCent(), x0 + geo->TPCXCent(), y0 + geo->TPCYCent());
      xaxis.SetPoint(1, z0 + geo->TPCZCent(), sz + x0 + geo->TPCXCent(), y0 + geo->TPCYCent());

      yaxis.SetPoint(0, z0 + geo->TPCZCent(), x0 + geo->TPCXCent(), y0 + geo->TPCYCent());
      yaxis.SetPoint(1, z0 + geo->TPCZCent(), x0 + geo->TPCXCent(), y0 + sz + geo->TPCYCent());

      zaxis.SetPoint(0, z0 + geo->TPCZCent(), x0 + geo->TPCXCent(), y0 + geo->TPCYCent());
      zaxis.SetPoint(1, z0 + sz + geo->TPCZCent(), x0 + geo->TPCXCent(), y0 + geo->TPCYCent());

      TPolyLine3D& xpoint = view->AddPolyLine3D(3, c, w, s);
      TPolyLine3D& ypoint = view->AddPolyLine3D(3, c, w, s);
      TPolyLine3D& zpoint = view->AddPolyLine3D(3, c, w, s);

      xpoint.SetPoint(0,
                      z0 - 0.05 * sz + geo->TPCZCent(),
                      0.95 * sz + x0 + geo->TPCXCent(),
                      y0 + geo->TPCYCent());
      xpoint.SetPoint(
        1, z0 + geo->TPCZCent(), 1.00 * sz + x0 + geo->TPCXCent(), y0 + geo->TPCYCent());
      xpoint.SetPoint(2,
                      z0 + 0.05 * sz + geo->TPCZCent(),
                      0.95 * sz + x0 + geo->TPCXCent(),
                      y0 + geo->TPCYCent());

      ypoint.SetPoint(0,
                      z0 - 0.05 * sz + geo->TPCZCent(),
                      x0 + geo->TPCXCent(),
                      0.95 * sz + y0 + geo->TPCYCent());
      ypoint.SetPoint(
        1, z0 + geo->TPCZCent(), x0 + geo->TPCXCent(), 1.00 * sz + y0 + geo->TPCYCent());
      ypoint.SetPoint(2,
                      z0 + 0.05 * sz + geo->TPCZCent(),
                      x0 + geo->TPCXCent(),
                      0.95 * sz + y0 + geo->TPCYCent());

      zpoint.SetPoint(0,
                      0.95 * sz + z0 + geo->TPCZCent(),
                      x0 - 0.05 * sz + geo->TPCXCent(),
                      y0 + geo->TPCYCent());
      zpoint.SetPoint(1,
                      1.00 * sz + z0 + geo->TPCZCent(),
                      x0 + 0.00 * sz + geo->TPCXCent(),
                      y0 + geo->TPCYCent());
      zpoint.SetPoint(2,
                      0.95 * sz + z0 + geo->TPCZCent(),
                      x0 + 0.05 * sz + geo->TPCXCent(),
                      y0 + geo->TPCYCent());

      TPolyLine3D& zleg = view->AddPolyLine3D(4, c, w, s);
      zleg.SetPoint(0,
                    z0 + 1.05 * sz + geo->TPCZCent(),
                    x0 - 0.05 * sz + geo->TPCXCent(),
                    y0 + 0.05 * sz + geo->TPCYCent());
      zleg.SetPoint(1,
                    z0 + 1.05 * sz + geo->TPCZCent(),
                    x0 + 0.05 * sz + geo->TPCXCent(),
                    y0 + 0.05 * sz + geo->TPCYCent());
      zleg.SetPoint(2,
                    z0 + 1.05 * sz + geo->TPCZCent(),
                    x0 - 0.05 * sz + geo->TPCXCent(),
                    y0 - 0.05 * sz + geo->TPCYCent());
      zleg.SetPoint(3,
                    z0 + 1.05 * sz + geo->TPCZCent(),
                    x0 + 0.05 * sz + geo->TPCXCent(),
                    y0 - 0.05 * sz + geo->TPCYCent());

      TPolyLine3D& yleg = view->AddPolyLine3D(5, c, w, s);
      yleg.SetPoint(0,
                    z0 + geo->TPCZCent(),
                    x0 - 0.05 * sz + geo->TPCXCent(),
                    y0 + 1.15 * sz + geo->TPCYCent());
      yleg.SetPoint(1,
                    z0 + geo->TPCZCent(),
                    x0 + 0.00 * sz + geo->TPCXCent(),
                    y0 + 1.10 * sz + geo->TPCYCent());
      yleg.SetPoint(2,
                    z0 + geo->TPCZCent(),
                    x0 + 0.00 * sz + geo->TPCXCent(),
                    y0 + 1.05 * sz + geo->TPCYCent());
      yleg.SetPoint(3,
                    z0 + geo->TPCZCent(),
                    x0 + 0.00 * sz + geo->TPCXCent(),
                    y0 + 1.10 * sz + geo->TPCYCent());
      yleg.SetPoint(4,
                    z0 + geo->TPCZCent(),
                    x0 + 0.05 * sz + geo->TPCXCent(),
                    y0 + 1.15 * sz + geo->TPCYCent());

      TPolyLine3D& xleg = view->AddPolyLine3D(7, c, w, s);
      xleg.SetPoint(0,
                    z0 - 0.05 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 + 0.05 * sz + geo->TPCYCent());
      xleg.SetPoint(1,
                    z0 - 0.00 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 + 0.00 * sz + geo->TPCYCent());
      xleg.SetPoint(2,
                    z0 + 0.05 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 + 0.05 * sz + geo->TPCYCent());
      xleg.SetPoint(3,
                    z0 - 0.00 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 + 0.00 * sz + geo->TPCYCent());
      xleg.SetPoint(4,
                    z0 - 0.05 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 - 0.05 * sz + geo->TPCYCent());
      xleg.SetPoint(5,
                    z0 - 0.00 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 + 0.00 * sz + geo->TPCYCent());
      xleg.SetPoint(6,
                    z0 + 0.05 * sz + geo->TPCZCent(),
                    x0 + 1.05 * sz + geo->TPCXCent(),
                    y0 - 0.05 * sz + geo->TPCYCent());
    }

  }
} // namespace
////////////////////////////////////////////////////////////////////////
