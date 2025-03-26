///
/// \file    HeaderPad.cxx
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: HeaderPad.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "EventDisplay/EVD/HeaderPad.h"
#include "EventDisplay/EVD/HeaderDrawer.h"
#include "nuevdb/EventDisplayBase/View2D.h"

#include "TPad.h"
#include "TText.h"

namespace gar {
  namespace evd {

    // clang didn't like these unused varaibles

    //static const int kRAW   = 0;
    //static const int kCALIB = 1;
    //static const int kPE =  2;
    //static const int kTNS = 3;

    //......................................................................

    HeaderPad::HeaderPad(const char* nm,
                         const char* ti,
                         double x1,
                         double y1,
                         double x2,
                         double y2,
                         const char* /*opt*/)
      : DrawingPad(nm, ti, x1, y1, x2, y2)
    {
      fView = new evdb::View2D();
    }

    //......................................................................

    HeaderPad::~HeaderPad()
    {
      if (fView != nullptr) {
        delete fView;
        fView = nullptr;
      }
    }

    //......................................................................

    void HeaderPad::Draw(const char* /* opt */)
    {
      fView->Clear();

      this->HeaderDraw()->Header(fView);

      this->Pad()->Clear();
      this->Pad()->cd();
      fView->Draw();
    }

    //......................................................................

  }
} // namespace
//////////////////////////////////////////////////////////////////////////
