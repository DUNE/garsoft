///
/// \file    Display3DPad.h
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
/// \version $Id: Display3DPad.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
#ifndef EVD_DISPLAY3DPAD_H
#define EVD_DISPLAY3DPAD_H

#include "Buttons.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGDimension.h"
#include "TGFrame.h" // For TGMainFrame, TGHorizontalFrame
#include "TGLabel.h"
#include "TGLayout.h" // For TGLayoutHints
#include "TGNumberEntry.h"
#include "TGTextView.h"
#include "TLine.h"
#include "TMath.h"
#include "TROOT.h"
#include "TRootEmbeddedCanvas.h"
#include "TString.h"

#include "EventDisplay/EVD/DrawingPad.h"
#include "EventDisplay/EVD/GeometryDrawer.h"
#include "EventDisplay/EVD/HeaderPad.h"
#include "EventDisplay/EVD/MCBriefPad.h"

class TH3F;

namespace gar {
  namespace evd {
    class RawDataDrawer;
    class RecoBaseDrawer;

    /// A drawing pad showing a 3D rendering of the detector
    class Display3DPad : public DrawingPad {
    public:
      Display3DPad(const char* nm,
                   const char* ti,
                   double x1,
                   double y1,
                   double x2,
                   double y2,
                   const char* opt);
      ~Display3DPad();

      void Draw();

    private:
      evdb::View3D* fView; ///< Collection of graphics objects to render
    };
  }
}
#endif
////////////////////////////////////////////////////////////////////////
