///
/// \file    Display3DView.h
/// \brief   A view showing a 3D rendering of the detector
/// \author  messier@indiana.edu
/// \version $Id: Display3DView.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
#ifndef EVD_DISPLAY3DVIEW_H
#define EVD_DISPLAY3DVIEW_H
#include "RQ_OBJECT.h"

#include "TGButton.h" // For TGCheckButton

#include "nutools/EventDisplayBase/Canvas.h"

#include "EventDisplay/HeaderPad.h"
#include "EventDisplay/MCBriefPad.h"

namespace gar {
namespace evd {
  class Display3DPad;

  /// View of event shoing the XZ and YZ readout planes
  class Display3DView : public evdb::Canvas {

  public:
    
    RQ_OBJECT("gar::evd::Display3DView")

  public:
    Display3DView(TGMainFrame* mf);
    ~Display3DView();
    
    const char* Description() const { return "3D Detector Display"; }
    const char* PrintTag()    const { return "duneNDMPD3d";               }
    void Draw(const char* opt="");
    void CloseWindow();
    void SetRawReco();
    void SetMCInfo();
    
  private:
    
    Display3DPad*  fDisplay3DPad; /// Pad showing 3D view of the detector
    HeaderPad*     fHeaderPad;    ///< Show header information
    MCBriefPad*    fMC;           ///< Short summary of MC event
    TGCheckButton* fMCOn;         ///< Display MC truth information
    TGRadioButton* fRawDraw;      ///< Draw Raw information only
    TGRadioButton* fRecoDraw;     ///< Draw Reconstructed information only
    
  };
}
}
#endif
////////////////////////////////////////////////////////////////////////
