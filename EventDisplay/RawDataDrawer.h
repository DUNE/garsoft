/// \file    RawDataDrawer.h
/// \brief   Class to aid in the rendering of RawData objects
/// \author  messier@indiana.edu
/// \version $Id: RawDataDrawer.h,v 1.2 2010/11/10 22:38:34 p-novaart Exp $
#ifndef EVD_RAWDATADRAWER_H
#define EVD_RAWDATADRAWER_H

#include <vector>
#ifndef __CINT__

#endif

#include  "TH1F.h"

class TVirtualPad;
namespace art  { class Event;  }
namespace evdb { class View3D; }

namespace gar {

  namespace raw { class RawDigit;  }
  namespace raw { class CaloRawDigit;  }

  namespace evd {

  /// Aid in the rendering of RawData objects
  class RawDataDrawer {
  public:

    RawDataDrawer();
    ~RawDataDrawer();

    /**
     * @brief Draws raw digit content in 2D wire plane representation
     * @param evt source for raw digits
     * @param view target rendered object
     * @param plane number of the plane to be drawn
     * @param bZoomToRoI whether to render only te region of interest
     *
     * This function performs pre-rendering of the raw digit content into a
     * 3D view of the yzt space
     * The material for rendering is created and sent to view object for actual
     * rendering.
     *
     */
    void RawDigit3D(art::Event const& evt,
                    evdb::View3D* view);

  private:

    // Fill a histogram with the charge as a function of time
    void FillQHisto(gar::raw::RawDigit const& dig,
                    TH1F*                     histo);

    /// Performs the 3D drawing
    void DrawRawDigit3D (art::Event   const& evt,
                         evdb::View3D*       view);

    /**
     * @brief Makes sure raw::RawDigit's are available for the current settings
     *
     */
    void GetRawDigits(art::Event                        const& evt,
                      std::vector<const raw::RawDigit*>      & digits);

    TH1F*  fEventTQHist; ///< Charge vs time for all digits in an event
    TH1F*  fDigitTQHist; ///< Charge vs time for a single channel in an event

  }; // class RawDataDrawer

}
}

#endif
////////////////////////////////////////////////////////////////////////
