///
/// \file    CalorView.h
/// \brief   A view showing calorimetric particle ID information.
/// \author  msoderbe@syr.edu
///

#ifndef EVD_CALORVIEW_H
#define EVD_CALORVIEW_H

#include "RQ_OBJECT.h"

#include "nuevdb/EventDisplayBase/Canvas.h"

namespace gar {
  namespace evd {

    class CalorPad;

    /// View showing calorimetric particle ID information

    class CalorView : public evdb::Canvas {

    public:
      RQ_OBJECT("gar::evd::CalorView")

    public:
      // Constructor, destructor.
      CalorView(TGMainFrame* mf);
      virtual ~CalorView();

      // Required methods.
      const char* Description() const { return "Calorimetric PID Display"; }
      const char* PrintTag() const { return "larcalor"; }
      void Draw(const char* opt = "");
      void CloseWindow();

    private:
      CalorPad* fDeDxPad; ///< Graphics pad for dEdx vs. Res. range
      CalorPad* fKEPad;   ///< Graphics pad for KE vs. Total range
    };
  }
}
#endif
////////////////////////////////////////////////////////////////////////
