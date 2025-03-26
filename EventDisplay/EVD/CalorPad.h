///
/// \file    CalorPad.h
/// \brief   Drawing pad showing calorimetric particle ID information
/// \author  msoderbe@syr.edu
///
#ifndef EVD_CALORPAD_H
#define EVD_CALORPAD_H

#include "EventDisplay/EVD/DrawingPad.h"
#include <string>

class TGraph;

namespace evdb {
  class View2D;
}

namespace gar {
  namespace evd {

    /// A drawing pad showing calorimetric particle ID information

    class CalorPad : public DrawingPad {
    public:
      // Constructor, destructor.
      CalorPad(const char* name,
               const char* title,
               double x1,
               double y1,
               double x2,
               double y2,
               int curvetype);
      ~CalorPad();

      // Methods.
      void Draw(const char* opt = 0);
      void DrawRefCurves();

    private:
      std::string fROOTfile;
      TGraph* dedx_range_pro; ///< proton template
      TGraph* dedx_range_ka;  ///< kaon template
      TGraph* dedx_range_pi;  ///< pion template
      TGraph* dedx_range_mu;  ///< muon template

      TGraph* ke_range_pro; ///< proton template
      TGraph* ke_range_ka;  ///< kaon template
      TGraph* ke_range_pi;  ///< pion template
      TGraph* ke_range_mu;  ///< muon template
      int fcurvetype;       //dEdx vs. Res. range, or Kinetic Energy vs. range

      evdb::View2D* fView; ///< Collection of graphics objects to render; text labels
    };
  }
}
#endif
////////////////////////////////////////////////////////////////////////
