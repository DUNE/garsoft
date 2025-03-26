///
/// \file    HeaderPad.h
/// \brief   Drawing pad for time or charge histograms
/// \author  messier@indiana.edu
/// \version $Id: HeaderPad.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
///
#ifndef EVD_HEADER_H
#define EVD_HEADER_H
#include "EventDisplay/EVD/DrawingPad.h"
namespace evdb {
  class View2D;
}
class TH1F;

namespace gar {
  namespace evd {
    class HeaderPad : public DrawingPad {
    public:
      HeaderPad(const char* nm,
                const char* ti,
                double x1,
                double y1,
                double x2,
                double y2,
                const char* opt);
      ~HeaderPad();
      void Draw(const char* opt = "");

    private:
      evdb::View2D* fView; ///< Collection of drawn objects
    };
  }
}
#endif
////////////////////////////////////////////////////////////////////////
