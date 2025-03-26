//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#ifndef EventDisplay3D_EvtDisplayUtils_hh
#define EventDisplay3D_EvtDisplayUtils_hh

#include <TApplication.h>
#include <TGTextBuffer.h>
#include <TObject.h>
#include <iostream>

namespace gar {
  namespace evd3d {
    class EventDisplay3DUtils {

    public:
      explicit EventDisplay3DUtils();
      void PrevEvent();
      void NextEvent();
      void GotoEvent();
      TGTextBuffer* fTbRun;
      TGTextBuffer* fTbEvt;

      int LogColor(double val, double minVal, double maxVal, double magScale = 5.0);

      int fColorBase;
      int fColorCount;
    };
  }
}
#endif /* EventDisplay3D_EvtDisplayUtils_hh */
