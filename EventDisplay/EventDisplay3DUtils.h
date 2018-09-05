//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#ifndef EventDisplay3D_EvtDisplayUtils_hh
#define EventDisplay3D_EvtDisplayUtils_hh

#include <TObject.h>
#include <TApplication.h>
#include <TGTextBuffer.h>
#include <iostream>

namespace gar{
  namespace evd
  {
    class EventDisplay3DUtils
    {

    public:
      explicit EventDisplay3DUtils();
      void PrevEvent();
      void NextEvent();
      void GotoEvent();
      TGTextBuffer *fTbRun;
      TGTextBuffer *fTbEvt;

      int LogColor(double val, double minVal, double maxVal, double magScale=5.0);

      int fColorBase;
      int fColorCount;

    };
  }
}
#endif /* EventDisplay3D_EvtDisplayUtils_hh */
