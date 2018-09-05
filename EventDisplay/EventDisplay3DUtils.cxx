//
// - Serves as the link between ROOT "events" (e.g. mouse-clicks) and the ART
//   event display service by providing a receiver slot for signals generated
//   by the ROOT events.  A ROOT dictionary needs to be generated for this.
//

#include "nutools/EventDisplayBase/NavState.h"
#include "EventDisplay3DUtils.h"

#include <string>
#include <cmath>

#include <TColor.h>

namespace gar{

  namespace evd
  {
    EventDisplay3DUtils::EventDisplay3DUtils()
    :fTbRun(0)
    ,fTbEvt(0)
    {
      Double_t s[]  = { 0.00, 0.20, 0.35, 0.43, 0.45, 0.50, 1.00};
      Double_t r[]  = { 0.00, 0.90, 1.00, 1.00, 1.00, 1.00, 0.00};
      Double_t g[]  = { 0.50, 0.10, 0.75, 0.90, 1.00, 1.00, 0.00};
      Double_t b[]  = { 0.30, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00};
      fColorCount = 200;
      unsigned int abc = sizeof(s)/sizeof(s[0]);

      fColorBase = TColor::CreateGradientColorTable(abc, s, r, g, b,
                                                    fColorCount);
    }

      void EventDisplay3DUtils::PrevEvent()
      {
        evdb::NavState::Set(evdb::kPREV_EVENT);
      }

      void EventDisplay3DUtils::NextEvent()
      {
        evdb::NavState::Set(evdb::kNEXT_EVENT);
      }

      void EventDisplay3DUtils::GotoEvent()
      {
        int run = std::stoi(fTbRun->GetString());
        int event = std::stoi(fTbEvt->GetString());
        evdb::NavState::SetTarget(run, event);
        evdb::NavState::Set(evdb::kGOTO_EVENT);
      }

      int EventDisplay3DUtils::LogColor(double value, double minVal, double maxVal, double magScale)
      {
        int nCol = fColorCount/2;
        double scale = std::pow(10.0,magScale);
        double nvalue = std::max(0.0,std::min((value-minVal)/(maxVal-minVal),1.0));
        double lValue = std::log10(1.0+scale*nvalue)/magScale;
        int iValue = nCol*lValue;
        iValue = std::max(0,std::min(iValue,nCol-1));
        return fColorBase + iValue;
      }
    }
  }
