///
/// \file    HeaderDrawer.cxx
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: HeaderDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#include "EventDisplay/HeaderDrawer.h"
#include "nutools/EventDisplayBase/Colors.h"
#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"

#include "TText.h"
#include "TTimeStamp.h"

#include "art/Framework/Principal/Event.h"

namespace gar {
namespace evd{

  HeaderDrawer::HeaderDrawer() { }

  //......................................................................

  HeaderDrawer::~HeaderDrawer() { }


  //......................................................................

  void HeaderDrawer::Text(std::string& titles,
			  std::string& runs,
			  std::string& events,
			  std::string& dates,
			  std::string& times)
  {
    titles = "GOAT";
    
    // get the event
    const art::Event* evt = evdb::EventHolder::Instance()->GetEvent();
    if(!evt) return;

    int run   = evt->run();
    int srun  = evt->subRun();
    int event = evt->id().event();

    unsigned int year, month, day, dayofweek;
    unsigned int hour, minute, second;
    int          nano;

    // get the time stamp.  art::Timestamp::value() returns a TimeValue_t which is a typedef to unsigned long long.
    // The conventional use is for the upper 32 bits to have the seconds since 1970 epoch and the lower 32 bits to be
    // the number of microseconds with the current second.
    unsigned long long int tsval = evt->time().value();
  
    // taking it apart
    // the masking isn't strictly necessary *if* "long" is truly 32bits
    // but this can vary w/ compiler/platform
    const unsigned long int mask32 = 0xFFFFFFFFUL;
    unsigned long int lup = ( tsval >> 32 ) & mask32;
    unsigned long int llo = tsval & mask32;
    TTimeStamp ts(lup, (int)llo);
    
    ts.GetDate(kTRUE,0,&year,&month,&day);
    ts.GetTime(kTRUE,0,&hour,&minute,&second);
    nano = ts.GetNanoSec();
    dayofweek = ts.GetDayOfWeek();
    char eventbuff[256];
    char runbuff[256];
    char datebuff[256];
    char timebuff[256];
  
    // Skip first one since ROOT returns these numbers starting from 1 not 0
    static const char* days[] = {"",
				 "Mon","Tue","Wed","Thu","Fri","Sat","Sun"
    };
    static const char* months[] = {"",
				   "Jan","Feb","Mar","Apr","May","Jun",
				   "Jul","Aug","Sep","Oct","Nov","Dec"
    };
  
    sprintf(runbuff,  "Run:   %d/%d",run,srun);
    runs = runbuff;
    
    sprintf(eventbuff,"Event: %d",event);
    events = eventbuff;
    
    sprintf(datebuff, "UTC %s %s %d, %d", 
            days[dayofweek],
            months[month],
            day,
            year);
    dates = datebuff;
    
    sprintf(timebuff, "%.2d:%.2d:%2.9f", 
            hour,
            minute,
            (float)second+(float)nano/1.0E9);
    times = timebuff;
  }

  void HeaderDrawer::Header(evdb::View2D* view)
  {
    std::string title;
    std::string run;
    std::string event;
    std::string date;
    std::string time;
    this->Text(title,run,event,date,time);

    art::ServiceHandle<evdb::Colors> colors;
    int c = colors->Foreground(0);
    
    TText& titlet = view->AddText(0.03,0.80, title.c_str());
    TText& runt   = view->AddText(0.04,0.60, run.c_str());
    TText& eventt = view->AddText(0.04,0.45, event.c_str());
    TText& datet  = view->AddText(0.04,0.25, date.c_str());
    TText& timet  = view->AddText(0.04,0.10, time.c_str());

    titlet.SetTextSize(0.13);
    titlet.SetTextFont(62);
    titlet.SetTextColor(c);
    
    runt.SetTextSize(0.12);
    runt.SetTextFont(42);
    runt.SetTextColor(c);
    
    eventt.SetTextSize(0.12);
    eventt.SetTextFont(42);
    eventt.SetTextColor(c);
    
    datet.SetTextSize(0.12);
    datet.SetTextFont(42);
    datet.SetTextColor(c);
    
    timet.SetTextSize(0.12);
    timet.SetTextFont(42);
    timet.SetTextColor(c);
  }
  void HeaderDrawer::Header(evdb::View3D* view)
  {
    std::string title;
    std::string run;
    std::string event;
    std::string date;
    std::string time;
    this->Text(title,run,event,date,time);

    art::ServiceHandle<evdb::Colors> colors;
    int c = colors->Foreground(0);
    
    TText& titlet = view->AddText(-0.98,-0.75, title.c_str());
    TText& runt   = view->AddText(-0.98,-0.80, run.c_str());
    TText& eventt = view->AddText(-0.98,-0.85, event.c_str());
    TText& datet  = view->AddText(-0.98,-0.90, date.c_str());
    TText& timet  = view->AddText(-0.98,-0.95, time.c_str());

    titlet.SetTextSize(0.05);
    titlet.SetTextFont(62);
    titlet.SetTextColor(c);
    
    runt.SetTextSize(0.04);
    runt.SetTextFont(42);
    runt.SetTextColor(c);

    eventt.SetTextSize(0.04);
    eventt.SetTextFont(42);
    eventt.SetTextColor(c);
    
    datet.SetTextSize(0.04);
    datet.SetTextFont(42);
    datet.SetTextColor(c);
    
    timet.SetTextSize(0.04);
    timet.SetTextFont(42);
    timet.SetTextColor(c);
  }
}
}// namespace
////////////////////////////////////////////////////////////////////////
