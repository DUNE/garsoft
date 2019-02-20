///
/// \file    HeaderDrawer.h
/// \brief   Render the objects from the Simulation package
/// \author  messier@indiana.edu
/// \version $Id: HeaderDrawer.h,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
///
#ifndef EVD_HEADERDRAWER_H
#define EVD_HEADERDRAWER_H
#include <string>
#include <vector>
#include <map>

namespace art  { class Event;  }
namespace evdb { class View2D; class View3D; }

namespace gar {
namespace evd {
  class HeaderDrawer {
  public:
    HeaderDrawer();
    ~HeaderDrawer();
    
    void Header(evdb::View2D* view);
    void Header(evdb::View3D* view);
    
    void Text(std::string& title,
	      std::string& run,
	      std::string& event,
	      std::string& date,
	      std::string& time);
  public:
  };
}
}
#endif
////////////////////////////////////////////////////////////////////////
