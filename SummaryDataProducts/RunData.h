////////////////////////////////////////////////////////////////////////
/// \file RunData.h
///
/// Definition of object to store run related information
///
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SD_RUNDATA_H
#define SD_RUNDATA_H

#include <string>

namespace gar {
  namespace sumdata {

    class RunData {

    public:
      RunData(); // Default constructor
      void aggregate(RunData const& other);

    private:
      std::string fDetName; ///< detector name
#ifndef __GCCXML__

    public:
      explicit RunData(std::string const& detectorName);
      std::string const& DetName() const;

#endif
    };
  }
} // gar

#ifndef __GCCXML__

inline std::string const& gar::sumdata::RunData::DetName() const
{
  return fDetName;
}

#endif

#endif // SD_RUNDATA_H
