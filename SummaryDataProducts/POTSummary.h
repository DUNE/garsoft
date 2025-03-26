////////////////////////////////////////////////////////////////////////
/// \file POTSummary.h
///
/// Definition of object to store pot related information
///
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef POTSUM_H
#define POTSUM_H

#include <iostream>

namespace gar {
  namespace sumdata {

    class POTSummary {
    public:
      POTSummary() = default;
      void aggregate(POTSummary const& other);

    public:
      double totpot;
      double totgoodpot;
      unsigned int totspills;
      unsigned int goodspills;

#ifndef __GCCXML__

      friend std::ostream& operator<<(std::ostream& o, POTSummary const& a);

    public:
      double const& TotalPOT() const;
      unsigned int const& TotalSpills() const;
#endif
    };

  } // namespace sumdata
} // namespace gar

#ifndef __GCCXML__

inline double const& gar::sumdata::POTSummary::TotalPOT() const
{
  return totpot;
}
inline unsigned int const& gar::sumdata::POTSummary::TotalSpills() const
{
  return totspills;
}

#endif

#endif //POTSUM_H
