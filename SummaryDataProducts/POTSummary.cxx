////////////////////////////////////////////////////////////////////////
/// \file POTSummary.cxx
///
/// Definition of object to store pot related information
///
/// \version $Id: RunData.h,v 1.1.1.1 2011/03/03 00:19:49 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SummaryDataProducts/POTSummary.h"

namespace gar {
  namespace sumdata {

    //----------------------------------------------------------------------------
    void POTSummary::aggregate(POTSummary const& other)
    {
      totpot += other.totpot;
      totgoodpot += other.totgoodpot;
      totspills += other.totspills;
      goodspills += other.goodspills;
    } // POTSummary::aggregate()

    //----------------------------------------------------------------------------
    std::ostream& operator<<(std::ostream& o, POTSummary const& a)
    {
      o << "This sub run has " << a.totspills << " total spills with an exposure of " << a.totpot
        << " POT" << std::endl
        << " with cuts on beam quality, there are " << a.goodspills
        << " good spills with an exposure of " << a.totgoodpot << std::endl;

      return o;
    }

  } // sumdata
} // end namespace
