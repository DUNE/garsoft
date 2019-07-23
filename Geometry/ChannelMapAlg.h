////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAPALG_H
#define GEO_CHANNELMAPALG_H

// Framework libraries
#include "cetlib_except/exception.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <vector>
#include <map>
#include <set>

#include "Geometry/GeometryCore.h"

namespace gar{
  namespace geo{
    
    /// Exception thrown on invalid wire number (e.g. NearestWireID())
    class InvalidChannelIDError: public cet::exception {
    public:
      InvalidChannelIDError(std::string const& cat): cet::exception(cat) {}
      
      InvalidChannelIDError(std::string const& cat,
                            int                bad_chan,
                            int                better_chan = -1)
      : cet::exception    (cat)
      , chan_number       (bad_chan)
      , better_chan_number(better_chan)
      {}
      
      int chan_number        = std::numeric_limits<int>::min(); ///< the invalid wire number
      int better_chan_number = std::numeric_limits<int>::min(); ///< a suggestion for a good wire number

    }; // class InvalidWireIDError
    
    
    class ChannelMapAlg{
      
    public:
      
      virtual ~ChannelMapAlg() = default;
      
      virtual void          Initialize(gar::geo::GeometryCore & geo)   = 0;
      virtual void          Uninitialize()                             = 0;
      virtual unsigned int  Nchannels()                          const = 0;
      virtual unsigned int  NearestChannel(float const* xyz)     const = 0;
      virtual unsigned int  GapChannelNumber()                   const = 0;
      virtual void          ChannelToPosition(unsigned int chan,
                                              float*       xyz)  const = 0;

    protected:
      
    };
  }
} // namespace gar

#endif // GEO_CHANNELMAPALG_H

