//
//  EnergyDeposit.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 1/30/17.
//

#include "SimulationDataProducts/LArDeposit.h"
#include <limits>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace sdp{

    //-------------------------------------------------
    LArDeposit::LArDeposit()
    : fTrackID (0)
    , fTime   (0.)
    , fEnergy (0.)
    , fX      (0.)
    , fY		  (0.)
    , fZ		  (0.)
    , fdX		  (0.)
    , fIsPrimary(true)
    {}

    //-------------------------------------------------
    // order the energy deposits by id, then time, then z position
    bool gar::sdp::LArDeposit::operator<(gar::sdp::LArDeposit const& b) const
    {
      if( std::abs(fTrackID) <  std::abs(b.TrackID()) )
        return true;
      else if( std::abs(fTrackID) == std::abs(b.TrackID()) ){
        if( fTime < b.Time() )
          return true;
        else if(fTime == b.Time() ){
          if(fZ < b.Z()) return true;
        }
      }

      return false;
    }

  } //sdp
} // gar
