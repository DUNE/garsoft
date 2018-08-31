//
//  CaloDeposit.cxx
//  garsoft-mrb
//
//  Created by Eldwan Brianne on 8/23/18.
//

#include "SimulationDataProducts/CaloDeposit.h"
#include <limits>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace sdp{

    //-------------------------------------------------
    CaloDeposit::CaloDeposit()
    : fTime   (std::numeric_limits<double>::max())
    , fEnergy (std::numeric_limits<float >::max())
    , fX      (std::numeric_limits<float >::max())
    , fY		  (std::numeric_limits<float >::max())
    , fZ		  (std::numeric_limits<float >::max())
    , fCaloID		(std::numeric_limits<float >::max())
    {}

      //-------------------------------------------------
      // order the hit deposits by time
      bool gar::sdp::CaloDeposit::operator<(gar::sdp::CaloDeposit const& b) const
      {
        if( fTime <  b.Time() )
        return true;

        return false;
      }

    } //sdp
  } // gar
