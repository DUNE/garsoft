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
    : fTime   (1.0E7)
    , fEnergy (1.0E7)
    , fX      (1.0E7)
    , fY		  (1.0E7)
    , fZ		  (1.0E7)
    , fCaloID		(1.0E7)
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
