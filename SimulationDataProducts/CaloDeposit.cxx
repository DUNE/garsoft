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
    : fTrackID(0)
    , fTime   (0.)
    , fEnergy (0.)
    , fX      (0.)
    , fY		  (0.)
    , fZ		  (0.)
    , fCaloID	(0)
    , fCellID (0)
    , fLayer(0)
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
