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
        , fCellID (0)
        , fStepLength(0.)
        {
            fPos[0] = 0.;
            fPos[1] = 0.;
            fPos[2] = 0.;
        }

        //-------------------------------------------------
        // order the hit deposits by time
        bool gar::sdp::CaloDeposit::operator<(gar::sdp::CaloDeposit const& b) const
        {
            if( fTime <  b.Time() )
            return true;

            return false;
        }

        void gar::sdp::CaloDeposit::operator+=(gar::sdp::CaloDeposit const& b)
        {
            //Simply adding contributions
            fEnergy += b.Energy();

            //Taking the earliest time
            if(b.Time() < fTime)
            fTime = b.Time();
        }

    } //sdp
} // gar
