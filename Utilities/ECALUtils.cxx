////////////////////////////////////////////////////////////////////////
/// \file  ECALUtils.cxx
/// \brief Interface to algorithm class ECAL specific
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#define SHIFT_CALOID 61
#define SHIFT_NEG_I 60
#define SHIFT_I 34
#define SHIFT_NEG_J 33
#define SHIFT_J 7
#define SHIFT_LAYER 0

#include "Utilities/ECALUtils.h"

#include "SimulationDataProducts/CaloDeposit.h"

#include <cmath>

namespace util {

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils()
    : fNeffPx(0),
    fRatio(0)
    {
    }

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils(double NeffPx, double ratio)
    : fNeffPx(NeffPx),
    fRatio(ratio)
    {
    }

    //----------------------------------------------------------------------
    ECALUtils::~ECALUtils()
    {
    }

    //----------------------------------------------------------------------------
    double ECALUtils::Saturate(const double unsat_px)
    {
        /*  protect against negative/zero input*/
        if ( unsat_px <= 0 ) return unsat_px;

        /* calculate saturated signal using simple exponential formula */
        double saturatedSignal = fNeffPx * ( 1. - std::exp( - unsat_px / fNeffPx ) );

        return saturatedSignal;
    }

    //----------------------------------------------------------------------------
    double ECALUtils::DeSaturate(const double sat_px)
    {
        /*  protect against negative/zero input*/
        if ( sat_px <= 0 ) return sat_px;

        if(sat_px < fRatio * fNeffPx)
        {
            /* desaturate using inverse function */
            float unSaturatedSignal = - fNeffPx * std::log(1 - sat_px / fNeffPx);
            return unSaturatedSignal;
        }
        else
        {
            /* desaturate using linear continuation function for hits above linearisation threshold */
            float unSaturatedSignal = 1/( 1 - fRatio ) * (sat_px - fRatio * fNeffPx) - fNeffPx * std::log( 1 - fRatio );
            return unSaturatedSignal;
        }
    }

} // util
