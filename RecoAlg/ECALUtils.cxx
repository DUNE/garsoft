////////////////////////////////////////////////////////////////////////
/// \file  ECALUtils.cxx
/// \brief Interface to algorithm class ECAL specific
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/ECALUtils.h"
#include <cmath>

namespace util {

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils()
    : fNeffPx(0)
    {
    }

    //----------------------------------------------------------------------
    ECALUtils::ECALUtils(double NeffPx)
    : fNeffPx(NeffPx)
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

        double ratio = 0.95;
        if(sat_px < ratio * fNeffPx)
        {
            /* desaturate using inverse function */
            float unSaturatedSignal = - fNeffPx * std::log(1 - sat_px / fNeffPx);
            return unSaturatedSignal;
        }
        else
        {
            /* desaturate using linear continuation function for hits above linearisation threshold */
            float unSaturatedSignal = 1/( 1 - ratio ) * (sat_px - ratio * fNeffPx) - fNeffPx * std::log( 1 - ratio );
            return unSaturatedSignal;
        }
    }

} // util
