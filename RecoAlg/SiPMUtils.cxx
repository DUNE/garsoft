////////////////////////////////////////////////////////////////////////
/// \file  SiPMUtils.cxx
/// \brief Interface to algorithm class SiPM specific
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#include "RecoAlg/SiPMUtils.h"
#include <cmath>

namespace util {

  //----------------------------------------------------------------------
  SiPMUtils::SiPMUtils() : fNeffPx(0) {}

  //----------------------------------------------------------------------
  SiPMUtils::SiPMUtils(double NeffPx) : fNeffPx(NeffPx) {}

  //----------------------------------------------------------------------
  SiPMUtils::~SiPMUtils() {}

  //----------------------------------------------------------------------------
  double SiPMUtils::Saturate(const double unsat_px)
  {
    /*  protect against negative/zero input*/
    if (unsat_px <= 0) return unsat_px;

    /* calculate saturated signal using simple exponential formula */
    double saturatedSignal = fNeffPx * (1. - std::exp(-unsat_px / fNeffPx));

    return saturatedSignal;
  }

  //----------------------------------------------------------------------------
  double SiPMUtils::DeSaturate(const double sat_px)
  {
    /*  protect against negative/zero input*/
    if (sat_px <= 0) return sat_px;

    double ratio = 0.95;
    if (sat_px < ratio * fNeffPx) {
      /* desaturate using inverse function */
      float unSaturatedSignal = -fNeffPx * std::log(1 - sat_px / fNeffPx);
      return unSaturatedSignal;
    }
    else {
      /* desaturate using linear continuation function for hits above linearisation threshold */
      float unSaturatedSignal =
        1 / (1 - ratio) * (sat_px - ratio * fNeffPx) - fNeffPx * std::log(1 - ratio);
      return unSaturatedSignal;
    }
  }

} // util
