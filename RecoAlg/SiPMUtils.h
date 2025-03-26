////////////////////////////////////////////////////////////////////////
/// \file  SiPMUtils.h
/// \brief Interface to SiPMReadoutSimAlg class for SiPM specific
///
/// \version $Id:  $
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////
#ifndef SiPMUtils_H
#define SiPMUtils_H

namespace util {

  class SiPMUtils {

  public:
    SiPMUtils();

    SiPMUtils(double NeffPx);

    ~SiPMUtils();

    double Saturate(const double unsat_px);

    double DeSaturate(const double sat_px);

  private:
    double fNeffPx;
  };
} // util

#endif // SiPMUtils_H
