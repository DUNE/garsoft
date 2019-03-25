////////////////////////////////////////////////////////////////////////
/// \file  ECALUtils.h
/// \brief Interface to ECALReadoutSimAlg class for ECAL specific
///
/// \version $Id:  $
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////
#ifndef ECALUTILS_H
#define ECALUTILS_H

namespace util {

    class ECALUtils {

    public:

        ECALUtils();

        ECALUtils(double NeffPx);

        ~ECALUtils();

        double Saturate(const double unsat_px);

        double DeSaturate(const double sat_px);

    private:
        double fNeffPx;
    };
} // util

#endif // ECALUTILS_H
