/** ****************************************************************************
* @file CaloRawDigit.cxx
* @brief Definition of basic calo raw digits
* @author eldwan.brianne@desy.de
* @see  CaloRawDigit.h raw.h
* Created on the 08/27/18
*
* ****************************************************************************/

#include <limits>

#include "RawDataProducts/CaloRawDigit.h"

// C/C++ standard libraries

namespace gar {
  namespace raw{

    //----------------------------------------------------------------------
    CaloRawDigit::CaloRawDigit()
    : fADC(0),
    fTime(0),
    fX(0),
    fY(0),
    fZ(0),
    fCaloID(-1)
    {}

    //----------------------------------------------------------------------
    CaloRawDigit::CaloRawDigit(double ADC, double time, double x, double y, double z, int id)
    : fADC(ADC),
    fTime(time),
    fX(x),
    fY(y),
    fZ(z),
    fCaloID(id)
    {}

    } // namespace raw
  } // gar
    ////////////////////////////////////////////////////////////////////////
