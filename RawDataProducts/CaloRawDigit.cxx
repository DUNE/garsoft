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
    fTime(std::make_pair(0., 0.)),
    fX(0.),
    fY(0.),
    fZ(0.),
    fCellID(0)
    {}

    //----------------------------------------------------------------------
    CaloRawDigit::CaloRawDigit(unsigned int ADC, float time, float x, float y, float z, long long int cellID)
    : fADC(ADC),
    fX(x),
    fY(y),
    fZ(z),
    fCellID(cellID)
    {
        fTime = std::make_pair(time, 0.);
    }

    //----------------------------------------------------------------------
    CaloRawDigit::CaloRawDigit(unsigned int ADC, std::pair<float, float> time, float x, float y, float z, long long int cellID)
    : fADC(ADC),
    fTime(time),
    fX(x),
    fY(y),
    fZ(z),
    fCellID(cellID)
    {}

    } // namespace raw
  } // gar
    ////////////////////////////////////////////////////////////////////////
