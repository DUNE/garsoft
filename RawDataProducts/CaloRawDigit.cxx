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
    : fCellID(0),
    fADC(0),
    fTime(0),
    fX(0),
    fY(0),
    fZ(0),
    fLayer(0),
    fCaloID(-1)
    {}

    //----------------------------------------------------------------------
    CaloRawDigit::CaloRawDigit(unsigned long long int cellID, double ADC, double time, double x, double y, double z, unsigned int layer, int id)
    : fCellID(cellID),
    fADC(ADC),
    fTime(time),
    fX(x),
    fY(y),
    fZ(z),
    fLayer(layer),
    fCaloID(id)
    {}

    } // namespace raw
  } // gar
    ////////////////////////////////////////////////////////////////////////
