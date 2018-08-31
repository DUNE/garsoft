/** ****************************************************************************
 * @file CaloRawDigit.h
 * @brief Definition of basic calo raw digits
 * @author eldwan.brianne@desy.de
 * @see  CaloRawDigit.cxx raw.h
 * Created on the 08/27/18
 *
 * ****************************************************************************/

#ifndef GAR_RAWDATA_CALORAWDIGIT_H
#define GAR_RAWDATA_CALORAWDIGIT_H

// C/C++ standard libraries
#include <stdint.h> // uint32_t
#include <cstdlib> // size_t
#include <vector>
#include "RawDataProducts/RawTypes.h"
#include "RtypesCore.h"  // ULong64_t

/// Raw data description and utilities
namespace gar {
  namespace raw {

    class CaloRawDigit {

    public:

        /// Default constructor: an empty raw digit with zeros put in for paraneters and an invalid channel
      CaloRawDigit();

#ifndef __GCCXML__
    public:

      CaloRawDigit(double ADC, double time, double x, double y, double z, int id);

      /// Reference to the compressed ADC count vector
      double       ADC()     const;
      /// Timestmap
      double       Time()    const;
      /// X position
      double       X()    const;
      /// Y position
      double       Y()    const;
      /// Z position
      double       Z()    const;
      /// Calorimeter ID
      double       CaloId()    const;

#endif // !__GCCXML__
    private:
      double fADC; ///< energy of the hit in ADC
      double fTime; ///< time of the hit
      double fX; ///< x of the hit
      double fY; ///< y of the hit
      double fZ; ///< z of the hit
      int fCaloID; ///< id of the hit

    }; // class CaloRawDigit


  } // namespace raw
} // gar

//------------------------------------------------------------------------------
//--- inline implementation
//---
#ifndef __GCCXML__

inline double            gar::raw::CaloRawDigit::ADC()    const { return fADC;   }
inline double            gar::raw::CaloRawDigit::Time()        const { return fTime;       }
inline double            gar::raw::CaloRawDigit::X()        const { return fX;       }
inline double            gar::raw::CaloRawDigit::Y()        const { return fY;       }
inline double            gar::raw::CaloRawDigit::Z()        const { return fZ;       }
inline double            gar::raw::CaloRawDigit::CaloId()        const { return fCaloID;       }

#endif // !__GCCXML__

#endif // GAR_RAWDATA_CALORAWDIGIT_H

////////////////////////////////////////////////////////////////////////
