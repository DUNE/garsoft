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

      CaloRawDigit(unsigned int ADC, float time, float x, float y, float z, unsigned int id, unsigned long long int cellID, unsigned int layer);

      /// Reference to the compressed ADC count vector
      unsigned int           ADC()     const;
      /// Timestmap
      float                  Time()    const;
      /// X position
      float                  X()    const;
      /// Y position
      float                  Y()    const;
      /// Z position
      float                  Z()    const;
      /// Calorimeter ID
      unsigned int           CaloID()    const;

      unsigned long long int CellID() const;

      unsigned int           Layer() const;

#endif // !__GCCXML__
    private:

      unsigned int fADC; ///< energy of the hit in ADC
      float fTime; ///< time of the hit
      float fX; ///< x of the hit
      float fY; ///< y of the hit
      float fZ; ///< z of the hit
      unsigned int fCaloID; ///< id of the hit
      unsigned long long int fCellID; ///< cellID of the hit based on 64 bits
      unsigned int fLayer; ///< layer of the hit

    }; // class CaloRawDigit


  } // namespace raw
} // gar

//------------------------------------------------------------------------------
//--- inline implementation
//---
#ifndef __GCCXML__

inline unsigned int            gar::raw::CaloRawDigit::ADC()           const { return fADC;   }
inline float                   gar::raw::CaloRawDigit::Time()          const { return fTime;       }
inline float                   gar::raw::CaloRawDigit::X()             const { return fX;       }
inline float                   gar::raw::CaloRawDigit::Y()             const { return fY;       }
inline float                   gar::raw::CaloRawDigit::Z()             const { return fZ;       }
inline unsigned int            gar::raw::CaloRawDigit::CaloID()        const { return fCaloID;       }
inline unsigned long long int  gar::raw::CaloRawDigit::CellID()        const { return fCellID;       }
inline unsigned int            gar::raw::CaloRawDigit::Layer()         const { return fLayer;       }

#endif // !__GCCXML__

#endif // GAR_RAWDATA_CALORAWDIGIT_H

////////////////////////////////////////////////////////////////////////
