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

    typedef long long int CellID_t;

    class CaloRawDigit {

    public:

        /// Default constructor: an empty raw digit with zeros put in for paraneters and an invalid channel
      CaloRawDigit();

#ifndef __GCCXML__
    public:

      CaloRawDigit(unsigned int ADC, float time, float x, float y, float z, CellID_t cellID);
      CaloRawDigit(unsigned int ADC, std::pair<float, float> time, float x, float y, float z, CellID_t cellID);
      CaloRawDigit(std::pair<unsigned int, unsigned int> ADC, std::pair<float, float> time, float x, float y, float z, CellID_t cellID);

      /// Reference to the compressed ADC count vector
      std::pair<unsigned int, unsigned int> ADC() const;
      /// Timestmap
      std::pair<float, float> Time() const;
      /// X position
      float X() const;
      /// Y position
      float Y() const;
      /// Z position
      float Z() const;
      /// cellID
      CellID_t CellID() const;


#endif // !__GCCXML__
    private:

      std::pair<unsigned int, unsigned int> fADC; ///< energy of the hit in ADC for both SiPM
      std::pair<float, float> fTime; ///< time of the hit for both SiPM
      float fX; ///< x of the hit
      float fY; ///< y of the hit
      float fZ; ///< z of the hit
      CellID_t fCellID;              ///< cellID1 of the hit based on 64 bits

    }; // class CaloRawDigit


  } // namespace raw
} // gar

//------------------------------------------------------------------------------
//--- inline implementation
//---
#ifndef __GCCXML__

inline std::pair<unsigned int, unsigned int>            gar::raw::CaloRawDigit::ADC()           const {return fADC;   }
inline std::pair<float, float>                          gar::raw::CaloRawDigit::Time()          const {return fTime;  }
inline float                                            gar::raw::CaloRawDigit::X()             const {return fX;	 }
inline float                                            gar::raw::CaloRawDigit::Y()             const {return fY;	 }
inline float                                            gar::raw::CaloRawDigit::Z()             const {return fZ;	 }
inline gar::raw::CellID_t                               gar::raw::CaloRawDigit::CellID()        const {return fCellID;}

#endif // !__GCCXML__

#endif // GAR_RAWDATA_CALORAWDIGIT_H

////////////////////////////////////////////////////////////////////////
