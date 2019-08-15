//
//  CaloHit.h
//
//  Created by Eldwan Brianne on 08/29/18.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h

#include <iostream>

#include "Geometry/Geometry.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "IDNumberGen.h"



namespace gar {
  namespace rec {

    class CaloHit {

    public:
      CaloHit();

      // let the compiler provide the dtor

      const unsigned int GetLayer() const;
      const unsigned int GetCellLengthScale() const;

      bool operator< (const CaloHit &rhs) const;

    private:
      static IDNumberGen::IDNumber const FirstNumber = 100200000;
      IDNumberGen::IDNumber fIDnumero;


      float                            fEnergy;      ///< energy of the calo hit in GeV
      float                            fPosition[3]; ///< position of the calo hit in cm
      float                            fTime;        ///< time of the calo hit in ns
      long long int                    fCellID;      ///< cellID
      CLHEP::Hep3Vector                fPositionVector;

#ifndef __GCCXML__

    public:

      CaloHit(float energy, float time, float *pos, long long int cellID);

      bool operator==(const CaloHit& rhs) const;
      bool operator!=(const CaloHit& rhs) const;
      IDNumberGen::IDNumber getIDNumber() const;

      const float*                  Position()  const;
      float                         Energy()    const;
      float                         Time()      const;
      long long int                 CellID()    const;
      const CLHEP::Hep3Vector&      GetPositionVector() const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::CaloHit const& h);

#endif

    };

    inline float                         gar::rec::CaloHit::Energy()                   const { return fEnergy;      }
    inline const float*                  gar::rec::CaloHit::Position()                 const { return &fPosition[0]; }
    inline float                         gar::rec::CaloHit::Time()                     const { return fTime;       }
    inline long long int                 gar::rec::CaloHit::CellID()                   const { return fCellID;    }
    inline const CLHEP::Hep3Vector&      gar::rec::CaloHit::GetPositionVector()        const { return fPositionVector;    }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h */
