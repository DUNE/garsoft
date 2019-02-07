//
//  CaloHit.h
//
//  Created by Eldwan Brianne on 08/29/18.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h

#include <iostream>

namespace gar {
  namespace rec {

    class CaloHit {

    public:
      CaloHit();

      // let the compiler provide the dtor

    private:

      float                            fEnergy;      ///< energy of the calo hit in GeV
      float                            fPosition[3]; ///< position of the calo hit in cm
      float                            fTime;        ///< time of the calo hit in ns
      long long int                    fCellID;   ///< cellID

#ifndef __GCCXML__

    public:

      CaloHit(float energy, float time, float *pos, long long int cellID);

      const float*                  Position()  const;
      float                         Energy()    const;
      float                         Time()      const;
      long long int                 CellID()        const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::CaloHit const& h);

#endif

    };

    inline float                         gar::rec::CaloHit::Energy()    const { return fEnergy;      }
    inline const float*                  gar::rec::CaloHit::Position()  const { return &fPosition[0]; }
    inline float                         gar::rec::CaloHit::Time()      const { return fTime;       }
    inline long long int                 gar::rec::CaloHit::CellID()        const { return fCellID;    }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h */
