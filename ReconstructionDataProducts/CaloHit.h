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

      double        fEnergy;      ///< energy of the calo hit
      double        fPosition[3]; ///< position of the calo hit
      double        fTime;        ///< time of the calo hit
      int           fId;   ///< id of the calo hit

#ifndef __GCCXML__

    public:

      CaloHit(double energy, double time, double *pos, int id);

      const double*        Position()  const;
      double         Energy()    const;
      double         Time()      const;
      int            ID()        const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::CaloHit const& h);

#endif

    };

    inline double         gar::rec::CaloHit::Energy()    const { return fEnergy;      }
    inline const double*        gar::rec::CaloHit::Position()  const { return &fPosition[0]; }
    inline double         gar::rec::CaloHit::Time()      const { return fTime;       }
    inline int            gar::rec::CaloHit::ID()        const { return fId;    }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h */
