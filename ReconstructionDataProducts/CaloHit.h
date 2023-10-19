//
//  CaloHit.h
//
//  Created by Eldwan Brianne on 08/29/18.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h

#include <iostream>

#include "Geometry/GeometryGAr.h"

#include "TVector3.h"

#include "IDNumberGen.h"

namespace gar {
    namespace rec {

        class CaloHit {

        public:
            CaloHit();

            // let the compiler provide the dtor
            unsigned int GetCellLengthScale() const;

            bool operator< (const CaloHit &rhs) const;

        private:
            static gar::rec::IDNumber const FirstNumber = 100200000;
            gar::rec::IDNumber fIDnumero;


            float                            fEnergy;      ///< energy of the calo hit in GeV
            float                            fPosition[3]; ///< position of the calo hit in cm
            std::pair<float, float>          fTime;        ///< time of the calo hit in ns
            raw::CellID_t                    fCellID;      ///< cellID
            unsigned int                     fLayer;       ///< Layer

            #ifndef __GCCXML__

        public:

            CaloHit(float energy, float time, float *pos, raw::CellID_t cellID, unsigned int layer);

            CaloHit(float energy, std::pair<float, float> time, float *pos, raw::CellID_t cellID, unsigned int layer);

            bool operator==(const CaloHit& rhs) const;
            bool operator!=(const CaloHit& rhs) const;
            gar::rec::IDNumber getIDNumber() const;

            const float*                  Position()  const;
            float                         Energy()    const;
            std::pair<float, float>       Time()      const;
            raw::CellID_t                 CellID()    const;
            unsigned int                  Layer()     const;

            friend std::ostream& operator << (std::ostream & o, gar::rec::CaloHit const& h);

            #endif

        };

        inline float                         gar::rec::CaloHit::Energy()                   const { return fEnergy;       }
        inline const float*                  gar::rec::CaloHit::Position()                 const { return &fPosition[0]; }
        inline std::pair<float, float>       gar::rec::CaloHit::Time()                     const { return fTime;         }
        inline raw::CellID_t                 gar::rec::CaloHit::CellID()                   const { return fCellID;       }
        inline unsigned int                  gar::rec::CaloHit::Layer()                    const { return fLayer;        }

    } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_CaloHit_h */
