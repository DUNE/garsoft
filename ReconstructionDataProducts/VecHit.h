//
//  VecHit.h
//
//  Created by Tom Junk on Feb. 5, 2019.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_VecHit_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_VecHit_h

#include <iostream>
#include <TVector3.h>

namespace gar {
  namespace rec {

    class VecHit {
    
    public:
      VecHit();
      
      // let the compiler provide the dtor
      
    private:
      
      float        fPosition[3];   ///< position of the vector hit
      float        fDirection[3];  ///< direction cosines
      float        fLength;        ///< length of VH (centered on fPosition, half of the length is in either direction)


#ifndef __GCCXML__
      
    public:
      
      VecHit(float       *pos,
             float       *dir,
	     float       len);
      
      const float*        Position()   const;
      const float*        Direction()  const;
      float               Length()     const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::VecHit const& h);
      
#endif
      
    };
    
    inline const  float* VecHit::Position()   const { return &fPosition[0]; }
    inline const  float* VecHit::Direction()  const { return &fDirection[0]; }
    inline float         VecHit::Length()     const { return fLength;        }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_VecHit_h */
