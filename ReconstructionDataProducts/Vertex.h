//
//  Vertex.h
//
//  Created by Tom Junk on 7/11/16.
//  Not much track information is stored in the vertex data product -- use the associations with tracks to find them.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Vertex_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Vertex_h

#include "RtypesCore.h"
#include <stdint.h>
#include <iostream>
#include "IDNumberGen.h"



namespace gar {
  namespace rec {

    class Vertex {
    
    public:
      Vertex();
      
      // let the compiler provide the dtor
      
    private:
      static gar::rec::IDNumber const FirstNumber = 400000;
      gar::rec::IDNumber fIDnumero;
      
      float                     fPosition[3];    ///< position of vertex
      float                     fCovMat[3][3];   ///< uncertianties on the position
      double                    fTime;           ///< timestamp (ROOT's portable version of uint64_t)

#ifndef __GCCXML__
      
    public:
      
      Vertex(const float       *pos,
             const float       *covmat,
	     const double      time);

      bool operator==(const Vertex& rhs) const;
      bool operator!=(const Vertex& rhs) const;
      gar::rec::IDNumber getIDNumber() const;
      
      const float*        Position()  const;
      const float*        CovMat()    const;
      double              Time()      const;
      
#endif
      
    };
    
    inline const  float* Vertex::Position()  const { return &fPosition[0]; }
    inline const  float* Vertex::CovMat()    const { return &(fCovMat[0][0]); }
    inline double        Vertex::Time()      const { return fTime; }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Hit_h */
