//
//  Shower.hpp
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//  Additions by Tom Junk, 2018
//

#ifndef Shower_hpp
#define Shower_hpp

#include <stdio.h>
#include "RtypesCore.h"
#include <stdint.h>

namespace gar {
  namespace rec {
    
    class Shower {
      
    public:
      Shower();
      
        // let the compiler provide the dtor
      
    private:
      
      float fEnergy;    ///< energy of the shower
      float fVertex[3]; ///< Shower vertex position
      float fVtxDir[3]; ///< Shower vertex direction
      double  fTime; ///< Timestamp
      
#ifndef __GCCXML__
      
    public:
      
      Shower(float  energy,
             float *vtx,
             float *vtxDir,
	     double time);
      
      const float* Vertex() const;
      const float* VtxDir() const;
      float const& Energy() const;
      double    Time()   const;
      
#endif
      
    };
    
    inline const float* Shower::Vertex() const { return fVertex; }
    inline const float* Shower::VtxDir() const { return fVtxDir; }
    inline float const& Shower::Energy() const { return fEnergy; }
    inline double Shower::Time() const { return fTime; }
    
  } // rec
} // gar


#endif /* Shower_hpp */

