//
//  Shower.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 10/6/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//  Additions by Tom Junk, 2018

#include "ReconstructionDataProducts/Shower.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Shower::Shower()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    Shower::Shower(float  energy,
                   float *vtx,
                   float *vtxDir,
		   double time)
      : fEnergy(energy), fTime(time)
    {
      fVertex[0] = vtx[0];
      fVertex[1] = vtx[1];
      fVertex[2] = vtx[2];
      
      fVtxDir[0] = vtxDir[0];
      fVtxDir[1] = vtxDir[1];
      fVtxDir[2] = vtxDir[2];
      
      return;
    }
    
  } // rec
} // gar
