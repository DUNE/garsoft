//
//  VecHit.cxx
//
//  Created by Tom Junk on Feb. 5, 2019
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/VecHit.h"
#include "TMath.h"

namespace gar {
  namespace rec {
   
    //--------------------------------------------------------------------------
    VecHit::VecHit()
    {
      return;
    }

    //--------------------------------------------------------------------------
    VecHit::VecHit(float       *pos,
                   float       *dir,
  	           float       len)
    : fLength      (len   )
    {
      
      fPosition[0] = pos[0];
      fPosition[1] = pos[1];
      fPosition[2] = pos[2];

      fDirection[0] = dir[0];
      fDirection[1] = dir[1];
      fDirection[2] = dir[2];
      
      return;
    }
    
    //--------------------------------------------------------------------------
    std::ostream& operator<< (std::ostream& o, gar::rec::VecHit const& h)
    {
      
      o << "Vector Hit "
      << "\n\tposition = ("
      << h.Position()[0]
      << ", "
      << h.Position()[1]
      << ", "
      << h.Position()[2]
      << ")"
      << "\tdirection = ("
      << h.Direction()[0]
      << ", "
      << h.Direction()[1]
      << ", "
      << h.Direction()[2]
      << ")"
      << "\tLength = " << h.Length();
      
      return o;
    }

  } // rec
} // gar
