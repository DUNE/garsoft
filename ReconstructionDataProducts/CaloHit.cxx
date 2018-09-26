//
//  Hit.cxx
//
//  Created by Brian Rebel on 10/6/16.
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/CaloHit.h"

namespace gar {
  namespace rec {

    //--------------------------------------------------------------------------
    CaloHit::CaloHit()
    {
      return;
    }

    //--------------------------------------------------------------------------
    CaloHit::CaloHit(double energy, double time, double *pos, int id)
    : fEnergy  (energy  )
    , fTime   (time   )
    , fId     (id  )
    {

      fPosition[0] = pos[0];
      fPosition[1] = pos[1];
      fPosition[2] = pos[2];

      return;
    }

    //--------------------------------------------------------------------------
    std::ostream& operator<< (std::ostream& o, gar::rec::CaloHit const& h)
    {

      o << "CaloHit "
      << "\n\tposition = ("
      << h.Position()[0]
      << ", "
      << h.Position()[1]
      << ", "
      << h.Position()[2]
      << ")"
      << "\n\tenergy = "
      << h.Energy()
      << "\n\t time: "
      << h.Time()
      << " id: "
      << h.ID();

      return o;
    }

  } // rec
} // gar
