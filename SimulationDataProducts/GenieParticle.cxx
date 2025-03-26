//
//  GenieParticle.cxx
//  garsoft-mrb
//
//  Created by E. Brianne
//

#include "SimulationDataProducts/GenieParticle.h"
#include <limits>

namespace gar {
  namespace sdp {

    //-------------------------------------------------
    GenieParticle::GenieParticle()
      : fInteractionIndex(0)
      , fIndex(0)
      , fPdg(0)
      , fStatus(0)
      , fName("")
      , fFirstMother(-1)
      , fLastMother(-1)
      , fFirstDaughter(-1)
      , fLastDaughter(-1)
      , fPx(0.)
      , fPy(0.)
      , fPz(0.)
      , fE(0.)
      , fMass(0.)
      , fRescatterCode(0)
    {}

  } //sdp
} // gar
