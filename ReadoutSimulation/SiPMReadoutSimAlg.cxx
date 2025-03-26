//
//  ECALReadoutSimAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//
#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/SiPMReadoutSimAlg.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------------
    SiPMReadoutSimAlg::SiPMReadoutSimAlg(CLHEP::HepRandomEngine& engine,
                                         fhicl::ParameterSet const& /* pset */)
      : fEngine(engine), fDetProp(nullptr)
    {
      return;
    }

    //----------------------------------------------------------------------------
    SiPMReadoutSimAlg::~SiPMReadoutSimAlg()
    {
      return;
    }

  } // end rosim
} // gar
