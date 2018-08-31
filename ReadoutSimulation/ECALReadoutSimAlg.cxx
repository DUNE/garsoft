//
//  ECALReadoutSimAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//
#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/ECALReadoutSimAlg.h"

namespace gar {
  namespace rosim{

    //----------------------------------------------------------------------------
    ECALReadoutSimAlg::ECALReadoutSimAlg(CLHEP::HepRandomEngine      & engine,
                                       fhicl::ParameterSet    const& pset)
    : fEngine (engine)
    , fDetProp(nullptr)
    {
      return;
    }

    //----------------------------------------------------------------------------
    ECALReadoutSimAlg::~ECALReadoutSimAlg()
    {
      return;
    }

  } // end rosim
} // gar
