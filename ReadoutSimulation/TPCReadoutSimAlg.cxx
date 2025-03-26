//
//  TPCReadoutSimAlg.cxx
//
//  Created by Brian Rebel on 2/14/17.
//
#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/TPCReadoutSimAlg.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------------
    TPCReadoutSimAlg::TPCReadoutSimAlg(CLHEP::HepRandomEngine& engine,
                                       fhicl::ParameterSet const& /* pset */)
      : fEngine(engine), fDetProp(nullptr)
    {
      return;
    }

    //----------------------------------------------------------------------------
    TPCReadoutSimAlg::~TPCReadoutSimAlg()
    {
      return;
    }

  } // end rosim
} // gar
