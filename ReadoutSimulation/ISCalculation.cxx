////////////////////////////////////////////////////////////////////////
/// \file  ISCalculation.cxx
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "ReadoutSimulation/ISCalculation.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------
    ISCalculation::ISCalculation() : fEnergyDeposit(0), fNumIonElectrons(0), fNumScintPhotons(0) {}

    //----------------------------------------------------------------------
    ISCalculation::~ISCalculation() {}

  }
} // gar