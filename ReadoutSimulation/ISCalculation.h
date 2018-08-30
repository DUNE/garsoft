////////////////////////////////////////////////////////////////////////
/// \file  ISCalculation.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef ROSIMISCALCULATION_H
#define ROSIMISCALCULATION_H

#include "SimulationDataProducts/EnergyDeposit.h"

namespace gar {
  namespace rosim{

    class ISCalculation{

    public:

      ISCalculation();
      virtual ~ISCalculation();

      virtual void   Initialize()                        = 0;
      virtual void   Reset()                             = 0;
      virtual void   CalculateIonizationAndScintillation(const gar::sdp::EnergyDeposit* dep) = 0;
      virtual double EnergyDeposit()              const  = 0;
      virtual int    NumberIonizationElectrons()  const  = 0;
      virtual int    NumberScintillationPhotons() const  = 0;
      virtual double StepSizeLimit()              const  = 0;

    protected:

      double fEnergyDeposit;   ///< total energy deposited in the step
      int    fNumIonElectrons; ///< number of ionization electrons for this step
      int    fNumScintPhotons; ///< number of scintillation photons for this step

    };
  }
} // gar

#endif // ROSIMISCALCULATION_H
