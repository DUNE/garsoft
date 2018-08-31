////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationNEST.h
/// \brief Interface to algorithm class for a specific calculation of
///        ionization electrons and scintillation photons using NEST
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef ROSIMISCALCULATIONNEST_H
#define ROSIMISCALCULATIONNEST_H

#include "ReadoutSimulation/ISCalculation.h"
#include "GArG4/NestAlg.h"

// forward declaration
namespace CLHEP { class HepRandomEngine; }

namespace gar {
  namespace rosim {

    class ISCalculationNEST : public ISCalculation {

    public:

      ISCalculationNEST(CLHEP::HepRandomEngine& engine);
      virtual ~ISCalculationNEST();

      void   Initialize();
      void   Reset();
      void   CalculateIonizationAndScintillation(const gar::sdp::EnergyDeposit* dep);
      double EnergyDeposit()              const { return fEnergyDeposit;   }
      int    NumberIonizationElectrons()  const { return fNumIonElectrons; }
      int    NumberScintillationPhotons() const { return fNumScintPhotons; }
      double StepSizeLimit()              const { return fStepSize;        }

    private:

      NestAlg* fNest;                  ///< the fast optical simulation process
      double   fStepSize;              ///< maximum step to take
      CLHEP::HepRandomEngine& fEngine; ///< random engine
    };
  }
} // gar

#endif // ROSIMISCALCULATIONNEST_H
