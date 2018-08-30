////////////////////////////////////////////////////////////////////////
/// \file  ISCalculationSeparate.h
/// \brief Interface to algorithm class for a specific calculation of
///        ionization electrons and scintillation photons assuming there
///        is no correlation between the two
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef ROSIMISCALCULATIONSEPARATE_H
#define ROSIMISCALCULATIONSEPARATE_H

#include <map>

#include "ReadoutSimulation/ISCalculation.h"

// forward declaration
namespace CLHEP { class HepRandomEngine; }

namespace gar {
  namespace rosim {

    class ISCalculationSeparate : public ISCalculation {

    public:

      ISCalculationSeparate(CLHEP::HepRandomEngine&);
      virtual ~ISCalculationSeparate();

      void   Initialize();
      void   Reset();
      void   CalculateIonizationAndScintillation(const gar::sdp::EnergyDeposit* dep);
      double EnergyDeposit()              const { return fEnergyDeposit;       }
      int    NumberIonizationElectrons()  const { return fNumIonElectrons;     }
      int    NumberScintillationPhotons() const { return fNumScintPhotons;     }
      double StepSizeLimit()              const { return fStepSize;            }

    private:

      double          fStepSize;            ///< maximum step to take
      double 	   	    fEfield;              ///< value of electric field from GArProperties service
      double 	   	    fGeVToElectrons;      ///< conversion factor from GArProperties service
      double 	   	    fRecombA;             ///< from GArG4Parameters service
      double 	   	    fRecombk;             ///< from GArG4Parameters service
      // clang says these are unused -- may need them later
      //bool   	   	    fScintByParticleType; ///< from GArProperties service
      //double 	   	    fScintYieldFactor;    ///< scintillation yield factor
    //G4EmSaturation* fEMSaturation;        ///< pointer to EM saturation

    };
  }
} // gar

#endif // ROSIMISCALCULATIONSEPARATE_H
