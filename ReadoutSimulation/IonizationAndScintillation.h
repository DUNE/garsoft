////////////////////////////////////////////////////////////////////////
/// \file  IonizationAndScintillation.h
/// \brief Singleton to access a unified treatment of ionization and
///        scintillation in LAr
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef ROSIMIONIZATIONANDSCINTILLATION_H
#define ROSIMIONIZATIONANDSCINTILLATION_H

#include <cstring>
#include <memory>

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/ISCalculation.h"
#include "SimulationDataProducts/EnergyDeposit.h"

#include "TH1.h"
#include "TH2.h"

namespace CLHEP { class HepRandomEngine; }

namespace gar {
  namespace rosim {

    // The Ionization and Scintillation singleton
    class IonizationAndScintillation
    {
    public:

      static IonizationAndScintillation* CreateInstance(CLHEP::HepRandomEngine       & engine,
        fhicl::ParameterSet     const& pset);
        static IonizationAndScintillation* Instance();

        // Method to reset the internal variables held in the ISCalculation
        // This method should be called at the start of any G4Step
        void Reset(const gar::sdp::EnergyDeposit* dep);
        void reconfigure(fhicl::ParameterSet const& pset);

        double EnergyDeposit()              const { return fISCalc->EnergyDeposit();              }
        int    NumberIonizationElectrons()  const { return fISCalc->NumberIonizationElectrons();  }
        int    NumberScintillationPhotons() const { return fISCalc->NumberScintillationPhotons(); }
        double StepSizeLimit()              const { return fISCalc->StepSizeLimit();              }

      private:

        IonizationAndScintillation(CLHEP::HepRandomEngine& engine);
        ~IonizationAndScintillation();

        rosim::ISCalculation*          fISCalc;             ///< object to calculate ionization and scintillation
        ///< produced by an energy deposition
        std::string                    fISCalculator;       ///< name of calculator to use, NEST or Separate
        const gar::sdp::EnergyDeposit* fEnergyDeposit;      ///< reference to the current energy deposit
        int                            fStepNumber;         ///< last StepNumber checked
        int                            fTrkID;              ///< last TrkID checked

        TH1F*                          fStepSize;           ///< histogram of the step sizes
        TH1F*                          fElectronsPerStep;   ///< histogram of electrons per step
        TH1F*                          fPhotonsPerStep;     ///< histogram of the photons per step
        TH1F*                          fElectronsPerLength; ///< histogram of electrons per dX
        TH1F*                          fPhotonsPerLength;   ///< histogram of the photons per dX
        TH1F*                          fEnergyPerStep;      ///< histogram of the energy deposited per step
        TH1F*                          fElectronsPerEDep;   ///< histogram of electrons per MeV deposited
        TH1F*                          fPhotonsPerEDep;     ///< histogram of photons per MeV deposited
        TH2F*                          fElectronsVsPhotons; ///< histogram of electrons vs photons per step
        CLHEP::HepRandomEngine&        fEngine;             ///< random engine
      };

    } // namespace rosim

  } // gar

  #endif // ROSIMIONIZATIONANDSCINTILLATION
