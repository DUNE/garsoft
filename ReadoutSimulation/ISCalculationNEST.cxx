////////////////////////////////////////////////////////////////////////
/// \file  ISCalcuationNEST.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using nest
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "ReadoutSimulation/ISCalculationNEST.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------------
    ISCalculationNEST::ISCalculationNEST(CLHEP::HepRandomEngine& engine)
      : ISCalculation(), fNest(0), fEngine(engine)
    {
      return;
    }

    //----------------------------------------------------------------------------
    ISCalculationNEST::~ISCalculationNEST()
    {
      if (fNest) delete fNest;

      return;
    }

    //----------------------------------------------------------------------------
    void ISCalculationNEST::Initialize()
    {
      // \todo should ideally make the yield factor passed to the NestAlg ctor a parameter
      if (!fNest) fNest = new NestAlg(1., fEngine);

      // Set the step size to small value if NEST is chosen, per Matthew Szydagis,
      // "because without delta rays, the yields are wrong.  The ICARUS model that is
      // in GArSoft uses a fudge factor to compensate, but NEST is "purer" -- no
      // fudge factor. "
      fStepSize = 0.05 * CLHEP::micrometer;

      return;
    }

    //----------------------------------------------------------------------------
    void ISCalculationNEST::Reset()
    {
      fEnergyDeposit = 0.;
      fNumIonElectrons = 0.;
      fNumScintPhotons = 0.;

      return;
    }

    //----------------------------------------------------------------------------
    void ISCalculationNEST::CalculateIonizationAndScintillation(const gar::sdp::EnergyDeposit* dep)
    {
      // get a const representation of the track for this step
      //const G4Track track(*(step->GetTrack()));

      //fNest->CalculateIonizationAndScintillation(track, *step);

      // compare the energy deposition of this step to what is in the fNest object
      if (fNest->EnergyDeposition() != dep->Energy() / CLHEP::GeV)
        MF_LOG_WARNING("ISCalculationNest")
          << "NEST and G4 step depositions do not agree!\n"
          << fNest->EnergyDeposition() << " vs " << dep->Energy() / CLHEP::GeV;

      // Nest uses Geant units, GArSoft assumes energy is in units of GeV here
      fEnergyDeposit = fNest->EnergyDeposition() / CLHEP::GeV;
      fNumIonElectrons = fNest->NumberIonizationElectrons();
      fNumScintPhotons = fNest->NumberScintillationPhotons();

      return;
    }

  } // namespace
} // gar
