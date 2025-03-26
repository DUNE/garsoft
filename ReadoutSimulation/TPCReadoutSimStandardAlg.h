//
//  TPCReadoutSimStandardAlg.h
//
//  Created by Brian Rebel on 2/14/17.
//

#ifndef GAR_READOUTSIM_TPCReadoutSimStandardAlg_h
#define GAR_READOUTSIM_TPCReadoutSimStandardAlg_h

#include "ReadoutSimulation/TPCReadoutSimAlg.h"

namespace fhicl {
  class ParameterSet;
}

namespace gar {
  namespace rosim {

    class TPCReadoutSimStandardAlg : public TPCReadoutSimAlg {

    public:
      TPCReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);
      virtual ~TPCReadoutSimStandardAlg();

      raw::RawDigit CreateRawDigit(unsigned int channel,
                                   std::vector<float> const& electrons,
                                   bool& todrop);
      void CreateNoiseDigits(std::vector<raw::RawDigit>& digits);

      void reconfigure(fhicl::ParameterSet const& pset);

    private:
      // AddNoiseToSignalDigits is for channels where signal is recorded.
      // Assume that the noise level is dependent on the amount of signal recorded
      void AddNoiseToADCs(std::vector<short>& adcs);
      short ElectronsToADCs(float electrons);

      // pre-generated noise -- in ADC counts
      std::vector<short> fNoiseVec;
    };

  }
}

#endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
