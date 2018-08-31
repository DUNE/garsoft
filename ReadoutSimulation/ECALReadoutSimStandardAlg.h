//
//  ECALReadoutSimStandardAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIM_ECALReadoutSimStandardAlg_h
#define GAR_READOUTSIM_ECALReadoutSimStandardAlg_h

#include "ReadoutSimulation/ECALReadoutSimAlg.h"
#include "ReadoutSimulation/ECALUtils.h"

namespace fhicl{
  class ParameterSet;
}

namespace gar{
  namespace rosim{

    class ECALReadoutSimStandardAlg : public ECALReadoutSimAlg {

    public:

      ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
        fhicl::ParameterSet    const& pset);
        virtual ~ECALReadoutSimStandardAlg();

        raw::CaloRawDigit CreateCaloRawDigit(sdp::CaloDeposit const &SimCaloHit);
        void DoPhotonStatistics(double &energy);
        void DoTimeSmearing(double &time);
        void CreateNoiseDigits(std::vector<raw::CaloRawDigit> & digits);

        void reconfigure(fhicl::ParameterSet const& pset);

      private:

        std::unique_ptr<ECALUtils> fECALUtils;
      };

    }
  }

  #endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
