//
//  ECALReadoutSimAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIMULATION_ECALReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_ECALReadoutSimAlg_hpp

#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandGauss.h"

#include "SimulationDataProducts/CaloDeposit.h"
#include "RawDataProducts/CaloRawDigit.h"

#include "DetectorInfo/DetectorProperties.h"

namespace fhicl {
  class ParameterSet;
}

namespace gar {

  namespace rosim {

    class ECALReadoutSimAlg{

    public:

      ECALReadoutSimAlg(CLHEP::HepRandomEngine      & engine,
                       fhicl::ParameterSet    const& pset);

      virtual ~ECALReadoutSimAlg();

      // Method to take simulated calo hits and turn them into CaloRawDigit hits
      virtual raw::CaloRawDigit CreateCaloRawDigit(sdp::CaloDeposit const &SimCaloHit) = 0;
      //Photon statistics
      virtual void DoPhotonStatistics(double &energy) = 0;

      virtual void CreateNoiseDigits(std::vector<raw::CaloRawDigit> & digits)      = 0;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

    protected:

      CLHEP::HepRandomEngine &           fEngine;   ///< random number engine
      bool                               fAddNoise; ///< flag to add noise or not
      bool                               fSaturation; ///< flag for sipm saturation or not
      bool                               fTimeSmearing; ///< flag for time smearing or not
      double                             fTimeResolution; ///< time resolution in ns
      int                                fADCSaturation; ///< limit of the ADC
      double                             fMeVtoMIP; ///< Conversion from MeV to MIP
      const detinfo::DetectorProperties* fDetProp;  ///< detector properties

    };

  } // end rosim

} // end gar


#endif /* GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp */
