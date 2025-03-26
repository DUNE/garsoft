//
//  ECALReadoutSimAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIMULATION_SiPMReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_SiPMReadoutSimAlg_hpp

#include "RawDataProducts/CaloRawDigit.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "DetectorInfo/DetectorProperties.h"
#include "Geometry/GeometryCore.h"

#include <map>

#include "CLHEP/Random/RandGauss.h"

namespace fhicl {
  class ParameterSet;
}

namespace gar {

  namespace rosim {

    class SiPMReadoutSimAlg {

    public:
      SiPMReadoutSimAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);

      virtual ~SiPMReadoutSimAlg();

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

      virtual void PrepareAlgo(const std::vector<art::Ptr<sdp::CaloDeposit>>& hitVector) = 0;

      virtual void DoDigitization() = 0;

      virtual std::vector<raw::CaloRawDigit*> GetDigitizedHits() const = 0;

    protected:
      CLHEP::HepRandomEngine& fEngine; ///< random number engine
      bool fAddNoise;                  ///< flag to add noise or not
      bool fSaturation;                ///< flag for sipm saturation or not
      bool fTimeSmearing;              ///< flag for time smearing or not

      const detinfo::DetectorProperties* fDetProp; ///< detector properties
      gar::geo::GeometryCore const* fGeo;          ///< geometry information
    };

  } // end rosim

} // end gar

#endif /* GAR_READOUTSIMULATION_SiPMReadoutSimAlg_hpp */
