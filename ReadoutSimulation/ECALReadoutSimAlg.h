//
//  ECALReadoutSimAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIMULATION_ECALReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_ECALReadoutSimAlg_hpp

#include "SimulationDataProducts/CaloDeposit.h"
#include "RawDataProducts/CaloRawDigit.h"

#include "DetectorInfo/DetectorProperties.h"
#include "Geometry/GeometryCore.h"

#include <map>

namespace fhicl {
    class ParameterSet;
}

namespace gar {

    namespace rosim {

        class ECALReadoutSimAlg{

        public:

            ECALReadoutSimAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);

            virtual ~ECALReadoutSimAlg();

            virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

            virtual void PrepareAlgo(const std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector) = 0;

            virtual void DoDigitization() = 0;

            virtual std::vector<raw::CaloRawDigit*> GetDigitizedHits() = 0;

        protected:

            CLHEP::HepRandomEngine&                  fEngine;   ///< random number engine
            bool                                     fAddNoise; ///< flag to add noise or not
            float                                    fNoiseMeV; ///< Electronic noise to add in MeV
            bool                                     fSaturation; ///< flag for sipm saturation or not
            bool                                     fTimeSmearing; ///< flag for time smearing or not
            float                                    fPosResolution; ///< position resolution in mm

            const detinfo::DetectorProperties*       fDetProp;  ///< detector properties
            gar::geo::GeometryCore const*            fGeo;        ///< geometry information

        };

    } // end rosim

} // end gar


#endif /* GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp */
