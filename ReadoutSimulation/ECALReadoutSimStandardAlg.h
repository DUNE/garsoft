//
//  ECALReadoutSimStandardAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIM_ECALReadoutSimStandardAlg_h
#define GAR_READOUTSIM_ECALReadoutSimStandardAlg_h

#include "ReadoutSimulation/ECALReadoutSimAlg.h"
#include "Utilities/ECALUtils.h"

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

        void CreateCaloRawDigits(std::vector<sdp::CaloDeposit> CaloVec, std::vector<raw::CaloRawDigit> &digCol);
        void DoPhotonStatistics(float &energy);
        void DoTimeSmearing(float &time);
        void AddElectronicNoise(float &energy);

        void reconfigure(fhicl::ParameterSet const& pset);

      private:

        std::unique_ptr<util::ECALUtils> fECALUtils;

        float fECALLayerThickness;
        float fEndcapStartXPosition;
        float fECALRinner;
        float fECALRouter;
        float fECALPVThickness;
        float fAbsorberThickness;
        float fActiveMatThickness;

        float fTPCOriginX;
        float fTPCOriginY;
        float fTPCOriginZ;
      };

    }
  }

  #endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
