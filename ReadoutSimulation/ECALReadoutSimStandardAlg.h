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

        void CreateCaloRawDigits(std::vector<sdp::CaloDeposit> CaloVec, std::vector<raw::CaloRawDigit> &digCol);
        sdp::CaloDeposit Segmentation(sdp::CaloDeposit const &SimCaloHit);
        void AddSubHits(std::map<unsigned long long int, std::vector<sdp::CaloDeposit> > m_SimCaloHits, std::vector<sdp::CaloDeposit> &SimCaloHitVec);
        void DoPhotonStatistics(float &energy);
        void DoTimeSmearing(float &time);
        void CreateNoiseDigits(std::vector<raw::CaloRawDigit> & digits);

        void reconfigure(fhicl::ParameterSet const& pset);

      private:

        std::unique_ptr<ECALUtils> fECALUtils;

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

        //MultiSegmentation
        bool fMultiSeg;
        unsigned int fMultiSegLayer;
        float fMultiSegCellSize;
      };

    }
  }

  #endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
