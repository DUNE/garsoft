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
        void DoPositionSmearing(float &x, float &y, float &z);
        void AddElectronicNoise(float &energy);

        void reconfigure(fhicl::ParameterSet const& pset);

        bool isTile(long long int &cID);
        void FillSimHitMap(long long int &cID, sdp::CaloDeposit const &SimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit> &m_SimHits);

      private:

        std::unique_ptr<util::ECALUtils> fECALUtils;
        // double fTileSize;   // clang complains about these two being unused so comment them out
        // double fStripWidth;  
      };

    }
  }

  #endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
