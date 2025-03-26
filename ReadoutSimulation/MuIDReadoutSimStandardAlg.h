//
//  MuIDReadoutSimStandardAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIM_MuIDReadoutSimStandardAlg_h
#define GAR_READOUTSIM_MuIDReadoutSimStandardAlg_h

#include "ReadoutSimulation/SiPMReadoutSimAlg.h"

#include "RecoAlg/SiPMUtils.h"

#include <unordered_map>

namespace fhicl {
  class ParameterSet;
}

namespace gar {
  namespace rosim {

    class MuIDReadoutSimStandardAlg : public SiPMReadoutSimAlg {

    public:
      MuIDReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);

      virtual ~MuIDReadoutSimStandardAlg();

      void reconfigure(fhicl::ParameterSet const& pset);

      void ClearLists();

      void PrepareAlgo(const std::vector<art::Ptr<sdp::CaloDeposit>>& hitVector);

      void DoDigitization();

      std::vector<raw::CaloRawDigit*> GetDigitizedHits() const { return m_DigitHitVec; }

    protected:
      float DoPhotonStatistics(float x, float y, float z, float energy) const;

      float DoTimeSmearing(float time) const;

      float AddElectronicNoise(float energy) const;

      raw::CaloRawDigit* DoStripDigitization(float x,
                                             float y,
                                             float z,
                                             float energy,
                                             float time,
                                             raw::CellID_t cID) const;

      std::pair<std::array<double, 3>, bool> CalculatePosition(float x,
                                                               float y,
                                                               float z,
                                                               raw::CellID_t cID) const;

      std::pair<float, float> DoLightPropagation(float x,
                                                 float y,
                                                 float z,
                                                 float time,
                                                 raw::CellID_t cID) const;

    private:
      std::unique_ptr<util::SiPMUtils> fSiPMUtils; ///<used for the SiPM saturation

      std::vector<const sdp::CaloDeposit*> m_SimCaloHitVec; ///<used to store the simulated hits

      std::vector<raw::CaloRawDigit*> m_DigitHitVec; ///<vector of digitized hits

      TGeoManager* fGeoManager;
    };

  }
}

#endif /* GAR_READOUTSIM_MuIDReadoutSimStandardAlg_h */
