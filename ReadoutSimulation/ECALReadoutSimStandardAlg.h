//
//  ECALReadoutSimStandardAlg.h
//
//  Created by Eldwan Brianne on 8/27/18.
//

#ifndef GAR_READOUTSIM_ECALReadoutSimStandardAlg_h
#define GAR_READOUTSIM_ECALReadoutSimStandardAlg_h

#include "ReadoutSimulation/ECALReadoutSimAlg.h"

#include "Utilities/ECALUtils.h"

#include <unordered_map>

namespace fhicl{
    class ParameterSet;
}

namespace gar{
    namespace rosim{

        class ECALReadoutSimStandardAlg : public ECALReadoutSimAlg {

        public:

            ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);

            virtual ~ECALReadoutSimStandardAlg();

            void reconfigure(fhicl::ParameterSet const& pset);

            void ClearLists();

            void PrepareAlgo(const std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

            void DoDigitization();

            std::vector< std::shared_ptr<raw::CaloRawDigit> > GetDigitizedHits() { return m_DigitHitVec; }

        protected:

            float DoPhotonStatistics(float x, float y, float z, float energy) const;

            float DoTimeSmearing(float time) const;

            float AddElectronicNoise(float energy) const;

            std::shared_ptr<raw::CaloRawDigit> DoStripDigitization(float x, float y, float z, float energy, float time, long long int cID) const;

            std::array<double, 3U> CalculateStripPosition(float x, float y, float z, long long int cID) const;

            std::pair<float, float> DoLightPropagation(float x, float y, float z, float time, long long int cID) const;

            void FillSimCaloHitMap(const sdp::CaloDeposit *const pSimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit*>& m_SimCaloHits) const;

        private:

            std::unique_ptr<util::ECALUtils> fECALUtils; ///<used for the SiPM saturation

            std::list<const sdp::CaloDeposit*> m_SimCaloHitList; ///<used to store the simulated hits

            std::vector< std::shared_ptr<raw::CaloRawDigit> > m_DigitHitVec; ///<vector of digitized hits

            TGeoManager* fGeoManager;
        };

    }
}

#endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
