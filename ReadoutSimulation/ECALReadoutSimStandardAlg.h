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

            void DoPhotonStatistics(float &energy, const long long int& cID);

            void DoTimeSmearing(float &time);

            void DoPositionSmearing(float &x, float &y, float &z, const long long int& cID, bool isStripX);

            void AddElectronicNoise(float &energy);

            bool isStripDirectionX(const long long int& cID);

            void FillSimCaloHitMap(const sdp::CaloDeposit *const pSimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit*>& m_SimCaloHits);

            double CalculateStripCenter(const long long int& cID);

            std::vector<raw::CaloRawDigit*> GetDigitizedHits() { return m_DigitHitVec; }

        private:

            std::unique_ptr<util::ECALUtils> fECALUtils; ///<used for the SiPM saturation

            std::list<const sdp::CaloDeposit*> m_SimCaloHitList; ///<used to store the simulated hits

            std::vector<raw::CaloRawDigit*> m_DigitHitVec; ///<vector of digitized hits
        };

    }
}

#endif /* GAR_READOUTSIM_TPCReadoutSimStandardAlg_h */
