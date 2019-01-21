//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"
#include "RawDataProducts/raw.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandBinomial.h"

#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TString.h"

namespace gar {
    namespace rosim{

        //----------------------------------------------------------------------------
        ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
        fhicl::ParameterSet    const& pset)
        : ECALReadoutSimAlg(engine, pset)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            this->reconfigure(pset);

            return;
        }

        //----------------------------------------------------------------------------
        ECALReadoutSimStandardAlg::~ECALReadoutSimStandardAlg()
        {
            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::reconfigure(fhicl::ParameterSet const& pset)
        {
            fAddNoise = pset.get<bool>("AddNoise", false);
            fNoiseMeV = pset.get<float>("NoiseMeV", 0.1);
            fSaturation = pset.get<bool>("Saturation", false);
            fTimeSmearing = pset.get<bool>("TimeSmearing", false);
            fTimeResolution = pset.get<double>("TimeResolution", 1.); // default time resolution in ns
            fADCSaturation = pset.get<int>("ADCSaturation", 4096); // default to 12-bit ADC's
            fMeVtoMIP = pset.get<double>("MeVtoMIP", 0.814); // default MeVtoMIP for 5 mm Scint
            fPosResolution = pset.get<double>("PositionResolution", 30); // default 3 cm

            fECALUtils = std::make_unique<util::ECALUtils>(fDetProp->EffectivePixel(), 0.95);

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::CreateCaloRawDigits(std::vector<sdp::CaloDeposit> SimCaloVec, std::vector<raw::CaloRawDigit> &digCol)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "CreateCaloRawDigits()";

            std::unordered_map<long long int, sdp::CaloDeposit> m_TileSimHits;
            std::unordered_map<long long int, sdp::CaloDeposit> m_StripSimHits;

            for(auto const &SimCaloHit : SimCaloVec)
            {
                long long int cellID = SimCaloHit.CellID();

                if(isTile(cellID))
                this->FillSimHitMap(cellID, SimCaloHit, m_TileSimHits);
                else{
                    // Naively add the simhits in the same strips
                    this->FillSimHitMap(cellID, SimCaloHit, m_StripSimHits);
                }
            }

            //Treating tiled hits
            for(auto const &SimCaloHit : m_TileSimHits)
            {
                float energy = SimCaloHit.second.Energy();
                float time = SimCaloHit.second.Time();
                float x = SimCaloHit.second.X();
                float y = SimCaloHit.second.Y();
                float z = SimCaloHit.second.Z();
                long long int cellID = SimCaloHit.second.CellID();

                if(fAddNoise)
                AddElectronicNoise(energy);

                DoPhotonStatistics(energy);

                if(fTimeSmearing)
                DoTimeSmearing(time);

                raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                digCol.emplace_back(digithit);
            }

            //Treating strips "naively smear energy and position"
            for(auto const &SimCaloHit : m_StripSimHits)
            {
                float energy = SimCaloHit.second.Energy();
                float time = SimCaloHit.second.Time();
                float x = SimCaloHit.second.X();
                float y = SimCaloHit.second.Y();
                float z = SimCaloHit.second.Z();
                long long int cellID = SimCaloHit.second.CellID();

                if(fAddNoise)
                AddElectronicNoise(energy);

                DoPhotonStatistics(energy);

                if(fTimeSmearing)
                DoTimeSmearing(time);

                DoPositionSmearing(x, y, z);

                raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                digCol.emplace_back(digithit);
            }
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoPhotonStatistics(float& energy)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
            CLHEP::RandBinomial BinomialRand(fEngine);

            //convertion from GeV to MIP
            energy /= fMeVtoMIP * CLHEP::MeV / CLHEP::GeV;
            //conversion to px
            energy *= fDetProp->LightYield();
            //Saturation
            if(fSaturation)
            energy = fECALUtils->Saturate(energy);

            //Binomial Smearing
            double prob = energy / fDetProp->EffectivePixel();
            energy += BinomialRand.shoot(fDetProp->EffectivePixel(), prob);

            //Convertion to ADC
            if(energy > 0)
            energy *= fDetProp->SiPMGain();
            else
            energy = 0.;

            if(energy > fADCSaturation)
            energy = fADCSaturation;

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoTimeSmearing(float& time)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
            CLHEP::RandGauss GausRand(fEngine);
            time = GausRand.fire(time, fTimeResolution);
            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::AddElectronicNoise(float& energy)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddElectronicNoise()";

            //Random noise shooting (Gaussian electronic noise)
            CLHEP::RandGauss GausRand(fEngine);
            energy += GausRand.shoot(0., fNoiseMeV * CLHEP::MeV / CLHEP::GeV);

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoPositionSmearing(float &x, float &y, float &z)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPositionSmearing()";

            //Random position smearing (Gaussian)
            CLHEP::RandGauss GausRand(fEngine);

            //Depends on the strip (X or Y orientation (in local frame))
            x = GausRand.fire(x, fPosResolution);
            y = GausRand.fire(y, fPosResolution);

            return;
        }

        //----------------------------------------------------------------------------
        bool ECALReadoutSimStandardAlg::isTile(long long int& cID)
        {
            int det_id = fGeo->getIDbyCellID(cID, "system");
            int layer = fGeo->getIDbyCellID(cID, "layer");
            bool isTile = fGeo->isTile(det_id, layer);

            return isTile;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::FillSimHitMap(long long int& cID, sdp::CaloDeposit const& SimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit>& m_SimHits)
        {
            if(m_SimHits.count(cID) == 0)
            {
                sdp::CaloDeposit hit = sdp::CaloDeposit(SimCaloHit.TrackID(), SimCaloHit.Time(), SimCaloHit.Energy(), const_cast<double*>(SimCaloHit.Pos()), cID);
                m_SimHits.emplace(cID, hit);
            }
            else{
                m_SimHits.at(cID) += SimCaloHit;
            }
        }

    } // rosim
} // gar
