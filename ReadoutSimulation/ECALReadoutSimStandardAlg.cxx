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
#include "Geometry/LocalTransformation.h"
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
            fPosResolution = pset.get<double>("PositionResolution", 3); // default 3 cm

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

                if(fGeo->isTile(cellID))
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
                this->AddElectronicNoise(energy);

                this->DoPhotonStatistics(energy, cellID);

                if(fTimeSmearing)
                this->DoTimeSmearing(time);

                raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                digCol.emplace_back(digithit);
            }

            //Treating strips "smear energy and naively smeared position"
            for(auto const &SimCaloHit : m_StripSimHits)
            {
                float energy = SimCaloHit.second.Energy();
                float time = SimCaloHit.second.Time();
                float x = SimCaloHit.second.X();
                float y = SimCaloHit.second.Y();
                float z = SimCaloHit.second.Z();
                long long int cellID = SimCaloHit.second.CellID();

                if(fAddNoise)
                this->AddElectronicNoise(energy);

                this->DoPhotonStatistics(energy, cellID);

                if(fTimeSmearing)
                this->DoTimeSmearing(time);

                this->DoPositionSmearing(x, y, z, cellID, this->isStripDirectionX(cellID));

                raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                digCol.emplace_back(digithit);
            }
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoPhotonStatistics(float& energy, const long long int& cID)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
            CLHEP::RandBinomial BinomialRand(fEngine);

            //convertion from GeV to MIP
            if(fGeo->isTile(cID))
            energy /= fMeVtoMIP * 2 * CLHEP::MeV / CLHEP::GeV; // Tiles are twice thicker than strips (TO DO: Lookup table for conversion factor)
            else
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

            if(energy > fDetProp->IntercalibrationFactor() * fADCSaturation)
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
        void ECALReadoutSimStandardAlg::DoPositionSmearing(float &x, float &y, float &z, const long long int& cID, bool isStripX)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPositionSmearing()";

            //Random position smearing (Gaussian)
            CLHEP::RandGauss GausRand(fEngine);

            //Find the volume path
            TVector3 point(x, y, z);
            std::string name = fGeo->VolumeName(point);
            auto const& path = fGeo->FindVolumePath(name);
            if (path.empty()) {
                throw cet::exception("ECALReadoutSimStandardAlg")
                << "DoPositionSmearing(): can't find volume '" << name << "'\n";
            }

            //Change to local frame
            gar::geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
            std::array<double, 3U> world{ {x, y, z} }, local;
            trans.WorldToLocal(world.data(), local.data());

            // std::cout << "Node " << name << std::endl;
            // trans.Matrix().Print();
            // std::cout << "world: x " << world[0] << " y " << world[1] << " z " << world[2] << std::endl;
            // std::cout << "local: x " << local[0] << " y " << local[1] << " z " << local[2] << std::endl;

            //Depends on the strip (X or Y orientation (in local frame))
            //The origin is at the center of the tile/strip
            double x_smeared, y_smeared, z_smeared = 0.;

            if(isStripX)
            {
                x_smeared = GausRand.fire(local[0], fPosResolution); //# in cm
                y_smeared = this->CalculateStripCenter(cID);
            } else {
                x_smeared = this->CalculateStripCenter(cID);
                y_smeared = GausRand.fire(local[1], fPosResolution); //# in cm
            }

            std::array<double, 3U> local_back{ {x_smeared, y_smeared, z_smeared} }, world_back;
            //Change back to World frame
            trans.LocalToWorld(local_back.data(), world_back.data());

            // std::cout << "local_back: x " << local_back[0] << " y " << local_back[1] << " z " << local_back[2] << std::endl;
            // std::cout << "world_back: x " << world_back[0] << " y " << world_back[1] << " z " << world_back[2] << std::endl;

            x = world_back[0];
            y = world_back[1];
            z = world_back[2];

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::FillSimHitMap(const long long int& cID, sdp::CaloDeposit const& SimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit>& m_SimHits)
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

        //----------------------------------------------------------------------------
        bool ECALReadoutSimStandardAlg::isStripDirectionX(const long long int& cID)
        {
            int cellX = fGeo->getIDbyCellID(cID, "cellX");
            // int cellY = fGeo->getIDbyCellID(cID, "cellY");

            if(cellX == 0)
            {
                return true;
            }
            else{
                return false;
            }
        }

        //----------------------------------------------------------------------------
        double ECALReadoutSimStandardAlg::CalculateStripCenter(const long long int& cID)
        {
            int cellX = fGeo->getIDbyCellID(cID, "cellX");
            int cellY = fGeo->getIDbyCellID(cID, "cellY");

            double stripSize = fGeo->getStripWidth() / CLHEP::cm; //strip width is returned in mm

            if(this->isStripDirectionX(cID))
            {
                return ( (cellY + 0.5) * stripSize );
            } else {
                return ( (cellX + 0.5) * stripSize );
            }
        }

    } // rosim
} // gar
