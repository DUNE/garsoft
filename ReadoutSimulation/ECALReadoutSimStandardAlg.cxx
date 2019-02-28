//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/Geometry.h"
#include "Geometry/LocalTransformation.h"
#include "CoreUtils/ServiceUtil.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandBinomial.h"

namespace gar {
    namespace rosim{

        //----------------------------------------------------------------------------
        ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset)
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

        void ECALReadoutSimStandardAlg::ClearLists()
        {
            m_SimCaloHitList.clear();
            m_DigitHitVec.clear();
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::PrepareAlgo(const std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
        {
            //Clear the lists
            ClearLists();

            //Loop over all hits
            for (std::vector< art::Ptr<sdp::CaloDeposit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<sdp::CaloDeposit> hitPtr = *iter;
                const sdp::CaloDeposit *hit = hitPtr.get();
                m_SimCaloHitList.push_back(hit);
            }

            //Sort the list by time
            m_SimCaloHitList.sort();

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoDigitization()
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoDigitization()";

            std::unordered_map<long long int, sdp::CaloDeposit*> m_TileSimHits;
            std::unordered_map<long long int, sdp::CaloDeposit*> m_StripSimHits;

            for (const sdp::CaloDeposit *const pSimCaloHit : m_SimCaloHitList)
            {
                if(fGeo->isTile(pSimCaloHit->CellID()))
                {
                    // Naively add the simhits in the same tile
                    this->FillSimCaloHitMap(pSimCaloHit, m_TileSimHits);
                }
                else
                {
                    // Naively add the simhits in the same strips
                    this->FillSimCaloHitMap(pSimCaloHit, m_StripSimHits);
                }
            }

            //Treating tiled hits
            for(auto it : m_TileSimHits)
            {
                const sdp::CaloDeposit *SimCaloHit = it.second;

                float energy = SimCaloHit->Energy();
                float time = SimCaloHit->Time();
                float x = SimCaloHit->X();
                float y = SimCaloHit->Y();
                float z = SimCaloHit->Z();
                long long int cellID = SimCaloHit->CellID();

                if(fAddNoise)
                this->AddElectronicNoise(energy);

                this->DoPhotonStatistics(energy, cellID);

                if(fTimeSmearing)
                this->DoTimeSmearing(time);

                raw::CaloRawDigit *digithit = new raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                m_DigitHitVec.push_back(digithit);
            }

            //Treating strips "smear energy and naively smeared position"
            for(auto it : m_StripSimHits)
            {
                const sdp::CaloDeposit *SimCaloHit = it.second;

                float energy = SimCaloHit->Energy();
                float time = SimCaloHit->Time();
                float x = SimCaloHit->X();
                float y = SimCaloHit->Y();
                float z = SimCaloHit->Z();
                long long int cellID = SimCaloHit->CellID();

                if(fAddNoise)
                this->AddElectronicNoise(energy);

                this->DoPhotonStatistics(energy, cellID);

                if(fTimeSmearing)
                this->DoTimeSmearing(time);

                this->DoPositionSmearing(x, y, z, cellID, this->isStripDirectionX(cellID));

                raw::CaloRawDigit *digithit = new raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
                m_DigitHitVec.push_back(digithit);
            }

            m_TileSimHits.clear();
            m_StripSimHits.clear();
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::DoPhotonStatistics(float& energy, const long long int& cID)
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
            CLHEP::RandBinomial BinomialRand(fEngine);

            //convertion from GeV to MIP
            float energy_mip = 0.;
            if(fGeo->isTile(cID))
            energy_mip = energy / (fMeVtoMIP * 2 * CLHEP::MeV / CLHEP::GeV); // Tiles are twice thicker than strips (TO DO: Lookup table for conversion factor)
            else
            energy_mip = energy / (fMeVtoMIP * CLHEP::MeV / CLHEP::GeV);

            //conversion to px
            float pixel = energy_mip * fDetProp->LightYield();

            //Saturation
            float sat_pixel = 0.;
            if(fSaturation)
            sat_pixel = fECALUtils->Saturate(pixel);
            else
            sat_pixel = pixel;

            //Binomial Smearing
            double prob = sat_pixel / fDetProp->EffectivePixel();
            float smeared_px = BinomialRand.shoot(fDetProp->EffectivePixel(), prob);

            //Convertion to ADC
            if(energy > 0)
            energy = smeared_px * fDetProp->SiPMGain();
            else
            energy = 0.;

            if(energy > fDetProp->IntercalibrationFactor() * fADCSaturation)
            energy = fDetProp->IntercalibrationFactor() * fADCSaturation;

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
            if (path.empty())
            {
                throw cet::exception("ECALReadoutSimStandardAlg") << "DoPositionSmearing(): can't find volume '" << name << "'\n";
            }

            //Change to local frame
            gar::geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
            std::array<double, 3U> world{ {x, y, z} }, local;
            trans.WorldToLocal(world.data(), local.data());

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

            x = world_back[0];
            y = world_back[1];
            z = world_back[2];

            return;
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::FillSimCaloHitMap(const sdp::CaloDeposit *const pSimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit*>& m_SimCaloHits)
        {
            if(m_SimCaloHits.count(pSimCaloHit->CellID()) == 0)
            {
                sdp::CaloDeposit *hit = new sdp::CaloDeposit(pSimCaloHit->TrackID(), pSimCaloHit->Time(), pSimCaloHit->Energy(), const_cast<double*>(pSimCaloHit->Pos()), pSimCaloHit->CellID());

                //Add the new made hit to the map
                m_SimCaloHits.emplace(pSimCaloHit->CellID(), hit);

                std::vector<const sdp::CaloDeposit*> pSimVec;
                pSimVec.push_back(pSimCaloHit);
            }
            else
            {
                //The sim hit already exist... adding naively the sim hits (adding energy)
                *(m_SimCaloHits.at(pSimCaloHit->CellID())) += *pSimCaloHit;
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
