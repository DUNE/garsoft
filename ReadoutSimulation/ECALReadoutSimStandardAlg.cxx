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

#include "TGeoManager.h"

namespace gar {
    namespace rosim{

        //----------------------------------------------------------------------------
        ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset)
        : ECALReadoutSimAlg(engine, pset)
        {
            fGeo = gar::providerFrom<geo::Geometry>();
            fGeoManager = fGeo->ROOTGeoManager();
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
            fSaturation = pset.get<bool>("Saturation", false);
            fTimeSmearing = pset.get<bool>("TimeSmearing", false);

            fECALUtils = std::make_unique<util::ECALUtils>(fDetProp->EffectivePixel());

            return;
        }

        //----------------------------------------------------------------------------
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

                float new_energy = this->DoPhotonStatistics(x, y, z, energy);
                float new_time = time;

                if(fTimeSmearing)
                new_time = this->DoTimeSmearing(time);

                std::shared_ptr<raw::CaloRawDigit> digihit = std::make_shared<raw::CaloRawDigit>( raw::CaloRawDigit(static_cast<unsigned int>(new_energy), new_time, x, y, z, cellID) );

                LOG_DEBUG("ECALReadoutSimStandardAlg") << "digihit " << digihit.get()
                << " with cellID " << cellID
                << " has energy " << static_cast<unsigned int>(new_energy)
                << " time " << new_time << " ns"
                << " pos (" << x << ", " <<  y << ", " << z << ")";

                m_DigitHitVec.push_back(digihit);
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

                std::shared_ptr<raw::CaloRawDigit> digihit = this->DoStripDigitization(x, y, z, energy, time, cellID);
                m_DigitHitVec.push_back(digihit);
            }

            m_TileSimHits.clear();
            m_StripSimHits.clear();
        }

        //----------------------------------------------------------------------------
        std::shared_ptr<raw::CaloRawDigit> ECALReadoutSimStandardAlg::DoStripDigitization(float x, float y, float z, float energy, float time, long long int cID) const
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoStripDigitization()";

            //Calculate the light propagation along the strip
            std::pair<float, float> times = this->DoLightPropagation(x, y, z, time, cID);

            //Calculate the SiPM energy
            float new_energy = this->DoPhotonStatistics(x, y, z, energy);

            //Calculate the position of the strip
            std::array<double, 3U> pos = this->CalculateStripPosition(x, y, z, cID);

            //make the shared ptr
            std::shared_ptr<raw::CaloRawDigit> digihit = std::make_shared<raw::CaloRawDigit>( raw::CaloRawDigit(static_cast<unsigned int>(new_energy), times, pos[0], pos[1], pos[2], cID) );

            LOG_DEBUG("ECALReadoutSimStandardAlg") << "digihit " << digihit.get()
            << " with cellID " << cID
            << " has energy " << static_cast<unsigned int>(new_energy)
            << " time (" << times.first << ", " << times.second << ")"
            << " pos (" << pos[0] << ", " <<  pos[1] << ", " << pos[2] << ")";

            return digihit;
        }

        //----------------------------------------------------------------------------
        float ECALReadoutSimStandardAlg::DoPhotonStatistics(float x, float y, float z, float energy) const
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
            CLHEP::RandBinomial BinomialRand(fEngine);

            TVector3 point(x, y, z);
            float factor = fGeo->GetSensVolumeThickness(point) / 0.5;

            //Convert energy from GeV to MIP
            float energy_mip = energy / (fDetProp->MeVtoMIP() * factor * CLHEP::MeV / CLHEP::GeV);

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

            float smeared_px_noise = smeared_px;
            //Add noise
            if(fAddNoise)
            smeared_px_noise = this->AddElectronicNoise(smeared_px);

            //Convertion to ADC
            float ADC = 0.;

            if(smeared_px_noise > 0)
            ADC = smeared_px_noise * fDetProp->SiPMGain();

            if(ADC > fDetProp->IntercalibrationFactor() * fDetProp->ADCSaturation())
            ADC = fDetProp->IntercalibrationFactor() * fDetProp->ADCSaturation();

            return ADC;
        }

        //----------------------------------------------------------------------------
        float ECALReadoutSimStandardAlg::DoTimeSmearing(float time) const
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
            CLHEP::RandGauss GausRand(fEngine);
            float smeared_time = time + GausRand.fire(0., fDetProp->TimeResolution());
            return smeared_time;
        }

        //----------------------------------------------------------------------------
        float ECALReadoutSimStandardAlg::AddElectronicNoise(float energy) const
        {
            LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddElectronicNoise()";

            //Random noise shooting (Gaussian electronic noise)
            CLHEP::RandGauss GausRand(fEngine);
            float smeared_energy = energy + GausRand.shoot(0., fDetProp->NoisePx());
            return smeared_energy;
        }

        //----------------------------------------------------------------------------
        std::array<double, 3U> ECALReadoutSimStandardAlg::CalculateStripPosition(float x, float y, float z, long long int cID) const
        {
            //Find the volume path
            TVector3 point(x, y, z);
            std::string name = fGeo->VolumeName(point);
            auto const& path = fGeo->FindVolumePath(name);
            if (path.empty())
            {
                throw cet::exception("ECALReadoutSimStandardAlg") << "CalculateStripPosition(): can't find volume '" << name << "'\n";
            }

            //Change to local frame
            gar::geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
            std::array<double, 3U> world{ {x, y, z} }, local;
            trans.WorldToLocal(world.data(), local.data());

            //Use the segmentation algo to get the position
            TGeoNode *node = fGeoManager->FindNode(x, y, z);//Node in cm...
            G4ThreeVector pos = fGeo->position(node, cID);//returns in cm

            std::array<double, 3U> local_back{ {pos.x(), pos.y(), 0.} }, world_back;
            //Change back to World frame
            trans.LocalToWorld(local_back.data(), world_back.data());

            return world_back;
        }

        //----------------------------------------------------------------------------
        std::pair<float, float> ECALReadoutSimStandardAlg::DoLightPropagation(float x, float y, float z, float time, long long int cID) const
        {
            //Find the volume path
            TVector3 point(x, y, z);
            std::string name = fGeo->VolumeName(point);
            auto const& path = fGeo->FindVolumePath(name);
            if (path.empty())
            {
                throw cet::exception("ECALReadoutSimStandardAlg") << "DoLightPropagation(): can't find volume '" << name << "'\n";
            }

            //Change to local frame
            gar::geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
            std::array<double, 3U> world{ {x, y, z} }, local;
            trans.WorldToLocal(world.data(), local.data());

            //Calculate light propagation along the strip
            std::pair<float, float> times = fGeo->CalculateLightPropagation(point, local, cID);

            //t1 is left SiPM, t2 is right SiPM
            float time1 = time + times.first;
            float time2 = time + times.second;

            //Smear the times
            float smeared_time1 = this->DoTimeSmearing(time1);
            float smeared_time2 = this->DoTimeSmearing(time2);

            return std::make_pair(smeared_time1, smeared_time2);
        }

        //----------------------------------------------------------------------------
        void ECALReadoutSimStandardAlg::FillSimCaloHitMap(const sdp::CaloDeposit *const pSimCaloHit, std::unordered_map<long long int, sdp::CaloDeposit*>& m_SimCaloHits) const
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


    } // rosim
} // gar
