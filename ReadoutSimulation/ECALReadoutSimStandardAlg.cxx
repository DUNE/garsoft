//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/GeometryGAr.h"
#include "Geometry/LocalTransformation.h"
#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "CLHEP/Random/RandBinomial.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "TGeoManager.h"

#include "canvas/Persistency/Common/Ptr.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------------
    ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine,
                                                         fhicl::ParameterSet const& pset)
      : SiPMReadoutSimAlg(engine, pset)
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();
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

      fSiPMUtils = std::make_unique<util::SiPMUtils>(fDetProp->EffectivePixel());

      return;
    }

    //----------------------------------------------------------------------------
    void ECALReadoutSimStandardAlg::ClearLists()
    {
      m_SimCaloHitVec.clear();
      m_DigitHitVec.clear();
    }

    //----------------------------------------------------------------------------
    void ECALReadoutSimStandardAlg::PrepareAlgo(
      const std::vector<art::Ptr<sdp::CaloDeposit>>& hitVector)
    {
      //Clear the lists
      ClearLists();

      //Loop over all hits
      for (std::vector<art::Ptr<sdp::CaloDeposit>>::const_iterator iter = hitVector.begin(),
                                                                   iterEnd = hitVector.end();
           iter != iterEnd;
           ++iter) {
        art::Ptr<sdp::CaloDeposit> hitPtr = *iter;
        const sdp::CaloDeposit* hit = hitPtr.get();
        m_SimCaloHitVec.emplace_back(hit);
      }

      //Sort the list by time
      std::sort(m_SimCaloHitVec.begin(),
                m_SimCaloHitVec.end(),
                [](const sdp::CaloDeposit* rha, const sdp::CaloDeposit* rhb) {
                  return rha->Time() < rhb->Time();
                });

      return;
    }

    //----------------------------------------------------------------------------
    void ECALReadoutSimStandardAlg::DoDigitization()
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoDigitization()";

      //Treating tiled hits
      for (auto const& it : m_SimCaloHitVec) {
        float energy = it->Energy();
        float time = it->Time();
        float x = it->X();
        float y = it->Y();
        float z = it->Z();
        std::array<double, 3> point = {x, y, z};
        raw::CellID_t cellID = it->CellID();

        raw::CaloRawDigit* digihit = nullptr;

        if (fGeo->isTile(point, cellID)) {
          digihit = this->DoTileDigitization(x, y, z, energy, time, cellID);
        }
        else {
          digihit = this->DoStripDigitization(x, y, z, energy, time, cellID);
        }

        if (nullptr != digihit) { m_DigitHitVec.emplace_back(digihit); }
        else {
          MF_LOG_DEBUG("ECALReadoutSimStandardAlg")
            << "Could not digitize the simulated hit " << it;
        }
      }
    }

    //----------------------------------------------------------------------------
    raw::CaloRawDigit* ECALReadoutSimStandardAlg::DoTileDigitization(float x,
                                                                     float y,
                                                                     float z,
                                                                     float energy,
                                                                     float time,
                                                                     raw::CellID_t cID) const
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTileDigitization()";

      float new_energy = this->DoPhotonStatistics(x, y, z, energy);
      float new_time = time;

      if (fTimeSmearing) new_time = this->DoTimeSmearing(time);

      //Calculate the position of the tile
      std::pair<std::array<double, 3>, bool> calc_pos = this->CalculatePosition(x, y, z, cID);
      //Check if need to drop the hit
      if (calc_pos.second) { return nullptr; }

      std::array<double, 3> pos = calc_pos.first;

      raw::CaloRawDigit* digihit = new raw::CaloRawDigit(
        static_cast<unsigned int>(new_energy), new_time, pos[0], pos[1], pos[2], cID);

      MF_LOG_DEBUG("ECALReadoutSimStandardAlg")
        << "Tile digihit " << digihit << " with cellID " << cID << " has energy "
        << static_cast<unsigned int>(new_energy) << " time " << new_time << " ns"
        << " pos (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")";

      return digihit;
    }

    //----------------------------------------------------------------------------
    raw::CaloRawDigit* ECALReadoutSimStandardAlg::DoStripDigitization(float x,
                                                                      float y,
                                                                      float z,
                                                                      float energy,
                                                                      float time,
                                                                      raw::CellID_t cID) const
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoStripDigitization()";

      //Calculate the light propagation along the strip
      std::pair<float, float> times = this->DoLightPropagation(x, y, z, time, cID);

      //Calculate the SiPM energy
      float new_energy = this->DoPhotonStatistics(x, y, z, energy);

      //Calculate the position of the strip
      std::pair<std::array<double, 3>, bool> calc_pos = this->CalculatePosition(x, y, z, cID);
      //Check if need to drop the hit
      if (calc_pos.second) { return nullptr; }

      std::array<double, 3> pos = calc_pos.first;

      //make the shared ptr
      raw::CaloRawDigit* digihit = new raw::CaloRawDigit(
        static_cast<unsigned int>(new_energy), times, pos[0], pos[1], pos[2], cID);

      MF_LOG_DEBUG("ECALReadoutSimStandardAlg")
        << "Strip digihit " << digihit << " with cellID " << cID << " has energy "
        << static_cast<unsigned int>(new_energy) << " time (" << times.first << ", " << times.second
        << ")"
        << " pos (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")";

      return digihit;
    }

    //----------------------------------------------------------------------------
    float ECALReadoutSimStandardAlg::DoPhotonStatistics(float x,
                                                        float y,
                                                        float z,
                                                        float energy) const
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoPhotonStatistics()";
      CLHEP::RandBinomial BinomialRand(fEngine);

      TVector3 point(x, y, z);
      float factor = fGeo->GetSensVolumeThickness(point) / 0.5;

      //Convert energy from GeV to MIP
      float energy_mip = energy / (fDetProp->MeVtoMIP() * factor * CLHEP::MeV / CLHEP::GeV);

      //conversion to px
      float pixel = energy_mip * fDetProp->LightYield();

      //Saturation
      float sat_pixel = 0.;
      if (fSaturation)
        sat_pixel = fSiPMUtils->Saturate(pixel);
      else
        sat_pixel = pixel;

      //Binomial Smearing
      double prob = sat_pixel / fDetProp->EffectivePixel();
      float smeared_px = BinomialRand.shoot(fDetProp->EffectivePixel(), prob);

      float smeared_px_noise = smeared_px;
      //Add noise
      if (fAddNoise) smeared_px_noise = this->AddElectronicNoise(smeared_px);

      //Convertion to ADC
      float ADC = 0.;

      if (smeared_px_noise > 0) ADC = smeared_px_noise * fDetProp->SiPMGain();

      if (ADC > fDetProp->IntercalibrationFactor() * fDetProp->ADCSaturation())
        ADC = fDetProp->IntercalibrationFactor() * fDetProp->ADCSaturation();

      return ADC;
    }

    //----------------------------------------------------------------------------
    float ECALReadoutSimStandardAlg::DoTimeSmearing(float time) const
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
      CLHEP::RandGauss GausRand(fEngine);
      float smeared_time = time + GausRand.fire(0., fDetProp->TimeResolution());
      return smeared_time;
    }

    //----------------------------------------------------------------------------
    float ECALReadoutSimStandardAlg::AddElectronicNoise(float energy) const
    {
      MF_LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddElectronicNoise()";

      //Random noise shooting (Gaussian electronic noise)
      CLHEP::RandGauss GausRand(fEngine);
      float smeared_energy = energy + GausRand.shoot(0., fDetProp->NoisePx());
      return smeared_energy;
    }

    //----------------------------------------------------------------------------
    std::pair<std::array<double, 3>, bool>
    ECALReadoutSimStandardAlg::CalculatePosition(float x, float y, float z, raw::CellID_t cID) const
    {
      bool drop = false;

      //Use the segmentation algo to get the position
      std::array<double, 3> point = {x, y, z};
      TGeoNode* node = fGeo->FindNode(point); //Node in cm...
      std::string nodename = node->GetName();

      std::array<double, 3> pointLocal;
      gar::geo::LocalTransformation<TGeoHMatrix> trans;
      fGeo->WorldToLocal(point, pointLocal, trans);

      std::array<double, 3> pointLocal_back = fGeo->GetPosition(node, cID); //returns in cm
      std::array<double, 3> point_back;
      fGeo->LocalToWorld(pointLocal_back, point_back, trans);
      TGeoNode* new_node = fGeo->FindNode(point_back); //Node in cm...
      std::string newnodename = new_node->GetName();

      //get the base of the node names (slice can be added in the new node)
      std::string base_node_name = nodename.substr(0, nodename.find("_layer_") + 9);
      std::string base_new_node_name = newnodename.substr(0, newnodename.find("_layer_") + 9);

      if (base_new_node_name != base_node_name) {
        MF_LOG_DEBUG("ECALReadoutSimStandardAlg")
          << "CalculatePosition()"
          << " isTile " << fGeo->isTile(point, cID) << "\n"
          << " Strip length " << fGeo->getStripLength(point, cID) << "\n"
          << " Local Point before new position ( " << pointLocal[0] << ", " << pointLocal[1] << ", "
          << pointLocal[2] << " ) in node " << nodename << "\n"
          << " Local Point after new position ( " << pointLocal_back[0] << ", "
          << pointLocal_back[1] << ", " << pointLocal_back[2] << " ) in node " << newnodename
          << "\n"
          << " Dropping the hit ";

        //Drop the hit
        drop = true;
      }

      return std::make_pair(point_back, drop);
    }

    //----------------------------------------------------------------------------
    std::pair<float, float> ECALReadoutSimStandardAlg::DoLightPropagation(float x,
                                                                          float y,
                                                                          float z,
                                                                          float time,
                                                                          raw::CellID_t cID) const
    {
      //Find the volume path
      std::array<double, 3> point = {x, y, z};
      std::array<double, 3> pointLocal;
      gar::geo::LocalTransformation<TGeoHMatrix> trans;
      fGeo->WorldToLocal(point, pointLocal, trans);

      //Calculate light propagation along the strip
      std::pair<float, float> times = fGeo->CalculateLightPropagation(point, pointLocal, cID);

      //t1 is left SiPM, t2 is right SiPM
      float time1 = time + times.first;
      float time2 = time + times.second;

      //Smear the times
      float smeared_time1 = this->DoTimeSmearing(time1);
      float smeared_time2 = this->DoTimeSmearing(time2);

      return std::make_pair(smeared_time1, smeared_time2);
    }

  } // rosim
} // gar
