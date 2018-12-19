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

        fECALUtils = std::make_unique<util::ECALUtils>(fDetProp->EffectivePixel(), 0.95);

        fECALLayerThickness = fGeo->GetECALLayerThickness(); // in cm
        fEndcapStartXPosition = fGeo->GetECALEndcapStartPosition(); // in cm
        fECALRinner = fGeo->GetECALInnerBarrelRadius(); // in cm
        fECALRouter = fGeo->GetECALOuterBarrelRadius(); // in cm
        fECALPVThickness = fGeo->GetPVThickness(); // in cm

        fAbsorberThickness = fGeo->GetECALAbsorberThickness();
        fActiveMatThickness = fGeo->GetECALActiveMatThickness();

        fTPCOriginX = fGeo->TPCXCent();
        fTPCOriginY = fGeo->TPCYCent();
        fTPCOriginZ = fGeo->TPCZCent();

        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::CreateCaloRawDigits(std::vector<sdp::CaloDeposit> SimCaloVec, std::vector<raw::CaloRawDigit> &digCol)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "CreateCaloRawDigits()";

        for(auto const &SimCaloHit : SimCaloVec)
        {
          float energy = SimCaloHit.Energy();
          float time = SimCaloHit.Time();
          float x = SimCaloHit.X();
          float y = SimCaloHit.Y();
          float z = SimCaloHit.Z();
          long long int cellID = SimCaloHit.CellID();

          if(fAddNoise)
          this->AddElectronicNoise(energy);

          this->DoPhotonStatistics(energy);

          if(fTimeSmearing)
          this->DoTimeSmearing(time);

          raw::CaloRawDigit digithit = raw::CaloRawDigit(static_cast<unsigned int>(energy), time, x, y, z, cellID);
          digCol.emplace_back(digithit);
        }
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::DoPhotonStatistics(float &energy)
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
      void ECALReadoutSimStandardAlg::DoTimeSmearing(float &time)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "DoTimeSmearing()";
        CLHEP::RandGauss GausRand(fEngine);
        time = GausRand.fire(time, fTimeResolution);
        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::AddElectronicNoise(float &energy)
      {
        LOG_DEBUG("ECALReadoutSimStandardAlg") << "AddElectronicNoise()";

        //Random noise shooting (Gaussian electronic noise)
        CLHEP::RandGauss GaussRand(fEngine);
        energy += GaussRand.shoot(0., fNoiseMeV * CLHEP::MeV / CLHEP::GeV);

        return;
      }

    } // rosim
  } // gar
