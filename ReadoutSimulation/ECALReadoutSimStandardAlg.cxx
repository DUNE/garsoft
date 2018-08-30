//
//  ECALReadoutSimStandardAlg.cxx
//
//  Created by Eldwan Brianne on 8/27/18
//

#include "fhiclcpp/ParameterSet.h"

#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "CoreUtils/ServiceUtil.h"
#include "RawDataProducts/raw.h"
#include "SimulationDataProducts/CaloDeposit.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandBinomial.h"

namespace gar {
  namespace rosim{

    //----------------------------------------------------------------------------
    ECALReadoutSimStandardAlg::ECALReadoutSimStandardAlg(CLHEP::HepRandomEngine      & engine,
      fhicl::ParameterSet    const& pset)
      : ECALReadoutSimAlg(engine, pset)
      {
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
        fTimeResolution = pset.get<double>("TimeResolution", 1); // default time resolution in ns
        fADCSaturation = pset.get<int>("ADCSaturation", 4095); // default to 12-bit ADC's
        fMeVtoMIP = pset.get<double>("MeVtoMIP", 0.814); // default MeVtoMIP for 5 mm Scint

        fECALUtils = std::make_unique<ECALUtils>(fDetProp->EffectivePixel(), 0.95);

        return;
      }

      //----------------------------------------------------------------------------
      raw::CaloRawDigit ECALReadoutSimStandardAlg::CreateCaloRawDigit(sdp::CaloDeposit const &SimCaloHit)
      {
        double energy = SimCaloHit.Energy();
        double x = SimCaloHit.X();
        double y = SimCaloHit.Y();
        double z = SimCaloHit.Z();
        double time = SimCaloHit.Time();
        int id = SimCaloHit.CaloID();

        this->DoPhotonStatistics(energy);

        if(fTimeSmearing)
        this->DoTimeSmearing(time);

        return raw::CaloRawDigit(energy, time, x, y, z, id);
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::DoPhotonStatistics(double &energy)
      {
        CLHEP::RandBinomial BinomialRand(fEngine);

        //convertion to MIP
        energy /= fMeVtoMIP;
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
      void ECALReadoutSimStandardAlg::DoTimeSmearing(double &time)
      {
        CLHEP::RandGauss GausRand(fEngine);
        time = GausRand.fire(time, fTimeResolution);
        return;
      }

      //----------------------------------------------------------------------------
      void ECALReadoutSimStandardAlg::CreateNoiseDigits(std::vector<raw::CaloRawDigit> & digits)
      {
        if(!fAddNoise) return;

        //TO DO

        return;
      }

    } // rosim
  } // gar
