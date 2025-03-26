//
//  TPCReadoutSimStandardAlg.cxx
//
//  Created by Brian Rebel on 2/14/17.
//

#include "fhiclcpp/ParameterSet.h"

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "RawDataProducts/raw.h"
#include "ReadoutSimulation/TPCReadoutSimStandardAlg.h"

namespace gar {
  namespace rosim {

    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::TPCReadoutSimStandardAlg(CLHEP::HepRandomEngine& engine,
                                                       fhicl::ParameterSet const& pset)
      : TPCReadoutSimAlg(engine, pset)
    {
      this->reconfigure(pset);

      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

      fNoiseVec.clear();
      std::vector<double> noisevecdbl(fNoiseVecSize, 0);
      if (fNoiseSpectrum == 0) {
        CLHEP::RandGauss GaussRand(fEngine);
        //std::cout << "noise amplitude: " << fNoiseAmplitude << std::endl;
        GaussRand.fireArray(fNoiseVecSize, &noisevecdbl[0], 0.0, fNoiseAmplitude);
        for (int i = 0; i < fNoiseVecSize; ++i) {
          long inoise = lrint(noisevecdbl[i]);
          if (inoise > 32767) inoise = 32767;   // to be stored in a signed short.
          if (inoise < -32767) inoise = -32767; // to be stored in a signed short.
          fNoiseVec.push_back((short)inoise);
          //std::cout << "inoise: " << i << " " << inoise << std::endl;
        }
      }

      return;
    }

    //----------------------------------------------------------------------------
    TPCReadoutSimStandardAlg::~TPCReadoutSimStandardAlg()
    {
      return;
    }

    //----------------------------------------------------------------------------
    raw::RawDigit TPCReadoutSimStandardAlg::CreateRawDigit(unsigned int channel,
                                                           std::vector<float> const& electrons,
                                                           bool& todrop)
    {
      // make a vector to store the adc information
      auto numTicks = fDetProp->NumberTimeSamples();
      std::vector<short> adcs(numTicks, 0);

      // check that the size of the electrons vector is the same as the
      // adc vector, they better be
      if (adcs.size() != electrons.size()) {
        throw cet::exception("TPCReadoutSimStandardAlg")
          << "The vectors for adc and electrons have different sizes, no idea "
          << "how to handle that, bail";
      }

      for (size_t e = 0; e < electrons.size(); ++e) {
        adcs[e] = this->ElectronsToADCs(electrons[e]);
        //        if(electrons[e] > 0)
        //          LOG_VERBATIM("TPCReadoutSim")
        //          << "ADC: "
        //          << adcs[e]
        //          << " electrons "
        //          << electrons[e]
        //          << " tdc "
        //          << e;
      }

      // add noise to the adc vector if we want that
      if (fAddNoise) this->AddNoiseToADCs(adcs);

      // check for zero suppression
      int retblocks = 1;
      gar::raw::Compress_t cflag = gar::raw::kNone;
      if (fCompressType == 1) {
        retblocks = gar::raw::Compress(
          adcs, raw::kZeroSuppression, fZSThreshold, fZSTicksBefore, fZSTicksAfter);
        cflag = gar::raw::kZeroSuppression;
      }

      todrop = (retblocks == 0);
      return raw::RawDigit(channel, numTicks, adcs, cflag, 0);
    }

    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
      fAddNoise = pset.get<bool>("AddNoise", false);
      fNoiseSpectrum = pset.get<int>("NoiseSpectrum", 0);
      fNoiseVecSize = pset.get<int>("NoiseVecSize", 500000);
      fNoiseAmplitude = pset.get<float>("NoiseAmplitude", 0);
      fCompressType = pset.get<int>("CompressType", 0);
      fZSThreshold = pset.get<int>("ZSThreshold", 5);
      fZSTicksBefore = pset.get<unsigned int>("ZSTicksBefore", 5);
      fZSTicksAfter = pset.get<unsigned int>("ZSTicksAfter", 5);
      fPedestal = pset.get<int>("Pedestal", 0);
      fADCSaturation = pset.get<int>("ADCSaturation", 4095); // default to 12-bit ADC's
      return;
    }

    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::CreateNoiseDigits(std::vector<raw::RawDigit>& digits)
    {
      if (!fAddNoise) return;

      // we'll need the geometry
      auto geo = gar::providerFrom<geo::GeometryGAr>();

      // figure out which channels we need to make pure noise
      std::set<unsigned int> channels;
      for (auto const& dig : digits)
        channels.insert(dig.Channel());

      // now loop over the channels in the detector and make a noise digit
      // for any channel that does not already have a rawdigit
      for (unsigned int c = 0; c < geo->NChannels(); ++c) {

        // do nothing if we already have a digit for this channel
        if (channels.count(c) > 0) continue;

        //std::cout << "Making noise digit for channel: " << c << std::endl;
        //std::cout << "ZS Threshold: " << fZSThreshold << std::endl;
        std::vector<short> adcs(fDetProp->NumberTimeSamples(), 0);
        this->AddNoiseToADCs(adcs);
        // check for zero suppression
        gar::raw::Compress_t cflag = gar::raw::kNone;
        int retblocks = 1;
        if (fCompressType == 1) {
          retblocks = gar::raw::Compress(
            adcs, raw::kZeroSuppression, fZSThreshold, fZSTicksBefore, fZSTicksAfter);
          cflag = gar::raw::kZeroSuppression;
        }
        //std::cout << "retblocks: " << retblocks << std::endl;
        auto numTicks = fDetProp->NumberTimeSamples();
        if (retblocks) digits.emplace_back(c, numTicks, adcs, cflag, 0);

      } // end make noise digits

      return;
    }

    // TODO:  Figure out how to parameterize noise that's correlated between channels
    // This method just works one channel at a time.
    // TODO: Include more noise models
    //----------------------------------------------------------------------------
    void TPCReadoutSimStandardAlg::AddNoiseToADCs(std::vector<short>& adcs)
    {
      // start sampling the pregenerated noise vector from a random spot
      CLHEP::RandFlat FlatRand(fEngine);
      double xnvi = ((double)fNoiseVec.size() - 1) * FlatRand.fire();
      if (xnvi < 0) xnvi = 0; // just to be extra sure
      size_t i = (size_t)xnvi;
      if (i >= fNoiseVec.size()) { i = fNoiseVec.size() - 1; }
      //std::cout << "In AddNoiseToADCs: " << i << " " << fNoiseVec.size() << " " << adcs.size() << std::endl;
      for (size_t j = 0; j < adcs.size(); ++j) {
        adcs[j] += fNoiseVec[i]; // we check the bounds of i already
        ++i;
        if (i >= fNoiseVec.size()) i = 0;
        if (adcs[j] > fADCSaturation) adcs[j] = fADCSaturation;
      }

      return;
    }

    //----------------------------------------------------------------------------
    short TPCReadoutSimStandardAlg::ElectronsToADCs(float electrons)
    {
      //LOG_DEBUG("TPCReadoutSimStandard")
      //<< "Electrons to ADC: "
      //<< fDetProp->ElectronsToADC()
      //<< " electrons: "
      //<< electrons
      //<< " product: "
      //<< fDetProp->ElectronsToADC() * electrons;

      int tmpadc = (int)(fDetProp->ElectronsToADC() * electrons);
      if (tmpadc > fADCSaturation) tmpadc = fADCSaturation;
      return (short)(tmpadc);
    }

  } // rosim
} // gar
