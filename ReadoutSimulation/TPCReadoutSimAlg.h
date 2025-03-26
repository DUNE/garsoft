//
//  TPCReadoutSimAlg.h
//
//  Created by Brian Rebel on 2/14/17.
//

#ifndef GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp
#define GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp

#include <vector>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"

#include "DetectorInfo/DetectorProperties.h"
#include "RawDataProducts/RawDigit.h"

#include <unordered_set>

namespace fhicl {
  class ParameterSet;
}

namespace gar {

  namespace rosim {

    // helper struct for tracking the information needed to make raw digits
    struct edepIDE {

      edepIDE();
      edepIDE(float numE, unsigned int chan, unsigned short tdc, size_t edepIdx, float edepWeight)
        : NumElect(numE), Channel(chan), TDC(tdc)
      {
        edepLocs.insert(edepIdx);
        edepWeights.push_front(edepWeight);
      }

      //------------------------------------------------------------------------
      bool operator<(edepIDE const& b) const
      {
        if (Channel < b.Channel)
          return true;
        else if (Channel == b.Channel) {
          if (TDC < b.TDC) return true;
        }

        return false;
      }

      //------------------------------------------------------------------------
      bool operator==(edepIDE const& b) const
      {
        if (Channel == b.Channel && TDC == b.TDC) return true;

        return false;
      }

      //------------------------------------------------------------------------
      bool operator!=(edepIDE const& b) const { return !(*this == b); }

      //------------------------------------------------------------------------
      void operator+=(edepIDE const& b)
      {
        if (Channel != b.Channel || TDC != b.TDC) {
          MF_LOG_WARNING("IonizationReadout") << "Attempting to add edepIDE with different "
                                              << "Channels: " << Channel << " / " << b.Channel
                                              << " or TDCs: " << TDC << " / " << b.TDC << " bail";
          return;
        }

        NumElect += b.NumElect;

        for (auto const e : b.edepLocs)
          edepLocs.insert(e);
        for (auto const w : b.edepWeights)
          edepWeights.push_front(w);

        return;
      }

      float NumElect;
      unsigned int Channel;
      unsigned short TDC;
      std::unordered_set<size_t> edepLocs; // keeps insert order, loads from front
      std::deque<float>
        edepWeights; // fraction of parent edep associated to digit, use deque for efficient front loading
    };

    class TPCReadoutSimAlg {

    public:
      TPCReadoutSimAlg(CLHEP::HepRandomEngine& engine, fhicl::ParameterSet const& pset);

      virtual ~TPCReadoutSimAlg();

      // Method to take IDEs and turn them into RawDigits
      virtual raw::RawDigit CreateRawDigit(unsigned int channel,
                                           std::vector<float> const& electrons,
                                           bool& todrop) = 0;
      virtual void CreateNoiseDigits(std::vector<raw::RawDigit>& digits) = 0;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

    protected:
      CLHEP::HepRandomEngine& fEngine;             ///< random number engine
      bool fAddNoise;                              ///< flag to add noise or not
      int fNoiseSpectrum;                          ///< 0: Gaussian white noise; more to come
      float fNoiseAmplitude;                       ///< noise amplitdue
      int fNoiseVecSize;                           ///< how much noise to pre-generate
      int fCompressType;                           ///< Switch to compress raw digits
      int fZSThreshold;                            ///< for ZS Compression, threshold (upwards)
      unsigned int fZSTicksBefore;                 ///< for ZS Compression, # samples before
      unsigned int fZSTicksAfter;                  ///< for ZS Compression, # samples after
      int fPedestal;                               ///< Raw Digit Pedestal
      const detinfo::DetectorProperties* fDetProp; ///< detector properties
      int fADCSaturation;                          ///< limit of the ADC

      // AddNoiseToADCs is foradding noise to recorded signal.  Assume that the
      // noise level is dependent on the amount of signal recorded
      virtual void AddNoiseToADCs(std::vector<short>& adcs) = 0;
      virtual short ElectronsToADCs(float electrons) = 0;
    };

  } // end rosim

} // end gar

#endif /* GAR_READOUTSIMULATION_TPCReadoutSimAlg_hpp */
