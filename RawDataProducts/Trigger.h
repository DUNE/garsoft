#ifndef TRIGGER_H
#define TRIGGER_H

// C++ includes
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace gar {
  namespace raw {

    class Trigger {

    public:
      /// Default ctor
      Trigger()
      {
        fTriggerNumber = std::numeric_limits<unsigned int>::max();

        fTriggerTime = std::numeric_limits<double>::max();
        fBeamGateTime = std::numeric_limits<double>::max();

        fTriggerBits = 0x0;
      }

    private:
      unsigned int fTriggerNumber; ///< Trigger counter
      double fTriggerTime;         ///< Trigger time w.r.t. electronics clock T0
      double fBeamGateTime;        ///< BeamGate time w.r.t. electronics clock T0
      unsigned int fTriggerBits;   ///< Trigger bits ... dedicated bit-by-bit function available

    public:
#ifndef __GCCXML__

      /// Alternative constructor
      Trigger(unsigned int counter, double trigger_time, double beamgate_time, uint32_t bits)
        : fTriggerNumber(counter)
        , fTriggerTime(trigger_time)
        , fBeamGateTime(beamgate_time)
        , fTriggerBits(bits)
      {}

      /// Trigger number
      unsigned int TriggerNumber() const { return fTriggerNumber; }
      /// Trigger time w.r.t. electronics clock T0 in us
      double TriggerTime() const { return fTriggerTime; }
      /// BeamGate time w.r.t. electronics clock T0 in us
      double BeamGateTime() const { return fBeamGateTime; }
      /// Trigger Bits
      unsigned int TriggerBits() const { return fTriggerBits; }
      /// Accessor to specific bit
      bool Triggered(const unsigned char bit) const;

#endif
    };
  }
} // gar

#endif
