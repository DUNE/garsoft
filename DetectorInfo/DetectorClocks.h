////////////////////////////////////////////////////////////////////////
// \file
//
// \brief pure virtual base interface for detector clocks
//
// \author jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_IDETCLOCKS_H
#define DETINFO_IDETCLOCKS_H

#include "DetectorInfo/ElecClock.h"

namespace gar {
  namespace detinfo {

    enum InheritConfigType_t {
      kG4RefTime = 0,
      kTriggerOffsetTPC,
      kFramePeriod,
      kClockSpeedTPC,
      kClockSpeedTrigger,
      kClockSpeedExternal,
      kDefaultTrigTime,
      kDefaultBeamTime,
      kDefaultSpillLength,
      kInheritConfigTypeMax
    };

    class DetectorClocks {

    public:
      DetectorClocks(const DetectorClocks&) = delete;
      DetectorClocks(DetectorClocks&&) = delete;
      DetectorClocks& operator=(const DetectorClocks&) = delete;
      DetectorClocks& operator=(DetectorClocks&&) = delete;
      virtual ~DetectorClocks() = default;

      virtual double TriggerOffsetTPC() const = 0;

      /// Given Geant4 time [ns], returns relative time [ns] w.r.t. electronics time T0
      virtual double G4ToElecTime(double g4_time) const = 0;

      /// Trigger electronics clock time in [ns]
      virtual double TriggerTime() const = 0;

      /// Beam gate electronics clock time in [ns]
      virtual double BeamGateTime() const = 0;

      /// Duration of spill [ns]
      virtual double SpillLength() const = 0;

      virtual std::vector<std::string> ConfigNames() const = 0;
      virtual std::vector<double> ConfigValues() const = 0;

      //
      // Getters of TPC ElecClock
      //
      /// Borrow a const TPC clock with time set to Trigger time [ns]
      virtual const ElecClock& TPCClock() const = 0;

      /// Create a TPC clock for a given time [ns] from clock counting start
      virtual ElecClock TPCClock(double time) const = 0;

      /// Create a TPC clock for a given sample/frame number in TPC clock frequency
      virtual ElecClock TPCClock(unsigned int sample, unsigned int frame) const = 0;

      //
      // Getters of Trigger ElecClock
      //
      /// Borrow a const Trigger clock with time set to Trigger time [ns]
      virtual const detinfo::ElecClock& TriggerClock() const = 0;

      /// Create a Trigger clock for a given time [ns] from clock counting start
      virtual detinfo::ElecClock TriggerClock(double time) const = 0;

      /// Create a Trigger clock for a given sample/frame number in Trigger clock frequency
      virtual detinfo::ElecClock TriggerClock(unsigned int sample, unsigned int frame) const = 0;

      //
      // Getters of External ElecClock
      //
      /// Borrow a const Trigger clock with time set to External Time [ns]
      virtual const detinfo::ElecClock& ExternalClock() const = 0;

      /// Create a External clock for a given time [ns] from clock counting start
      virtual detinfo::ElecClock ExternalClock(double time) const = 0;

      /// Create a External clock for a given sample/frame number in External clock frequency
      virtual detinfo::ElecClock ExternalClock(unsigned int sample, unsigned int frame) const = 0;

      //
      // Getters for time [us] w.r.t. trigger given information from waveform
      //

      /// Given TPC time-tick (waveform index), returns time [ns] w.r.t. trigger time stamp
      virtual double TPCTick2TrigTime(double tick) const = 0;

      /// Given TPC time-tick (waveform index), returns time [ns] w.r.t. beam gate time
      virtual double TPCTick2BeamTime(double tick) const = 0;

      /// Given External time-tick (waveform index), sample and frame number, returns time [ns] w.r.t. trigger time stamp
      virtual double ExternalTick2TrigTime(double tick, size_t sample, size_t frame) const = 0;

      /// Given External time-tick (waveform index), sample and frame number, returns time [ns] w.r.t. beam gate time stamp
      virtual double ExternalTick2BeamTime(double tick, size_t sample, size_t frame) const = 0;

      //
      // Getters for time [tdc] (electronics clock counting ... in double precision)
      //

      /// Given TPC time-tick (waveform index), returns electronics clock count [tdc]
      virtual double TPCTick2TDC(double tick) const = 0;

      /// Given G4 time [ns], returns corresponding TPC electronics clock count [tdc]
      virtual double TPCG4Time2TDC(double g4time) const = 0;

      /// Given External time-tick (waveform index), sample and frame number, returns time electronics clock count [tdc]
      virtual double ExternalTick2TDC(double tick, size_t sample, size_t frame) const = 0;

      /// Given G4 time [ns], returns corresponding External electronics clock count [tdc]
      virtual double ExternalG4Time2TDC(double g4time) const = 0;

      //
      // Getters for time [ns] (electronics clock counting ... in double precision)
      //
      /// Given TPC time-tick (waveform index), returns electronics clock [ns]
      virtual double TPCTick2Time(double tick) const = 0;

      /// Given External time-tick (waveform index), sample and frame number, returns electronics clock [ns]
      virtual double ExternalTick2Time(double tick, size_t sample, size_t frame) const = 0;

      //
      // Getters for time [ticks] (waveform index number)
      //

      /// Given electronics clock count [tdc] returns TPC time-tick
      virtual double TPCTDC2Tick(double tdc) const = 0;

      /// Given G4 time returns electronics clock count [tdc]
      virtual double TPCG4Time2Tick(double g4time) const = 0;

    protected:
      DetectorClocks() = default;

    }; // class DetectorClocks

  } //namespace detinfo
} // gar

#endif
