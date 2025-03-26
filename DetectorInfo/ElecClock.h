/**
 * \file ElecClock.hh
 *
 * \ingroup DetectorClocksService
 *
 * \brief Class def header for a class ElecClock
 *
 * @author kterao
 */

/** \addtogroup TimeService

    @{*/
#ifndef ElecClock_H
#define ElecClock_H

#include "DetectorInfo/ClockConstants.h"
#include "DetectorInfo/DetectorClocksException.h"
#include <iostream>

namespace gar {
  namespace detinfo {
    /**
     \class ElecClock
     Electronics clock class which holds double precision time in nanosecond,
     frequency in Hz, and frame period in seconds.
     */
    class ElecClock {

    public:
      /// Default constructor
      ElecClock(double time = 0, double frame_period = kTIME_MAX, double frequency = 1e9)
        : fTime(time), fFramePeriod(frame_period), fFrequency(frequency)
      {
        if (fFrequency <= 0)
          throw detinfo::DetectorClocksException("Negative frequency is prohibited!");
      }

      /// Default destructor
      ~ElecClock(){};

    protected:
      //-- Attribute variables --//

      double fTime;        ///< Time in nano-second
      double fFramePeriod; ///< Frame period in micro-second
      double fFrequency;   ///< Clock speed in MHz

    public:
      //-- Setters --//

      void SetTime(double time) { fTime = time; }

      void SetTime(int sample, int frame) { fTime = Time(sample, frame); }

      void SetTime(unsigned int sample, unsigned int frame) { SetTime(int(sample), int(frame)); }

      void SetTime(int ticks) { fTime = Time(ticks, 0); }

      void SetTime(unsigned int ticks) { SetTime(int(ticks)); }

      //-- Getters --//

      /// Time (stored) in nano-second
      double Time() const { return fTime; }

      /// Given sample & frame number in *this* clock, returns double precision time [ns]
      double Time(int sample, int frame) const
      {
        return (sample * 1.e3 / fFrequency + frame * fFramePeriod);
      }

      /// Given time [ns] w.r.t. electronics clock counting, return discretized time in double precision
      double Time(double time) const { return Time(Sample(time), Frame(time)); }

      /// Given time in ticks w.r.t. electronics clock counting, return discretized time in double precision
      double Time(int ticks) const { return ticks * 1.e3 / fFrequency; }

      /// Frequency in MHz
      double Frequency() const { return fFrequency; }

      /// A single frame period in micro-second
      double FramePeriod() const { return fFramePeriod; }

      /// # of Ticks
      int Ticks() const { return Ticks(fTime); }

      /// Given time [ns] w.r.t. electronics clock T0, return # of ticks
      int Ticks(double time) const { return (int)(time * 1.e-3 * fFrequency); }

      /// Given sample & frame returns # ticks
      int Ticks(int sample, int frame) const { return sample + frame * FrameTicks(); }

      /// Sample number
      int Sample() const { return Sample(fTime); }

      /// Given time [ns] w.r.t. electronics clock T0, return sample number
      int Sample(double time) const
      {
        return (int)((time - Frame(time) * fFramePeriod) * 1.e-3 * fFrequency);
      }

      /// Given ticks w.r.t. electronics clock T0, return sample number
      int Sample(int tick) const { return (tick % (int)(FrameTicks())); }

      /// Frame number
      int Frame() const { return Frame(fTime); }

      /// Given time [ns] w.r.t. electronics clock T0, return frame number
      int Frame(double time) const { return (int)(time / fFramePeriod); }

      /// Given ticks w.r.t. electronics clock T0, return frame number
      int Frame(int tick) const { return (tick / (int)(FrameTicks())); }

      /// Number ticks in a frame
      unsigned int FrameTicks() const { return (int)(fFramePeriod * 1.e-3 * fFrequency); }

      /// A single tick period in nano-second, frequency is in MHz
      double TickPeriod() const { return 1.e3 / fFrequency; }

    public:
      //-- operators --//
      // Remove these for now as they do not appear to be used.
      // If they are used again they will have to be updated so that the units
      // in the calculations are consistently ns
      //      ElecClock& operator++()                        { fTime += 1./fFrequency; return *this;}
      //      ElecClock  operator++(int)                     {ElecClock tmp(*this); operator++(); return tmp;}
      //      ElecClock& operator--()                        { fTime -= 1./fFrequency; return *this;}
      //      ElecClock  operator--(int)                     {ElecClock tmp(*this); operator--(); return tmp;}
      //      ElecClock& operator+=(const double &rhs)       { fTime += rhs; return *this;}
      //      ElecClock& operator-=(const double &rhs)       { fTime -= rhs; return *this;}
      //      ElecClock& operator+=(const float  &rhs)       { fTime += (double)rhs; return *this;}
      //      ElecClock& operator-=(const float  &rhs)       { fTime -= (double)rhs; return *this;}
      //      ElecClock& operator+=(const int &rhs)          { fTime += Time(rhs); return *this;}
      //      ElecClock& operator-=(const int &rhs)          { fTime -= Time(rhs); return *this;}
      //      ElecClock& operator+=(const unsigned int &rhs) { fTime += Time((int)rhs); return *this;}
      //      ElecClock& operator-=(const unsigned int &rhs) { fTime -= Time((int)rhs); return *this;}
      //      ElecClock& operator+=(const ElecClock& rhs)    { fTime += rhs.Time(); return *this;}
      //      ElecClock& operator-=(const ElecClock& rhs)    { fTime -= rhs.Time(); return *this;}

      inline ElecClock operator+(const ElecClock& rhs)
      {
        return ElecClock(fTime + rhs.Time(), fFramePeriod, fFrequency);
      }

      inline ElecClock operator-(const ElecClock& rhs)
      {
        return ElecClock(fTime - rhs.Time(), fFramePeriod, fFrequency);
      }

      inline bool operator<(const ElecClock& rhs) const { return fTime < rhs.Time(); }
      inline bool operator>(const ElecClock& rhs) const { return fTime > rhs.Time(); }
      inline bool operator<=(const ElecClock& rhs) const { return fTime <= rhs.Time(); }
      inline bool operator>=(const ElecClock& rhs) const { return fTime >= rhs.Time(); }
    };

  }
} // gar

#endif
/** @} */ // end of doxygen group
