
#ifndef CLOCKCONSTANTS_H
#define CLOCKCONSTANTS_H

#include <climits>
#include <limits>

namespace gar {
  namespace detinfo {

    /// Maximum time in nano-second
    const double kTIME_MAX = std::numeric_limits<double>::max();

    /// Default Frequency in MHz
    const double kDEFAULT_FREQUENCY = 1.e-6;

    /// Default Frame period in nano-second
    const double kDEFAULT_FRAME_PERIOD = 1.6e6;

    /// Default Optical clock speed in MHz
    const double kDEFAULT_FREQUENCY_OPTICAL = 64.;

    /// Default TPC clock speed in MHz
    const double kDEFAULT_FREQUENCY_TPC = 2.;

    /// Default Trigger clock speed in MHz
    const double kDEFAULT_FREQUENCY_TRIGGER = 16.;

    /// Default External clock speed in MHz
    const double kDEFAULT_FREQUENCY_EXTERNAL = 31.25;

    /// Default G4 reference time in nano-second
    const double kDEFAULT_MC_CLOCK_T0 = 0;

    /// Default TPC readout start time offset from trigger in nano-second
    const double kDEFAULT_TRIG_OFFSET_TPC = 0.;

  }
} // gar

#endif
