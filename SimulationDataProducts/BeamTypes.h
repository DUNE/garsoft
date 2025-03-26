

#ifndef BEAMTYPES_H
#define BEAMTYPES_H

namespace gar {
  namespace sdp {

    /// Defines category of beams to be stored in sdp::BeamGateInfo
    enum BeamType_t {
      kUnknown = 0, ///< Unknown beam type
      kBNB,         ///< BNB
      kNuMI,        ///< NuMI
      kBeamTypeMax  ///< Max value of enum for iteration
    };
  }
} // gar

#endif
