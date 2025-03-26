////////////////////////////////////////////////////////////////////////
/// \file  AuxDetSimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits
/// and doing back tracking
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef SDP_AUXDETSIMCHANNEL_H
#define SDP_AUXDETSIMCHANNEL_H

// C/C++ standard libraries
#include <stdint.h> // C header (need to be compatible with Reflex)
#include <vector>

namespace gar {
  namespace sdp {

    /**
     * @brief MC truth information to make RawDigits and do back tracking
     *
     * This structure describes the true position of momentum of a MC particle
     * entering and exiting a scintillator cell (channel) of an auxiliary
     * scintillator detector.
     */
    class AuxDetIDE {

    public:
      AuxDetIDE();

      //constructor for IDEs applying G4 offset...
      AuxDetIDE(AuxDetIDE const&, int);

      int trackID;           ///< Geant4 supplied track ID
      float energyDeposited; ///< total energy deposited for this track ID and time
      float entryX;          ///< Entry position X of particle
      float entryY;          ///< Entry position Y of particle
      float entryZ;          ///< Entry position Z of particle
      float entryT;          ///< Entry time of particle
      float exitX;           ///< Exit position X of particle
      float exitY;           ///< Exit position Y of particle
      float exitZ;           ///< Exit position Z of particle
      float exitT;           ///< Exit time of particle
      float exitMomentumX;   ///< Exit X-Momentum of particle
      float exitMomentumY;   ///< Exit Y-Momentum of particle
      float exitMomentumZ;   ///< Exit Z-Momentum of particle

#ifndef __GCCXML__
      bool operator<(const AuxDetIDE& other) const;
      bool operator==(const AuxDetIDE& other) const;

#endif
    }; // class AuxDetIDE

    /**
     * @brief Collection of particles crossing one auxiliary detector cell
     *
     * This structure collects information (as sdp::AuxDetIDE) from all the MC
     * particles crossing a single auxiliary detector cell (channel).
     */
    class AuxDetSimChannel {

    public:
      /// Default constructor (invalid, empty data)
      AuxDetSimChannel();

    private:
      uint32_t fAuxDetID; ///< geo->AuxDet(auxDetID), integer used to retrieve AuxDetGeo objec
      uint32_t fAuxDetSensitiveID; ///< integer used to retrieve AuxDetSensitiveGeo object
      std::vector<sdp::AuxDetIDE> fAuxDetIDEs; ///< one sdp::AuxDetIDE for each G4 track id

#ifndef __GCCXML__
    public:
      AuxDetSimChannel(uint32_t inputAuxDetID, uint32_t inputAuxDetSensitiveID);

      /// Constructor: copies from the specified IDE vector
      AuxDetSimChannel(uint32_t inputAuxDetID,
                       const std::vector<sdp::AuxDetIDE>& inputAuxDetIDEs,
                       uint32_t inputAuxDetSensitiveID = 0);

      /// Constructor: moves data from the specified IDE vector
      AuxDetSimChannel(uint32_t inputAuxDetID,
                       std::vector<sdp::AuxDetIDE>&& inputAuxDetIDEs,
                       uint32_t inputAuxDetSensitiveID = 0);

      std::pair<int, int> MergeAuxDetSimChannel(const AuxDetSimChannel&, int);

      ///@name Getters
      ///@{
      uint32_t AuxDetID() const;
      uint32_t AuxDetSensitiveID() const;

      bool operator<(const AuxDetSimChannel& other) const;
      bool operator==(const AuxDetSimChannel& other) const;

      std::vector<sdp::AuxDetIDE> const& AuxDetIDEs() const;
      ///@}

#endif

    }; // class AuxDetSimChannel

  } // namespace sim
} // gar

#ifndef __GCCXML__

inline bool gar::sdp::AuxDetIDE::operator<(const AuxDetIDE& other) const
{
  return trackID < other.trackID;
}
inline bool gar::sdp::AuxDetIDE::operator==(const AuxDetIDE& other) const
{
  return other.trackID == trackID;
}
inline uint32_t gar::sdp::AuxDetSimChannel::AuxDetID() const
{
  return fAuxDetID;
}
inline uint32_t gar::sdp::AuxDetSimChannel::AuxDetSensitiveID() const
{
  return fAuxDetSensitiveID;
}
inline std::vector<gar::sdp::AuxDetIDE> const& gar::sdp::AuxDetSimChannel::AuxDetIDEs() const
{
  return fAuxDetIDEs;
}
#endif

#endif // SDP_AUXDETSIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
