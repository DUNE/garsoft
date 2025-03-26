////////////////////////////////////////////////////////////////////////
///
/// \file  SimulationDataProducts/AuxDetSimChannel.cxx
///
////////////////////////////////////////////////////////////////////////

// our header
#include "SimulationDataProducts/AuxDetSimChannel.h"

// C/C++ standard library
#include <limits> // std::numeric_limits<>
#include <stdexcept>

namespace gar {
  namespace sdp {

    // Default constructor
    //-------------------------------------------------
    AuxDetIDE::AuxDetIDE()
      : trackID(std::numeric_limits<int>::min())
      , energyDeposited(std::numeric_limits<float>::min())
      , entryX(std::numeric_limits<float>::min())
      , entryY(std::numeric_limits<float>::min())
      , entryZ(std::numeric_limits<float>::min())
      , entryT(std::numeric_limits<float>::min())
      , exitX(std::numeric_limits<float>::min())
      , exitY(std::numeric_limits<float>::min())
      , exitZ(std::numeric_limits<float>::min())
      , exitT(std::numeric_limits<float>::min())
      , exitMomentumX(std::numeric_limits<float>::min())
      , exitMomentumY(std::numeric_limits<float>::min())
      , exitMomentumZ(std::numeric_limits<float>::min())
    {}

    // Copy with offset constructor
    //-------------------------------------------------
    AuxDetIDE::AuxDetIDE(AuxDetIDE const& ide, int offset)
      : trackID(ide.trackID + offset)
      , energyDeposited(ide.energyDeposited)
      , entryX(ide.entryX)
      , entryY(ide.entryY)
      , entryZ(ide.entryZ)
      , entryT(ide.entryT)
      , exitX(ide.exitX)
      , exitY(ide.exitY)
      , exitZ(ide.exitZ)
      , exitT(ide.exitT)
      , exitMomentumX(ide.exitMomentumX)
      , exitMomentumY(ide.exitMomentumY)
      , exitMomentumZ(ide.exitMomentumZ)
    {}

    //----------------------------------------------------------------------------
    AuxDetSimChannel::AuxDetSimChannel()
      : fAuxDetID(std::numeric_limits<uint32_t>::max())
      , fAuxDetSensitiveID(std::numeric_limits<uint32_t>::max())
    {}

    //----------------------------------------------------------------------------
    AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID, uint32_t inputAuxDetSensitiveID)
      : fAuxDetID(inputAuxDetID), fAuxDetSensitiveID(inputAuxDetSensitiveID)
    {}

    //----------------------------------------------------------------------------
    AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID,
                                       std::vector<sdp::AuxDetIDE> const& inputAuxDetIDEs,
                                       uint32_t inputAuxDetSensitiveID)
      : fAuxDetID(inputAuxDetID)
      , fAuxDetSensitiveID(inputAuxDetSensitiveID)
      , fAuxDetIDEs(inputAuxDetIDEs)
    {}

    //----------------------------------------------------------------------------
    AuxDetSimChannel::AuxDetSimChannel(uint32_t inputAuxDetID,
                                       std::vector<sdp::AuxDetIDE>&& inputAuxDetIDEs,
                                       uint32_t inputAuxDetSensitiveID)
      : fAuxDetID(inputAuxDetID)
      , fAuxDetSensitiveID(inputAuxDetSensitiveID)
      , fAuxDetIDEs(inputAuxDetIDEs)
    {}

    //----------------------------------------------------------------------------
    std::pair<int, int> AuxDetSimChannel::MergeAuxDetSimChannel(const AuxDetSimChannel& chan,
                                                                int offset)
    {
      if (this->fAuxDetID != chan.AuxDetID() &&
          this->fAuxDetSensitiveID != chan.AuxDetSensitiveID())
        throw std::runtime_error(
          "ERROR AuxDetSimChannel Merge: Trying to merge different channels!");

      std::pair<int, int> range_trackID(std::numeric_limits<int>::max(),
                                        std::numeric_limits<int>::min());

      for (auto const& ide : AuxDetIDEs()) {
        this->fAuxDetIDEs.emplace_back(ide, offset);

        if (ide.trackID + offset < range_trackID.first) range_trackID.first = ide.trackID + offset;
        if (ide.trackID + offset > range_trackID.second)
          range_trackID.second = ide.trackID + offset;
      }

      return range_trackID;
    }

    //----------------------------------------------------------------------------
    bool AuxDetSimChannel::operator<(const sdp::AuxDetSimChannel& other) const
    {
      if (fAuxDetID < other.AuxDetID()) return true;

      return fAuxDetSensitiveID < other.AuxDetSensitiveID();
    }

    //----------------------------------------------------------------------------
    bool AuxDetSimChannel::operator==(const sdp::AuxDetSimChannel& other) const
    {
      return (fAuxDetID == other.AuxDetID() && fAuxDetSensitiveID == other.AuxDetSensitiveID());
    }

  } //namespace sim
} // gar