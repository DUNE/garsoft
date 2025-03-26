/// $Id: SimChannel.cxx,v 1.3 2010/03/26 20:08:36 brebel Exp $
///
/// \file  SimulationDataProducts/SimChannel.cxx
///
///
/// \author  seligman@nevis.columbia.edu
///
////////////////////////////////////////////////////////////////////////

#include <limits> // std::numeric_limits
#include <stdexcept>
#include <utility>

#include "SimulationDataProducts/SimChannel.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

//-------------------------------------------------
gar::sdp::IDE::IDE()
  : fTrackID(std::numeric_limits<int>::max())
  , fNumElectrons(std::numeric_limits<float>::max())
  , fChannel(std::numeric_limits<unsigned int>::max())
  , fTDC(std::numeric_limits<unsigned short>::max())
{}

//-------------------------------------------------
gar::sdp::IDE::IDE(sdp::IDE const& ide, int offset)
  : fTrackID(ide.TrackID() + offset)
  , fNumElectrons(ide.NumElectrons())
  , fChannel(ide.Channel())
  , fTDC(ide.TDC())
{}

//-------------------------------------------------
void gar::sdp::IDE::operator+=(gar::sdp::IDE const& b)
{
  if (fTrackID != b.TrackID() || fTDC != b.TDC()) {
    MF_LOG_WARNING("SimChannel") << "attempting to add IDEs with different trackIDs: " << fTrackID
                                 << "; " << b.TrackID() << " or TDCs " << fTDC << "; " << b.TDC()
                                 << " bail";
    return;
  }

  fNumElectrons += b.NumElectrons();

  return;
}

//-------------------------------------------------
bool gar::sdp::IDE::operator==(gar::sdp::IDE const& b) const
{
  return (fChannel == b.Channel() && std::abs(fTrackID) == std::abs(b.TrackID()) &&
          fTDC == b.TDC());
}

//-------------------------------------------------
// want to sort these by channel number, then TDC value, then
// TrackID contributing to charge at a given TDC
bool gar::sdp::IDE::operator<(gar::sdp::IDE const& b) const
{
  if (fChannel < b.Channel())
    return true;
  else if (fChannel == b.Channel()) {
    if (fTDC < b.TDC())
      return true;
    else if (fTDC == b.TDC()) {
      int absTrkID = std::abs(fTrackID);
      int absTrkIDB = std::abs(b.TrackID());

      if (absTrkID < absTrkIDB) return true;
    } // end if the same TDC value
  }   // end if the same channel

  return false;
}
