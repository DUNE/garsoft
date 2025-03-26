////////////////////////////////////////////////////////////////////////
///
///
/// \file  SimChannel.h
///
/// \brief object containing MC truth information necessary for making RawDigits
/// and doing back tracking
///
/// \author  brebel@fnal.gov
///
////////////////////////////////////////////////////////////////////////

#ifndef SIM_SIMCHANNEL_H
#define SIM_SIMCHANNEL_H

#include <set>
#include <stdint.h>
#include <string>
#include <vector>

namespace gar {
  namespace sdp {

    class IDE {
    public:
      IDE();

#ifndef __GCCXML__
      //constructor for IDEs applying G4 offset...
      IDE(IDE const&, int);
      IDE(int id, float numE, unsigned int channel, unsigned short tdc)
        : fTrackID(id), fNumElectrons(numE), fChannel(channel), fTDC(tdc)
      {}

      int const& TrackID() const { return fTrackID; }
      float const& NumElectrons() const { return fNumElectrons; }
      unsigned int const& Channel() const { return fChannel; }
      unsigned short const& TDC() const { return fTDC; }

      void operator+=(gar::sdp::IDE const& b);
      bool operator==(gar::sdp::IDE const& b) const;
      bool operator<(gar::sdp::IDE const& b) const;

#endif

      int fTrackID;          ///< Geant4 supplied track ID
      float fNumElectrons;   ///< total number of electrons for this track ID and time
      unsigned int fChannel; ///< channel number for these electrons
      unsigned short fTDC;   ///< TDC value of the readout
    };

  } // namespace sdp
} // gar

#endif // sdp_SIMCHANNEL_H

////////////////////////////////////////////////////////////////////////
