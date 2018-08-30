////////////////////////////////////////////////////////////////////////
// $Id: RawDrawingOption.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef RAWDRAWINGOPTIONS_H
#define RAWDRAWINGOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "canvas/Utilities/InputTag.h"

namespace gar {
namespace evd {
  /**
   * @brief Display parameters for the raw data
   *
   * Configuration parameters
   * -------------------------
   *
   * This is an incomplete list of the supported parameters:
   * - *RoIthresholds* (list of real numbers, default: empty): threshold in ADC
   *   counts for a tick on a wire to be "interesting", thus extending the
   *   region of interest to include it. The thresholds are specified one per
   *   plane; if no threshold is specified for a plane, the threshold for the
   *   last plane is used instead (therefore specifying just one threshold will
   *   apply the same threshold to all planes). If no threshold is specified
   *   at all, the value of 'MinSignal' parameter is used as threshold for all
   *   planes
   *
   */
  class RawDrawingOptions
  {
  public:
    RawDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~RawDrawingOptions();

    void reconfigure(fhicl::ParameterSet const& pset);

    int                fDrawRawOrReco;            ///< 0 = draw raw, 1 = draw reco
    int    	           fTicksPerPoint;            ///< number of ticks to include in one point
    int    	           fScaleDigitsByCharge;      ///< scale the size of the digit by the charge
    double 	           fMinSignal;                ///< minimum ADC count to display a time bin
    double             fStartTick;                ///< Starting tick for the display
    double 	           fTicks;                    ///< number of TDC ticks to display, ie # fTicks past fStartTick
    unsigned int       fMinChannelStatus;         ///< Display channels with this status and above
    unsigned int       fMaxChannelStatus;         ///< Display channels with this status and below
    unsigned int       fChannel;                  ///< Channel to display in time/charge histogram
    art::InputTag      fRawDataLabel;             ///< module label that made the raw digits, default is daq

    bool               fUncompressWithPed;        ///< Option to uncompress with pedestal. Turned off by default
    bool               fSeeBadChannels;           ///< Allow "bad" channels to be viewed

  };
}
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(gar::evd::RawDrawingOptions, LEGACY)
#endif
