////////////////////////////////////////////////////////////////////////
// $Id: EvdLayoutOptions.h,v 1.15 2010/08/30 21:33:24 spitz7 Exp $
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVDLAYOUTOPTIONS_H
#define EVDLAYOUTOPTIONS_H
#ifndef __CINT__
#include <string>
#include <vector>

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar {
  namespace evd {
    class EvdLayoutOptions {
    public:
      EvdLayoutOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~EvdLayoutOptions();

      void reconfigure(fhicl::ParameterSet const& pset);

      int fShowSideBar;           ///< 1 to show, 0 don't show
      int fAutoZoomInterest;      ///< Set the automatic zoom to the interest region
      int fPrintTotalCharge;      ///< Print out the total charge in an event
      int fShowEndPointSection;   ///< Show section corresponding to EndPoint finding
      int fShowEndPointMarkers;   ///< Draw EndPoint Markers if clicked.
      int fShowClusterSection;    ///< Show section to make clusters
      int fMakeClusters;          ///< Draw two lines to make clusters if clicked
      int fMakeSeeds;             ///< Draw two lines to make clusters if clicked
      int fEnableMCTruthCheckBox; ///< 1 to have the check box appear, 0 otherwise
    };
  }
} //namespace
#endif // __CINT__
DECLARE_ART_SERVICE(gar::evd::EvdLayoutOptions, LEGACY)
#endif
