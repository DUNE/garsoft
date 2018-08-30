////////////////////////////////////////////////////////////////////////
/// \file RawDrawingOption_plugin.cc
///
/// \version $Id: RecoDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// GArSoft includes
#include "EventDisplay/RawDrawingOptions.h"

#include <iostream>
#include <limits>

namespace gar {
namespace evd {

  //......................................................................
  RawDrawingOptions::RawDrawingOptions(fhicl::ParameterSet const& pset,
                                       art::ActivityRegistry& /* reg */)
  {
    this->reconfigure(pset);
  }

  //......................................................................
  RawDrawingOptions::~RawDrawingOptions()
  {
  }

  //......................................................................
  void RawDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fDrawRawOrReco            = pset.get< int          >("DrawRawOrReco",        0    );
    fScaleDigitsByCharge     	= pset.get< int          >("ScaleDigitsByCharge"        );
    fTicksPerPoint            = pset.get< int          >("TicksPerPoint"              );
    fMinSignal                = pset.get< double       >("MinimumSignal"              );
    fStartTick                = pset.get< double       >("StartTick",            0    );
    fTicks                    = pset.get< double       >("TotalTicks",           2048 );
    fRawDataLabel             = pset.get< art::InputTag>("RawDataLabel",        "daq");
    fMinChannelStatus         = pset.get< unsigned int >("MinChannelStatus",     0    );
    fMaxChannelStatus         = pset.get< unsigned int >("MaxChannelStatus",     std::numeric_limits<unsigned int>::max());
    fUncompressWithPed        = pset.get< bool         >("UncompressWithPed",    false);
    fSeeBadChannels           = pset.get< bool         >("SeeBadChannels",       false);
    fChannel                  = pset.get< unsigned int >("Channel",              0);

  }
}

namespace evd {

  DEFINE_ART_SERVICE(RawDrawingOptions)

} // namespace evd
}
////////////////////////////////////////////////////////////////////////
