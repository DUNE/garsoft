////////////////////////////////////////////////////////////////////////
/// \file EvdLayoutOptions_service.cc
///
/// \version $Id: EvdLayoutOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  andrzejs@fnal.gov

// Framework includes

/// GArSoft includes
#include "EventDisplay/EVD/EvdLayoutOptions.h"

#include <iostream>

namespace gar {
  namespace evd {

    //......................................................................
    EvdLayoutOptions::EvdLayoutOptions(fhicl::ParameterSet const& pset,
                                       art::ActivityRegistry& /* reg */)
    {
      this->reconfigure(pset);
    }

    //......................................................................
    EvdLayoutOptions::~EvdLayoutOptions() {}

    //......................................................................
    void EvdLayoutOptions::reconfigure(fhicl::ParameterSet const& pset)
    {
      fEnableMCTruthCheckBox = pset.get<int>("EnableMCTruthCheckBox");
    }
  }

  namespace evd {

    DEFINE_ART_SERVICE(EvdLayoutOptions)

  } // namespace evd
}
////////////////////////////////////////////////////////////////////////
