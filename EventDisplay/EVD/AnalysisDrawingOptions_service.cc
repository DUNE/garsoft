////////////////////////////////////////////////////////////////////////
/// \file AnalysisDrawingOption_plugin.cc
///
/// \version $Id: AnalysisDrawingOptions_plugin.cc,v 1.1 2010/11/11 18:11:22 p-novaart Exp $
/// \author  brebel@fnal.gov

// Framework includes

/// GArSoft includes
#include "EventDisplay/EVD/AnalysisDrawingOptions.h"

#include <iostream>

namespace gar {
  namespace evd {

    //......................................................................
    AnalysisDrawingOptions::AnalysisDrawingOptions(fhicl::ParameterSet const& pset,
                                                   art::ActivityRegistry& /* reg */)
    {
      this->reconfigure(pset);
    }

    //......................................................................
    AnalysisDrawingOptions::~AnalysisDrawingOptions() {}

    //......................................................................
    void AnalysisDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
    {
      fDrawCalorimetry = pset.get<int>("DrawCalorimetry");
      fDrawParticleID = pset.get<int>("DrawParticleID");
      fDrawShowerCalor = pset.get<int>("DrawShowerCalor");
      fCaloPlane = pset.get<int>("CaloPlane");
      fCalorimetryLabels = pset.get<std::vector<std::string>>("CalorimetryModuleLabels");
      fParticleIDLabels = pset.get<std::vector<std::string>>("ParticleIDModuleLabels");

      fCalorTemplateFileName = pset.get<std::string>("CalorTemplateFileName");
    }

  }

  namespace evd {

    DEFINE_ART_SERVICE(AnalysisDrawingOptions)

  } // namespace evd
}
////////////////////////////////////////////////////////////////////////
