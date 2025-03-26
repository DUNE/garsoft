#ifndef ANALYSISDRAWINGOPTIONS_H
#define ANALYSISDRAWINGOPTIONS_H
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
    class AnalysisDrawingOptions {
    public:
      AnalysisDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~AnalysisDrawingOptions();

      void reconfigure(fhicl::ParameterSet const& pset);

      int fDrawCalorimetry;
      int fDrawParticleID;
      int fDrawShowerCalor;
      int fCaloPlane;

      std::vector<std::string> fCalorimetryLabels; ///< module labels that produced calorimetry
      std::vector<std::string> fParticleIDLabels;  ///< module labels that produced particleid
      std::string fCalorTemplateFileName;          ///< files that have calorimetry template curves
    };
  }
} //namespace
#endif // __CINT__
DECLARE_ART_SERVICE(gar::evd::AnalysisDrawingOptions, LEGACY)
#endif
