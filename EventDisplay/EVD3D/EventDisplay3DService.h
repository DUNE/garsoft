//
// - Event display service used with the EventDisplay3D module to start up the
//   TApplication and allow forward, backward, and jump-to navigation of events
//   in the root input file.  This is a very much simplified version based on
//   Nova's event visplay service by Mark Messier.
//

#ifndef EventDisplay3DService_EventDisplay3DService_hh
#define EventDisplay3DService_EventDisplay3DService_hh

#ifndef __CINT__
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Provenance/ScheduleContext.h"

namespace gar
{
  namespace evd3d
  {
    class EventDisplay3DService
    {
    public:

      EventDisplay3DService(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

    private:

      void postBeginJobWorkers(art::InputSource* inputs,
        std::vector<art::Worker*> const& workers);
      void postProcessEvent(art::Event const&, art::ScheduleContext);

      private:
        art::InputSource* fInputSource; ///< Input source of events

      public:
      };
    }
  }

  #endif // __CINT__
  DECLARE_ART_SERVICE(gar::evd3d::EventDisplay3DService, LEGACY)
  #endif // EventDisplay3DService_EventDisplay3DService_hh
