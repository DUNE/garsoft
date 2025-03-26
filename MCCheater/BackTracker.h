////////////////////////////////////////////////////////////////////////
/// \file  BackTracker.h
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GAR_CHEAT_BACKTRACKER_H
#define GAR_CHEAT_BACKTRACKER_H

#include <map>
#include <vector>

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "fhiclcpp/ParameterSet.h"

#include "MCCheater/BackTrackerCore.h"

namespace gar {
  namespace cheat {

    class BackTracker : public BackTrackerCore {

    public:
      using provider_type = BackTrackerCore; ///< type of service provider

      BackTracker(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
      ~BackTracker();

      void beginJob();

      /// Returns a pointer to the geometry service provider
      provider_type const* provider() const { return static_cast<provider_type const*>(this); }

      // The Rebuild function rebuilds the stl containers needed for backtracking queries.
      // It is called automatically at the time registered by the constructor.  In this case,
      // the constructor registers PreProcessEvent.  See the art::ActivityRegistry class.
      // art wants to see a ScheduleContext in the argument, but we don't need that.
      void Rebuild(art::Event const& evt, art::ScheduleContext);

    protected:
    private:
    };
  }
}

DECLARE_ART_SERVICE(gar::cheat::BackTracker, LEGACY)

#endif // CHEAT_BACKTRACK_H
