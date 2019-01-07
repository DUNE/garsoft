////////////////////////////////////////////////////////////////////////
/// \file  BackTracker.h
/// \brief back track the reconstruction to the simulation
///
/// \version $Id: Geometry.h,v 1.16 2009/11/03 22:53:20 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GAR_CHEAT_BACKTRACKER_H
#define GAR_CHEAT_BACKTRACKER_H

#include <vector>
#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Persistency/Provenance/ScheduleContext.h"

#include "MCCheater/BackTrackerCore.h"


namespace gar{
  namespace cheat{
    
    class BackTracker: public BackTrackerCore{
      
    public:
      
      using provider_type = BackTrackerCore; ///< type of service provider
      
      BackTracker(fhicl::ParameterSet     const& pset,
                  ::art::ActivityRegistry      & reg);
      ~BackTracker();
      
      void beginJob();
      
      /// Returns a pointer to the geometry service provider
      provider_type const* provider() const { return static_cast<provider_type const*>(this); }

      // The Rebuild function rebuilds the various maps we need to answer backtracking queries.
      // It is called automatically before each event is processed. For jobs involving
      // Monte Carlo generation, this is too soon. So, we'll call rebuild after those data
      // products are put into the event in LArG4.  This is the least bad way of ensuring the
      // BackTracker works in jobs that combine MC production and reconstruction analysis based
      // on MC truth.  Don't copy this design pattern without talking to brebel@fnal.gov first
      void Rebuild(::art::Event const& evt, art::ScheduleContext);
      void RebuildNoSC(::art::Event const& evt);
      
    private:
      
    };
  } // namespace
}

DECLARE_ART_SERVICE(gar::cheat::BackTracker, LEGACY)

#endif // CHEAT_BACKTRACK_H
