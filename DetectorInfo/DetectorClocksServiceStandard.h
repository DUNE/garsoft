////////////////////////////////////////////////////////////////////////
// DetectorClocksServiceStandard.h
//
// Service interface for Detector Clock functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORCLOCKSSERVICESTANDARD_H
#define DETECTORCLOCKSSERVICESTANDARD_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Provenance/ScheduleContext.h"

#include "DetectorInfo/DetectorClocksStandard.h"
#include "DetectorInfo/DetectorClocksService.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    class DetectorClocksServiceStandard : public DetectorClocksService {
    public:
      DetectorClocksServiceStandard(fhicl::ParameterSet   const& pset,
                                    ::art::ActivityRegistry      & reg);
      
      virtual void  reconfigure(fhicl::ParameterSet const& pset) override;
      void          preBeginRun(::art::Run const& run);
      void          preProcessEvent(::art::Event const& evt, art::ScheduleContext);
      void          postOpenFile(std::string const& filename);
      
      virtual const provider_type* provider() const override { return fClocks.get();}
      
    private:
      
      std::unique_ptr<detinfo::DetectorClocksStandard> fClocks;
      
    };
  } //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::DetectorClocksServiceStandard,
                                   gar::detinfo::DetectorClocksService,
                                   LEGACY)

#endif // DETECTORCLOCKSSERVICESTANDARD_H
