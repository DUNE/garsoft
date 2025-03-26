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

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Provenance/ScheduleContext.h"
#include "fhiclcpp/ParameterSet.h"

#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "DetectorInfo/DetectorClocksStandardGAr.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {
    class DetectorClocksServiceStandardGAr : public DetectorClocksServiceGAr {
    public:
      DetectorClocksServiceStandardGAr(fhicl::ParameterSet const& pset,
                                       ::art::ActivityRegistry& reg);

      virtual void reconfigure(fhicl::ParameterSet const& pset) override;
      void preBeginRun(::art::Run const& run);
      void preProcessEvent(::art::Event const& evt, art::ScheduleContext);
      void postOpenFile(std::string const& filename);

      virtual const provider_type* provider() const override { return fClocks.get(); }

    private:
      std::unique_ptr<detinfo::DetectorClocksStandardGAr> fClocks;
    };
  } //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::DetectorClocksServiceStandardGAr,
                                   gar::detinfo::DetectorClocksServiceGAr,
                                   LEGACY)

#endif // DETECTORCLOCKSSERVICESTANDARD_H
