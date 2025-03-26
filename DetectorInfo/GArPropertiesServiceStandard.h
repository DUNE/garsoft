////////////////////////////////////////////////////////////////////////
// GArPropertiesServiceStandard.h
//
// Service interface for Utility LAr functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef GArPropertiesSERVICESTANDARD_H
#define GArPropertiesSERVICESTANDARD_H

#include "DetectorInfo/GArPropertiesService.h"
#include "DetectorInfo/GArPropertiesStandard.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {
    class GArPropertiesServiceStandard : public GArPropertiesService {
    public:
      // this enables art to print the configuration help:
      using Parameters =
        ::art::ServiceTable<gar::detinfo::GArPropertiesStandard::ConfigurationParameters_t>;

      GArPropertiesServiceStandard(fhicl::ParameterSet const& pset, ::art::ActivityRegistry& reg);

      virtual void reconfigure(fhicl::ParameterSet const& pset) override;
      void preBeginRun(const ::art::Run& run);

      virtual const provider_type* provider() const override { return fProp.get(); }

    private:
      std::unique_ptr<detinfo::GArPropertiesStandard> fProp;

    }; // class GArPropertiesServiceStandard
  }    //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::GArPropertiesServiceStandard,
                                   gar::detinfo::GArPropertiesService,
                                   LEGACY)

#endif // GArPropertiesSERVICESTANDARD_H
