////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef GArPropertiesSERVICE_H
#define GArPropertiesSERVICE_H

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/GArProperties.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {
    class GArPropertiesService {

    public:
      typedef detinfo::GArProperties provider_type;

    public:
      virtual ~GArPropertiesService() = default;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const detinfo::GArProperties* provider() const = 0;

    }; // class GArPropertiesService
  }    //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE(gar::detinfo::GArPropertiesService, LEGACY)
#endif // GArPropertiesSERVICE_H
