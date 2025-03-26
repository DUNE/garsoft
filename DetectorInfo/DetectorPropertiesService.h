////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORPROPERTIESSERVICE_H
#define DETECTORPROPERTIESSERVICE_H

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorProperties.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {
    class DetectorPropertiesService {

    public:
      typedef detinfo::DetectorProperties provider_type;

    public:
      virtual ~DetectorPropertiesService() = default;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const detinfo::DetectorProperties* provider() const = 0;

    }; // class DetectorPropertiesService
  }    //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE(gar::detinfo::DetectorPropertiesService, LEGACY)
#endif // DETECTORPROPERTIESSERVICE_H
