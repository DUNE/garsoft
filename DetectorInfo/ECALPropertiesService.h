////////////////////////////////////////////////////////////////////////
// DetectorPropertiesService.h
//
// Pure virtual service interface for DetectorProperties functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef ECALPropertiesSERVICE_H
#define ECALPropertiesSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "DetectorInfo/ECALProperties.h"
#include "CoreUtils/ServiceUtil.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    class ECALPropertiesService {

    public:
      typedef detinfo::ECALProperties provider_type;

    public:
      virtual ~ECALPropertiesService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const  detinfo::ECALProperties* provider() const = 0;

    }; // class ECALPropertiesService
  } //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE(gar::detinfo::ECALPropertiesService, LEGACY)
#endif // ECALPropertiesSERVICE_H
