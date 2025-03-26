////////////////////////////////////////////////////////////////////////
// DetectorClocksService.h
//
// Pure virtual service interface for DetectorClocks functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETECTORCLOCKSSERVICE_H
#define DETECTORCLOCKSSERVICE_H

#include "CoreUtils/ServiceUtil.h"
#include "DetectorInfo/DetectorClocks.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {
    class DetectorClocksServiceGAr {

    public:
      typedef detinfo::DetectorClocks provider_type;

    public:
      virtual ~DetectorClocksServiceGAr() = default;

      virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;
      virtual const detinfo::DetectorClocks* provider() const = 0;

    }; // class DetectorClocksServiceGAr
  }    //namespace detinfo
} //gar

DECLARE_ART_SERVICE_INTERFACE(gar::detinfo::DetectorClocksServiceGAr, LEGACY)

#endif // DETECTORCLOCKSSERVICE_H
