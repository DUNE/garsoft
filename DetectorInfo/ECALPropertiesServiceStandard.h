////////////////////////////////////////////////////////////////////////
// ECALPropertiesServiceStandard.h
//
// Service interface for Utility ECAL functions
//
//  eldwan.brianne@desy.de
//
////////////////////////////////////////////////////////////////////////
#ifndef ECALPropertiesSERVICESTANDARD_H
#define ECALPropertiesSERVICESTANDARD_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "DetectorInfo/ECALPropertiesStandard.h"
#include "DetectorInfo/ECALPropertiesService.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo{
    class ECALPropertiesServiceStandard : public ECALPropertiesService {
    public:

      // this enables art to print the configuration help:
      using Parameters = ::art::ServiceTable<gar::detinfo::ECALPropertiesStandard::ConfigurationParameters_t>;

      ECALPropertiesServiceStandard(fhicl::ParameterSet   const& pset,
                                   ::art::ActivityRegistry      & reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset) override;
      void   preBeginRun(const ::art::Run& run);

      virtual const  provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<detinfo::ECALPropertiesStandard> fProp;

    }; // class GArPropertiesServiceStandard
  } //namespace detinfo
} // gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::ECALPropertiesServiceStandard,
                                   gar::detinfo::ECALPropertiesService,
                                   LEGACY)

#endif // ECALPropertiesSERVICESTANDARD_H
