#ifndef MAGNETICFIELDSERVICEGAR_H
#define MAGNETICFIELDSERVICEGAR_H

#include "DetectorInfo/GArMagneticField.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"

namespace mag {
  class MagneticFieldServiceGAr : public MagneticFieldService {
  public:
    MagneticFieldServiceGAr(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

  private:
    void reconfigure(fhicl::ParameterSet const& pset);
    void preBeginRun(const art::Run& run);

    const provider_type* provider() const override { return &fProp; }

    mag::GArMagneticField fProp;
  };
}

DECLARE_ART_SERVICE_INTERFACE_IMPL(mag::MagneticFieldServiceGAr, mag::MagneticFieldService, SHARED)

#endif
