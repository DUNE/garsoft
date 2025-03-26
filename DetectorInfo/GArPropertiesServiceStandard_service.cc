////////////////////////////////////////////////////////////////////////
//
//  GArProperties_plugin
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <iostream>

// GArSoft includes
#include "DetectorInfo/GArPropertiesServiceStandard.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
gar::detinfo::GArPropertiesServiceStandard::GArPropertiesServiceStandard(
  fhicl::ParameterSet const& pset,
  ::art::ActivityRegistry& reg)
{
  fProp.reset(new detinfo::GArPropertiesStandard());

  this->reconfigure(pset);
  reg.sPreBeginRun.watch(this, &GArPropertiesServiceStandard::preBeginRun);
}

//----------------------------------------------
void gar::detinfo::GArPropertiesServiceStandard::preBeginRun(const ::art::Run& run)
{
  fProp->Update(run.id().run());
}

//------------------------------------------------
/// \todo these values should eventually come from a database
void gar::detinfo::GArPropertiesServiceStandard::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::GArPropertiesServiceStandard,
                                  gar::detinfo::GArPropertiesService)
