////////////////////////////////////////////////////////////////////////
//
//  GArProperties_plugin
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <iostream>

// GArSoft includes
#include "DetectorInfo/ECALPropertiesServiceStandard.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
gar::detinfo::ECALPropertiesServiceStandard::ECALPropertiesServiceStandard(fhicl::ParameterSet const& pset,
                                                                         ::art::ActivityRegistry &reg)
{
  fProp.reset(new detinfo::ECALPropertiesStandard());

  this->reconfigure(pset);
  reg.sPreBeginRun.watch(this, &ECALPropertiesServiceStandard::preBeginRun);
}

//----------------------------------------------
void gar::detinfo::ECALPropertiesServiceStandard::preBeginRun(const ::art::Run& run)
{
  fProp->Update(run.id().run());
}



//------------------------------------------------
/// \todo these values should eventually come from a database
void gar::detinfo::ECALPropertiesServiceStandard::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::ECALPropertiesServiceStandard,
                                  gar::detinfo::ECALPropertiesService)
