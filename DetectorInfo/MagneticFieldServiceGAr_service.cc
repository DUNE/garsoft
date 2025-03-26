////////////////////////////////////////////////////////////////////////
// \file MagneticFieldServiceStandard_service.cc
//
// \brief implementation of class for storing/accessing magnetic fields
//
// \author ebrianne@fnal.gov
//
////////////////////////////////////////////////////////////////////////

// C++ language includes

//nug4 includes
#include "DetectorInfo/GArMagneticField.h"
#include "DetectorInfo/MagneticFieldServiceGAr.h"

//Framework includes
#include "art/Framework/Principal/Run.h" // for Run
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/GlobalSignal.h"
#include "canvas/Persistency/Provenance/RunID.h"

//-----------------------------------------------
mag::MagneticFieldServiceGAr::MagneticFieldServiceGAr(fhicl::ParameterSet const& pset,
                                                      art::ActivityRegistry& reg)
  : fProp{pset}
{
  reg.sPreBeginRun.watch(this, &MagneticFieldServiceGAr::preBeginRun);
}

//----------------------------------------------
void mag::MagneticFieldServiceGAr::preBeginRun(const art::Run&) {}

//------------------------------------------------
void mag::MagneticFieldServiceGAr::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp.reconfigure(pset);
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(mag::MagneticFieldServiceGAr, mag::MagneticFieldService)
