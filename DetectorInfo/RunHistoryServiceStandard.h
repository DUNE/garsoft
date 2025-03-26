////////////////////////////////////////////////////////////////////////
// RunHistoryService_service.h
//
// Service interface for Run History functions
//
//  jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef RUNHISTORY_SERVICE_H
#define RUNHISTORY_SERVICE_H

#include "DetectorInfo/RunHistory.h"
#include "Utilities/IRunHistoryService.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace util class RunHistoryService : public IRunHistoryService {
public:
  RunHistoryService(fhicl::ParameterSet const& pset, ::art::ActivityRegistry& reg);

  virtual void reconfigure(fhicl::ParameterSet const& pset);
  void preBeginRun(const ::art::Run& run);

  virtual const provider_type* provider() const override { return fRH.get(); }

private:
  std::unique_ptr<detinfo::RunHistory> fRH;

}; // class RunHistoryService
} //namespace utils
DECLARE_ART_SERVICE_INTERFACE_IMPL(detinfo::RunHistoryService, util::IRunHistoryService, LEGACY)
#endif // GArProperties_SERVICE_H
