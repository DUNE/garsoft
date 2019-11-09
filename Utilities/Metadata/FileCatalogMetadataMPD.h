////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataMPD_service.cc.  
//
// Purpose:  Art service adds dune-specific per-job sam metadata.
//
//           This service does not have user-callable methods.  Simply
//           add to an art configuration in services.user block of job
//           file.
//
// Created:  7-Nov-2019
//  copied from DUNE for MPD by T. Junk
//  split from the original by T. Yang so we can include this
// Copied FileCatalogMetadataMicroBooNE_service.cc by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#ifndef FILECATALOGMETADATAMPD_H
#define FILECATALOGMETADATAMPD_H

#include <string>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

namespace util {

  // Class declaration.

  class FileCatalogMetadataMPD
  {
  public:

    FileCatalogMetadataMPD(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~FileCatalogMetadataMPD() = default;

    const std::string&  Generators()		const   { return fGenerators; }		  
    const std::string&  TriggerListVersion()	const   { return fTriggerListVersion; }	  
    const std::string&  BeamFluxID()		const   { return fBeamFluxID; }		  
    const std::string&  MCName()	        const   { return fMCName; }		  
    const std::string&  NeutrinoFlavors()	const   { return fNeutrinoFlavors; }	  
    const std::string&  Miscellaneous()         const   { return fMiscellaneous; }	  
    const std::string&  GeometryVersion()	const   { return fGeometryVersion; }	  
    const std::string&  Overlay()		const   { return fOverlay; }		  
    const std::string&  Run_Type()              const   { return fRunType; }	  
    const std::string&  Data_Tier()             const   { return fDataTier; }	  
    const std::string&  ApplName()              const   { return fApplName; }	  
    const std::string&  EnvVarVer()             const   { return fEnvVarVer; }	  
    const std::string&  Systematic()            const   { return fSystematic; }	  
    const std::string&  HornPolarity()          const   { return fHornPolarity; }	  


  private:

    // Callbacks.

    void postBeginJob();

    // Data members.

    std::string fGenerators;
    std::string fTriggerListVersion;
    std::string fBeamFluxID;
    std::string fMCName;
    std::string fNeutrinoFlavors;
    std::string fMiscellaneous;
    std::string fGeometryVersion;
    std::string fOverlay;
    std::string fRunType;
    std::string fDataTier;
    std::string fApplName;
    std::string fEnvVarVer;
    std::string fSystematic;
    std::string fHornPolarity;
  };


} // namespace util

DECLARE_ART_SERVICE(util::FileCatalogMetadataMPD, LEGACY)

#endif
