////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataMPD_service.cc.  
//
// Purpose:  Art service adds MPD-specific per-job sam metadata.
//
//           This service does not have user-callable methods.  Simply
//           add to an art configuration in services.user block of job
//           file.
//
// MPD version based on DUNE version, 7 Nov 2019, T. Junk
// DUNE version Created:  3-Dec-2014,  T. Yang
// Copied FileCatalogMetadataMicroBooNE_service.cc by H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include "FileCatalogMetadataMPD.h"
#include <stdlib.h>

  //--------------------------------------------------------------------
  // Constructor.

util::FileCatalogMetadataMPD::
  FileCatalogMetadataMPD(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    // Get parameters.

    fGenerators =         pset.get<std::string>("Generators","");
    fTriggerListVersion = pset.get<std::string>("TriggerListVersion","");
    fBeamFluxID =         pset.get<std::string>("BeamFluxID","");
    fMCName =             pset.get<std::string>("MCName","");
    fNeutrinoFlavors =    pset.get<std::string>("NeutrinoFlavors","");
    fMiscellaneous =      pset.get<std::string>("Miscellaneous","");
    fGeometryVersion =    pset.get<std::string>("GeometryVersion","");
    fOverlay =            pset.get<std::string>("Overlay","");
    fRunType =            pset.get<std::string>("run_type","");
    fDataTier =           pset.get<std::string>("data_tier","");
    fApplName =           pset.get<std::string>("ApplName","garsoft");
    fEnvVarVer =          pset.get<std::string>("EnvVarVer","GARSOFT_VERSION");
    fHornPolarity =       pset.get<std::string>("HornPolarity","");

    // Register for callbacks.

    reg.sPostBeginJob.watch(this, &FileCatalogMetadataMPD::postBeginJob);
  }

  //--------------------------------------------------------------------
  // PostBeginJob callback.
  // Insert per-job metadata via FileCatalogMetadata service.
  void util::FileCatalogMetadataMPD::postBeginJob()
  {
    // Get art metadata service.

    art::ServiceHandle<art::FileCatalogMetadata> mds;

    // Add metadata.

    if (fGenerators!="")         mds->addMetadata("DUNE.Generators",              fGenerators);
    if (fGenerators!="")         mds->addMetadata("DUNE_MC.Generators",           fGenerators);
    if (fTriggerListVersion!="") mds->addMetadata("DUNE_MC.trigger-list-version", fTriggerListVersion);
    if (fBeamFluxID!="")         mds->addMetadata("DUNE_MC.beam_flux_ID",         fBeamFluxID);
    if (fMCName!="")             mds->addMetadata("DUNE_MC.name",                 fMCName);
    if (fNeutrinoFlavors!="")    mds->addMetadata("DUNE_MC_neutrino_flavors",     fNeutrinoFlavors);
    if (fMiscellaneous!="")      mds->addMetadata("DUNE_MC.miscellaneous",        fMiscellaneous);
    if (fGeometryVersion!="")    mds->addMetadata("DUNE_MC.GeometryVersion",      fGeometryVersion);
    if (fOverlay!="")            mds->addMetadata("DUNE_MC.Overlay",              fOverlay);
    if (fSystematic!="")         mds->addMetadata("DUNE_MC.Systematic",           fSystematic);
    if (fRunType!="")            mds->addMetadata("run_type",                     fRunType);
    if (fDataTier!="")           mds->addMetadata("data_tier",                    fDataTier); // e.g. "reco"
    if (fApplName!="")           mds->addMetadata("application.family",           fApplName); // e.g. "garsoft"
    if (fDataTier!="")           mds->addMetadata("application.name",             fDataTier); // e.g. "reco"
    if (fEnvVarVer!="")          mds->addMetadata("application.version",          getenv(fEnvVarVer.c_str())); 
    if (fHornPolarity!="")       mds->addMetadata("DUNE_MC.HornPolarity",         fHornPolarity);
  }

DEFINE_ART_SERVICE(util::FileCatalogMetadataMPD)
