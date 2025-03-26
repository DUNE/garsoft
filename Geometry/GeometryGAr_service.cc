/**
* @file   GeometryGAr_service.cc
* @brief  art framework interface to geometry description - implementation file
* @author brebel@fnal.gov
* @see    GeometryGAr.h
*/

// class header
#include "Geometry/GeometryGAr.h"

// gar includes
#include "Geometry/ExptGeoHelperInterface.h"
#include "SummaryDataProducts/RunData.h"

// Framework includes
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <algorithm> // std::min()
#include <string>
#include <vector>

namespace gar {
  namespace geo {

    //......................................................................
    // Constructor.
    GeometryGAr::GeometryGAr(fhicl::ParameterSet const& pset, ::art::ActivityRegistry& reg)
      : GeometryCore(pset)
      , fRelPath(pset.get<std::string>("RelativePath", ""))
      , fNonFatalConfCheck(pset.get<bool>("SkipConfigurationCheck", false))
      , fSortingParameters(
          pset.get<fhicl::ParameterSet>("SortingParameters", fhicl::ParameterSet()))
      , fSegParameters(pset.get<fhicl::ParameterSet>("SegmentationAlgPars", fhicl::ParameterSet()))
    {
      // add a final directory separator ("/") to fRelPath if not already there
      if (!fRelPath.empty() && (fRelPath.back() != '/')) fRelPath += '/';

      // register a callback to be executed when a new run starts
      reg.sPreBeginRun.watch(this, &GeometryGAr::preBeginRun);

      //......................................................................
      // 5.15.12 BJR: use the gdml file for both the fGDMLFile and fROOTFile
      // variables as ROOT v5.30.06 is once again able to read in gdml files
      // during batch operation, in this case think of fROOTFile meaning the
      // file used to make the ROOT TGeoManager.  I don't want to remove
      // the separate variables in case ROOT breaks again
      std::string GDMLFileName = pset.get<std::string>("GDML");
      std::string ROOTFileName = pset.get<std::string>("GDML");

      // load the geometry
      LoadNewGeometry(GDMLFileName, ROOTFileName);

      FillGeometryConfigurationInfo(pset);

    } // GeometryGAr::Geometry()

    //------------------------------------------------------------------------
    void GeometryGAr::preBeginRun(::art::Run const& run)
    {
      sumdata::GeometryConfigurationInfo const inputGeomInfo = ReadConfigurationInfo(run);
      if (!CheckConfigurationInfo(inputGeomInfo)) {
        if (fNonFatalConfCheck) {
          // disable the non-fatal option if you need the details
          MF_LOG_WARNING("Geometry") << "Geometry used for " << run.id()
                                     << " is incompatible with the one configured in the job.";
        }
        else {
          throw cet::exception("Geometry")
            << "Geometry used for run " << run.id()
            << " is incompatible with the one configured in the job!"
            << "\n=== job configuration " << std::string(50, '=') << "\n"
            << fConfInfo << "\n=== run configuration " << std::string(50, '=') << "\n"
            << inputGeomInfo << "\n======================" << std::string(50, '=') << "\n";
        }
      }
    } // GeometryGAr::preBeginRun()

    //......................................................................
    void GeometryGAr::InitializeSegmentations()
    {
      // the channel map is responsible of calling the channel map configuration
      // of the TPC geometry
      ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureChannelMapAlg(
        fSortingParameters, this);

      if (!ChannelMap()) {
        throw cet::exception("ChannelMapLoadFail") << " failed to load new channel map";
      }

      // the channel map is responsible of calling the channel map configuration
      // of the ECAL geometry
      fECALSegParameters =
        fSegParameters.get<fhicl::ParameterSet>("ECALSegmentationAlgPars", fhicl::ParameterSet());
      if (not fECALSegParameters.is_empty()) {
        ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureECALSegmentationAlg(
          fECALSegParameters, this);

        if (!ECALSegmentationAlg()) {
          throw cet::exception("ECALSegmentationAlgLoadFailed")
            << " failed to load the ECAL segmentation";
        }
      }

      // the channel map is responsible of calling the channel map configuration
      // of the Tracker Sc geometry
      fMinervaSegParameters = fSegParameters.get<fhicl::ParameterSet>("MinervaSegmentationAlgPars",
                                                                      fhicl::ParameterSet());
      if (not fMinervaSegParameters.is_empty()) {
        ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureMinervaSegmentationAlg(
          fMinervaSegParameters, this);

        if (!MinervaSegmentationAlg()) {
          throw cet::exception("MinervaSegmentationAlgLoadFailed")
            << " failed to load the Minerva segmentation";
        }
      }

      // the channel map is responsible of calling the channel map configuration
      // of the MuID geometry
      fMuIDSegParameters =
        fSegParameters.get<fhicl::ParameterSet>("MuIDSegmentationAlgPars", fhicl::ParameterSet());
      if (not fMuIDSegParameters.is_empty()) {
        ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureMuIDSegmentationAlg(
          fMuIDSegParameters, this);

        if (!MuIDSegmentationAlg()) {
          throw cet::exception("MuIDSegmentationAlgLoadFailed")
            << " failed to load the MuID segmentation";
        }
      }

    } // GeometryGAr::InitializeSegmentation()

    //......................................................................
    void GeometryGAr::LoadNewGeometry(std::string const& gdmlfile,
                                      std::string const& /* rootfile */,
                                      bool bForceReload /* = false */)
    {
      // start with the relative path
      std::string GDMLFileName(fRelPath), ROOTFileName(fRelPath);

      // add the base file names
      ROOTFileName.append(gdmlfile);
      GDMLFileName.append(gdmlfile);

      // Search all reasonable locations for the GDML file that contains
      // the detector geometry.
      // cet::search_path constructor decides if initialized value is a path
      // or an environment variable
      std::string GDMLfile("");
      std::string ROOTfile("");

      if (fRelPath.empty()) {
        cet::search_path sp("FW_SEARCH_PATH");

        ;
        if (!sp.find_file(GDMLFileName, GDMLfile)) {
          throw cet::exception("Geometry") << "cannot find the gdml geometry file:"
                                           << "\n"
                                           << GDMLFileName << "\nbail ungracefully.\n";
        }

        if (!sp.find_file(ROOTFileName, ROOTfile)) {
          throw cet::exception("Geometry") << "cannot find the root geometry file:\n"
                                           << "\n"
                                           << ROOTFileName << "\nbail ungracefully.\n";
        }
      }
      else {
        GDMLfile = GDMLFileName;
        ROOTfile = GDMLFileName;
      }

      // initialize the geometry with the files we have found
      LoadGeometryFile(GDMLfile, ROOTfile, bForceReload);

      //Get detector parameters
      GetGeometryParameters();

      // now init the detector Segmentations (Channel Map, ECAL etc...)
      InitializeSegmentations();

      // Finish the geometry when it uses the segmentation
      FinalizeGeometryParameters();

      //Print the geometry parameters
      PrintGeometry();
    } // GeometryGAr::LoadNewGeometry()

    //......................................................................
    void GeometryGAr::FillGeometryConfigurationInfo(fhicl::ParameterSet const& config)
    {

      gar::sumdata::GeometryConfigurationInfo confInfo;
      confInfo.dataVersion = gar::sumdata::GeometryConfigurationInfo::DataVersion_t{2};

      // version 1+:
      confInfo.detectorName = DetectorName();

      // version 2+:
      confInfo.geometryServiceConfiguration = config.to_indented_string();
      fConfInfo = std::move(confInfo);

      MF_LOG_DEBUG("Geometry") << "Geometry configuration information:\n" << fConfInfo;

    } // Geometry::FillGeometryConfigurationInfo()

    //......................................................................
    bool GeometryGAr::CheckConfigurationInfo(
      gar::sumdata::GeometryConfigurationInfo const& other) const
    {
      MF_LOG_DEBUG("Geometry") << "New geometry information:\n" << other;
      return CompareConfigurationInfo(fConfInfo, other);
    } // Geometry::CheckConfigurationInfo()

    //......................................................................
    gar::sumdata::GeometryConfigurationInfo const& GeometryGAr::ReadConfigurationInfo(
      art::Run const& run)
    {

      try {
        return run.getProduct<gar::sumdata::GeometryConfigurationInfo>(
          art::InputTag{"GeometryGArConfigurationWriter"});
      }
      catch (art::Exception const& e) {
        throw art::Exception{e.categoryCode(),
                             "Can't read geometry configuration information.\n"
                             "Is `GeometryGArConfigurationWriter` service configured?\n",
                             e};
      }

    } // Geometry::ReadConfigurationInfo()

    //......................................................................
    bool GeometryGAr::CompareConfigurationInfo(gar::sumdata::GeometryConfigurationInfo const& A,
                                               gar::sumdata::GeometryConfigurationInfo const& B)
    {
      /*
          * Implemented criteria:
          *
          * * both informations must be valid
          * * the detector names must exactly match
          *
          */

      if (!A.isDataValid()) {
        MF_LOG_WARNING("Geometry") << "Geometry::CompareConfigurationInfo(): "
                                      "invalid version for configuration A:\n"
                                   << A;
        return false;
      }
      if (!B.isDataValid()) {
        MF_LOG_WARNING("Geometry") << "Geometry::CompareConfigurationInfo(): "
                                      "invalid version for configuration B:\n"
                                   << B;
        return false;
      }

      return true;
    } // CompareConfigurationInfo()

    DEFINE_ART_SERVICE(GeometryGAr)
  } // namespace geo
} // namespace gar
