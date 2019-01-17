/**
* @file   Geometry_service.cc
* @brief  art framework interface to geometry description - implementation file
* @author brebel@fnal.gov
* @see    Geometry.h
*/

// class header
#include "Geometry/Geometry.h"

// gar includes
#include "SummaryDataProducts/RunData.h"
#include "Geometry/ExptGeoHelperInterface.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <vector>
#include <string>

namespace gar {
    namespace geo {


        //......................................................................
        // Constructor.
        Geometry::Geometry(fhicl::ParameterSet   const& pset,
        ::art::ActivityRegistry      & reg)
        : GeometryCore(pset)
        , fRelPath          (pset.get< std::string       >("RelativePath",     ""   ))
        , fForceUseFCLOnly  (pset.get< bool              >("ForceUseFCLOnly" , false))
        , fSortingParameters(pset.get<fhicl::ParameterSet>("SortingParameters", fhicl::ParameterSet() ))
        , fSegParameters(pset.get<fhicl::ParameterSet>("SegmentationAlgPars", fhicl::ParameterSet() ))
        {
            // add a final directory separator ("/") to fRelPath if not already there
            if (!fRelPath.empty() && (fRelPath.back() != '/')) fRelPath += '/';

            // register a callback to be executed when a new run starts
            reg.sPreBeginRun.watch(this, &Geometry::preBeginRun);

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

        } // Geometry::Geometry()


        //------------------------------------------------------------------------
        void Geometry::preBeginRun(::art::Run const& run)
        {
            // if we are requested to stick to the configured geometry, do nothing
            if (fForceUseFCLOnly) return;

            // check here to see if we need to load a new geometry.
            // get the detector id from the run object
            std::vector< ::art::Handle<sumdata::RunData> > rdcol;
            run.getManyByType(rdcol);
            if (rdcol.empty()) {
                mf::LogWarning("LoadNewGeometry") << "cannot find sumdata::RunData object to grab detector name\n"
                << "this is expected if generating MC files\n"
                << "using default geometry from configuration file\n";
                return;
            }

            // if the detector name is still the same, everything is fine
            std::string newDetectorName = rdcol.front()->DetName();
            if (DetectorName() == newDetectorName) return;

            // check to see if the detector name in the RunData
            // object has not been set.
            std::string const nodetname("nodetectorname");
            if (newDetectorName == nodetname) {
                LOG_WARNING("Geometry") << "Detector name not set: " << newDetectorName;
            } // if no detector name stored
            else {
                // the detector name is specified in the RunData object
                SetDetectorName(newDetectorName);
            }

            LoadNewGeometry(DetectorName() + ".gdml", DetectorName() + ".gdml", true);
        } // Geometry::preBeginRun()


        //......................................................................
        void Geometry::InitializeChannelMap()
        {
            // the channel map is responsible of calling the channel map configuration
            // of the geometry
            ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureChannelMapAlg(fSortingParameters, this);

            if ( ! ChannelMap() ) {
                throw cet::exception("ChannelMapLoadFail")
                << " failed to load new channel map";
            }

        } // Geometry::InitializeChannelMap()

        //......................................................................
        void Geometry::InitializeSegmentation()
        {
            // the channel map is responsible of calling the channel map configuration
            // of the geometry
            ::art::ServiceHandle<geo::ExptGeoHelperInterface>()->ConfigureECALSegmentationAlg(fSegParameters, this);

            if ( ! ECALSegmentationAlg() ) {
                throw cet::exception("ECALSegmentationAlgLoadFailed")
                << " failed to load the ECAL segmentation";
            }

        } // Geometry::InitializeSegmentation()

        //......................................................................
        void Geometry::LoadNewGeometry(std::string const& gdmlfile,
        std::string const& /* rootfile */,
        bool               bForceReload /* = false */)
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

            if(fRelPath.empty()){
                cet::search_path sp("FW_SEARCH_PATH");

                ;
                if( !sp.find_file(GDMLFileName, GDMLfile) ) {
                    throw cet::exception("Geometry")
                    << "cannot find the gdml geometry file:"
                    << "\n" << GDMLFileName
                    << "\nbail ungracefully.\n";
                }

                if( !sp.find_file(ROOTFileName, ROOTfile) ) {
                    throw cet::exception("Geometry")
                    << "cannot find the root geometry file:\n"
                    << "\n" << ROOTFileName
                    << "\nbail ungracefully.\n";
                }
            }
            else{
                GDMLfile = GDMLFileName;
                ROOTfile = GDMLFileName;
            }

            // initialize the geometry with the files we have found
            LoadGeometryFile(GDMLfile, ROOTfile, bForceReload);

            // now update the channel map
            InitializeChannelMap();

            // now init the ECAL Segmentation
            InitializeSegmentation();

        } // Geometry::LoadNewGeometry()

        DEFINE_ART_SERVICE(Geometry)
    } // namespace geo
} // namespace gar
