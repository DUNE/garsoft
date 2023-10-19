# ======================================================================
#  garsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/garsoft/ups/setup_for_development
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


# use cmake 2.8 or later
cmake_minimum_required (VERSION 3.19 FATAL_ERROR)

project(garsoft VERSION 02.20.00)

find_package( cetbuildtools REQUIRED)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
                              $ENV{CETBUILDTOOLS_DIR}/Modules)

include(CetCMakeEnv)
cet_cmake_env()

# -D_GLIBCXX_USE_NANOSLEEP is needed to make std::this_thread::sleep_for available.
cet_set_compiler_flags(DIAGS CAUTIOUS 
                       WERROR
                       NO_UNDEFINED
                       EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
                       EXTRA_CXX_FLAGS -D_GLIBCXX_USE_NANOSLEEP ${ARTDFLAG} -Wno-variadic-macros)
cet_report_compiler_flags()

find_package(art REQUIRED EXPORT)
include_directories( $ENV{ART_INC} )

find_package(art_root_io REQUIRED EXPORT)
find_package(canvas_root_io REQUIRED)
find_package(cetlib REQUIRED EXPORT)
find_package(cetlib_except REQUIRED EXPORT)
find_package(fhiclcpp REQUIRED EXPORT)
find_package(messagefacility REQUIRED EXPORT)

find_package(nug4  REQUIRED EXPORT)
include_directories( $ENV{NUG4_INC} )
link_directories( $ENV{NUG4_FQ_DIR}/lib )

find_package(nugen REQUIRED EXPORT)
include_directories( $ENV{NUGEN_INC} )
find_package(nurandom REQUIRED EXPORT)
find_package(nusimdata REQUIRED EXPORT)
find_package(nutools REQUIRED EXPORT)
find_package(nuevdb REQUIRED EXPORT)

find_package(Boost REQUIRED COMPONENTS headers math_tr1 REQUIRED EXPORT) # boost::math::policies, boost::math::ellint*
find_package(CLHEP COMPONENTS Evaluator Geometry Random Vector REQUIRED EXPORT)
find_package(CRY REQUIRED EXPORT) # find module in nutools
find_package(GENIE REQUIRED EXPORT) # find module in nugen
include_directories( $ENV{GENIE_INC}/GENIE )
include_directories( $ENV{GENIE_INC} )

find_package(Geant4 REQUIRED EXPORT)
include_directories( $ENV{GEANT4_FQ_DIR}/include )
find_package(LibXml2 REQUIRED EXPORT)
find_package(PostgreSQL REQUIRED EXPORT)
find_package(ROOT COMPONENTS Core EG GenVector Geom Gpad Hist MathCore Matrix Physics RIO RooFitCore Tree REQUIRED EXPORT)
find_package(Range-v3 REQUIRED EXPORT)
find_package(SQLite3 REQUIRED EXPORT)
find_package(ifdhc REQUIRED EXPORT)
find_package(log4cpp REQUIRED EXPORT)
find_package(nlohmann_json REQUIRED EXPORT)
find_package(HDF5 REQUIRED EXPORT)

# Issue with capitalization here.  UPS product is called highfive and the package is called HighFive
# the product itself has a mixture
find_package( HighFive )
include_directories( $ENV{HIGHFIVE_INC} )

find_package( dunedaqdataformats )
include_directories( $ENV{DUNEDAQDATAFORMATS_INC} )
find_package( dunedetdataformats )
include_directories( $ENV{DUNEDETDATAFORMATS_INC} )

# find_package( GSL REQUIRED EXPORT )
find_ups_product( GSL )
find_package( EDepSim REQUIRED EXPORT )
find_package( garana REQUIRED EXPORT )
find_ups_product( pythia )

#message(STATUS "GARANA_DIR is $ENV{GARANA_DIR}")
#message(STATUS "GARANA_LIB is $ENV{GARANA_LIB}")

#message(STATUS "ROOTSYS is ${ROOTSYS}")
#message(STATUS "GSL_FQ_DIR is $ENV{GSL_FQ_DIR}")

# macros
include(ArtMake)
include(ArtDictionary)
include(BuildPlugins)
include(CetRootCint)

# source
add_subdirectory (Utilities)
add_subdirectory (EventDisplay)
add_subdirectory (EventGenerator)
add_subdirectory (CoreUtils)
add_subdirectory (DetectorInfo)
add_subdirectory (GArG4)
add_subdirectory (Geometry)
add_subdirectory (MCCheater)
add_subdirectory (RawDataProducts)
add_subdirectory (ReadoutSimulation)
add_subdirectory (ReconstructionDataProducts)
add_subdirectory (SimulationDataProducts)
add_subdirectory (SummaryDataProducts)
add_subdirectory (TestUtils)
add_subdirectory (Reco)
add_subdirectory (RecoAlg)
add_subdirectory (Ana)
add_subdirectory (Filters)
add_subdirectory (StandardRecord)
add_subdirectory (Prototypes)

message(STATUS "GSL_FQ_DIR is $ENV{GSL_FQ_DIR}")

# Pandora
option(BUILD_W_PANDORA "Build with Pandora" ON )

# to do -- get rid of the cet_find_library calls here to satisfy the migration script

if( BUILD_W_PANDORA )
find_ups_product( pandora )
cet_find_library( PANDORASDK NAMES PandoraSDK PATHS ENV PANDORA_LIB )
cet_find_library( PANDORAMONITORING NAMES PandoraMonitoring PATHS ENV PANDORA_LIB )
cet_find_library( LCPANDORACONTENT NAMES LCContent PATHS ENV PANDORA_LIB )
add_subdirectory (Pandora)
endif()

# tests
add_subdirectory (test)

# ups - table and config files
add_subdirectory(ups)
# packaging utility
include(UseCPack)
