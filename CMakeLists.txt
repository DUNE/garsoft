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

project(garsoft VERSION 02.16.00)

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

# these are minimum required versions, not the actual product versions
find_ups_product( art            )
find_ups_product( art_root_io    )
find_ups_product( canvas         )
find_ups_product( messagefacility)
find_ups_product( genie          )
#find_ups_product( dk2nu         )
find_ups_product( nutools        )
find_ups_product( nuevdb         )
find_ups_product( nusimdata      )
find_ups_product( nusystematics  )
find_ups_product( nug4           )
find_ups_product( nugen          )
find_ups_product( nurandom       )
find_ups_product( libwda         )
find_ups_product( postgresql     )
find_ups_geant4(  geant4         )
find_ups_boost  ()
find_ups_root   ()
find_ups_product( gsl )
find_ups_product( edepsim )
find_ups_product( garana )
find_ups_product( cetlib )

#message(STATUS "GARANA_DIR is $ENV{GARANA_DIR}")
#message(STATUS "GARANA_LIB is $ENV{GARANA_LIB}")

cet_find_library( EDEPSIMIO NAMES edepsim_io PATHS ENV EDEPSIM_LIB )

cet_find_library( GARANAPRODUCTS NAMES GaranaProducts PATHS ENV GARANA_LIB)

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

message(STATUS "GSL_FQ_DIR is $ENV{GSL_FQ_DIR}")

# Pandora
option(BUILD_W_PANDORA "Build with Pandora" ON )

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
