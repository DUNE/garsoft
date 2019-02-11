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
cmake_minimum_required (VERSION 2.8)

project(garsoft)

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_DIR}/Modules
                      $ENV{ART_DIR}/Modules
		                  $ENV{CETBUILDTOOLS_DIR}/Modules
		                  ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

# -D_GLIBCXX_USE_NANOSLEEP is needed to make std::this_thread::sleep_for available.
cet_set_compiler_flags(DIAGS CAUTIOUS WERROR
                       NO_UNDEFINED
                       EXTRA_FLAGS -pedantic
                       EXTRA_CXX_FLAGS -D_GLIBCXX_USE_NANOSLEEP ${ARTDFLAG}
                                       -Wno-unused-local-typedefs -Wno-variadic-macros)
cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( art            )
find_ups_product( canvas         )
find_ups_product( messagefacility)
find_ups_product( genie          )
#find_ups_product( dk2nu          )
find_ups_product( nutools        )
find_ups_product( nusimdata      )
find_ups_product( libwda         )
find_ups_product( postgresql     )
find_ups_product( cetbuildtools  )
find_ups_boost  ()
find_ups_root   ()

message(STATUS "ROOTSYS is ${ROOTSYS}")

# macros
include(ArtMake)
include(ArtDictionary)
include(BuildPlugins)
include(CetRootCint)

# source
add_subdirectory(EventDisplay)
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
add_subdirectory (Utilities)
add_subdirectory (Reco)
add_subdirectory (RecoAlg)
add_subdirectory (Ana)

# tests
add_subdirectory (test)

# ups - table and config files
add_subdirectory(ups)
# packaging utility
include(UseCPack)
