include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )
include_directories( $ENV{GSL_FQ_DIR}/include     )

set(GSLLIBS
$ENV{GSL_FQ_DIR}/lib/libgsl.so
$ENV{GSL_FQ_DIR}/lib/libgslcblas.so
)

MESSAGE(STATUS "GSLLIBS is ${GSLLIBS}")

cet_make( LIBRARIES
                       Geometry
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
	                     canvas
	                     nutools_MagneticField_MagneticField_service
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                       ${ROOT_BASIC_LIB_LIST}
                       ${CLHEP}
                       ${Boost_SYSTEM_LIBRARY}
                       ${GSLLIBS}
                       NO_DICTIONARY )

# Simulation base is included in making the dictionary for associations
# between the MC truth and reconstructed quantities when cheating

art_dictionary( DICTIONARY_LIBRARIES ReconstructionDataProducts
                                     RawDataProducts
                                     nusimdata_SimulationBase)

install_headers()
install_source()
