include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

art_make(BASENAME_ONLY
          LIBRARY_NAME GArAna
	     LIB_LIBRARIES nusimdata_SimulationBase
                       ReconstructionDataProducts
                       ${ART_FRAMEWORK_CORE}
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
                       canvas
                       Geometry
                       ${MF_MESSAGELOGGER}
                       ${MF_UTILITIES}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                       ${ROOT_BASIC_LIB_LIST}
                       ${CLHEP}
                       ${Boost_SYSTEM_LIBRARY}
      MODULE_LIBRARIES Geometry
                       GArAna
                       Utilities
                       ReadoutSimulation
                       SimulationDataProducts
                       RawDataProducts
                       ReconstructionDataProducts
                       MCCheater_BackTracker_service
                       nusimdata_SimulationBase
                       nutools_RandomUtils_NuRandomService_service
                       nutools_MagneticField_MagneticField_service
                       ${ART_FRAMEWORK_CORE}
                       ${ART_FRAMEWORK_IO_SOURCES}
                       ${ART_FRAMEWORK_PRINCIPAL}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
                       canvas
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                       ${ROOT_BASIC_LIB_LIST}
         )

install_headers()
install_fhicl()
install_source()
