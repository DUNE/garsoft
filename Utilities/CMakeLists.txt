
add_subdirectory(scripts)

include_directories( $ENV{NUTOOLS_DIR}/include   )

art_make(
  LIB_LIBRARIES
    ${ROOT_EG}
    ${ROOT_CORE}
    nusimdata_SimulationBase
    Geometry
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    ${MF_MESSAGELOGGER}
    ${ROOT_EG}
    ${ROOT_BASIC_LIB_LIST}
    ${CETLIB_EXCEPT}
    ${FHICLCPP}
    ${Boost_SYSTEM_LIBRARY}

  MODULE_LIBRARIES
    Utilities
    nusimdata_SimulationBase
    Geometry
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    ${MF_MESSAGELOGGER}
    ${ROOT_EG}
    ${ROOT_BASIC_LIB_LIST}
    ${FHICLCPP}
    ${CETLIB_EXCEPT}
  )


install_headers()
install_fhicl()
install_source()
