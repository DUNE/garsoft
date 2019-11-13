include_directories($ENV{GEANT4_FQ_DIR}/include)

art_make( NO_PLUGINS
          LIB_LIBRARIES Geometry
	                   		CoreUtils
                        ${MF_MESSAGELOGGER}

                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                  			${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
                  			${Boost_SYSTEM_LIBRARY}
        )


simple_plugin(DetectorClocksServiceStandard "service"
              DetectorInfo
              RawDataProducts
  	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Persistency_Common canvas
	      art_Persistency_Provenance canvas
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      art_Utilities canvas
              ${ART_FRAMEWORK_IO_ROOT_ROOTDB}
              ${MF_MESSAGELOGGER}
              ${FHICLCPP}
              ${CETLIB_EXCEPT}
	      ${SQLITE3}
              ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(DetectorPropertiesServiceStandard "service"
              DetectorInfo
              Geometry_Geometry_service
              Geometry
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      art_Persistency_Common canvas
	      art_Persistency_Provenance canvas
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      art_Utilities canvas
              ${ART_FRAMEWORK_IO_ROOT_ROOTDB}
              ${MF_MESSAGELOGGER}

              ${CETLIB_EXCEPT}
	      ${SQLITE3}
              ${ROOT_BASIC_LIB_LIST}
             )

simple_plugin(GArPropertiesServiceStandard "service"
              DetectorInfo
	      DetectorInfo_DetectorPropertiesServiceStandard_service
              art_Utilities canvas
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${CETLIB_EXCEPT}
             )

simple_plugin(ECALPropertiesServiceStandard "service"
             DetectorInfo
             DetectorInfo_DetectorPropertiesServiceStandard_service
             art_Utilities canvas
             ${ART_FRAMEWORK_SERVICES_REGISTRY}
             ${MF_MESSAGELOGGER}
             ${CETLIB_EXCEPT}
             )

simple_plugin(MPDMagneticField "service"
             art_Utilities canvas
             ${ART_FRAMEWORK_SERVICES_REGISTRY}
             ${MF_MESSAGELOGGER}
             ${CETLIB_EXCEPT}
	     ${ROOT_BASIC_LIB_LIST}
	     ${ROOT_GEOM}
             ${ROOT_CORE}
             ${G4GEOMETRY}
             ${G4GLOBAL}
             )

install_headers()
install_fhicl()
install_source()

#add_subdirectory(test)
