
art_make( NO_PLUGINS
          LIB_LIBRARIES Geometry
                        CoreUtils
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${CETLIB_EXCEPT}
                        ${FHICLCPP}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
                        ${Boost_SYSTEM_LIBRARY}
                        art_Utilities canvas
			nug4_MagneticFieldServices_MagneticFieldServiceStandard_service
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${MF_MESSAGELOGGER}
                        ${CETLIB_EXCEPT}
                        ${ROOT_GEOM}
                        ${ROOT_CORE}
                        ${G4GEOMETRY}
                        ${G4GLOBAL}
                        ${ART_FRAMEWORK_SERVICES_BASIC}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                        ${ART_ROOT_IO_TFILE_SUPPORT}
                        ${CLHEP}
        )


simple_plugin(DetectorClocksServiceStandardGAr "service"
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
              ${ART_ROOT_IO_ROOTDB}
              ${SQLITE3}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
             )

simple_plugin(DetectorPropertiesServiceStandardGAr "service"
              DetectorInfo
              Geometry_GeometryGAr_service
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
              ${ART_ROOT_IO_ROOTDB}
              ${SQLITE3}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
             )

simple_plugin(DetectorPropertiesServiceStandardGArToad "service"
              DetectorInfo
              Geometry_GeometryGAr_service
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
              ${ART_ROOT_IO_ROOTDB}
              ${SQLITE3}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
             )
simple_plugin(GArPropertiesServiceStandard "service"
              DetectorInfo
              DetectorInfo_DetectorPropertiesServiceStandardGAr_service
              art_Utilities canvas
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${CETLIB_EXCEPT}
             )

simple_plugin(ECALPropertiesServiceStandard "service"
             DetectorInfo
             DetectorInfo_DetectorPropertiesServiceStandardGAr_service
             art_Utilities canvas
             ${ART_FRAMEWORK_SERVICES_REGISTRY}
             ${MF_MESSAGELOGGER}
             ${CETLIB_EXCEPT}
             )

simple_plugin(MagneticFieldServiceGAr "service"
              DetectorInfo
              art_Utilities canvas
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${MF_MESSAGELOGGER}
              ${CETLIB_EXCEPT}
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
              ${ROOT_GEOM}
              ${ROOT_CORE}
              ${G4GEOMETRY}
              ${G4GLOBAL}
              ${ART_FRAMEWORK_SERVICES_BASIC}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_ROOT_IO_TFILE_SUPPORT}
              ${CLHEP}
              )

install_headers()
install_fhicl()
install_source()

#add_subdirectory(test)
