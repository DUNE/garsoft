include(CetRootCint)
set( PACKAGE EventDisplay3D )
FILE( GLOB src_files *.cxx )
cet_rootcint( ${PACKAGE} )

art_make_library(SOURCE    ${src_files} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
                 LIBRARIES SimulationDataProducts
 			   RawDataProducts
	                   DetectorInfo_DetectorPropertiesServiceStandardGAr_service
                           Geometry_GeometryGAr_service
			   ReconstructionDataProducts
                           nusimdata_SimulationBase
                           nuevdb_EventDisplayBase
			   nuevdb_EventDisplayBase_Colors_service
                           ${ART_FRAMEWORK_CORE}
			   ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_IO}
			   art_Persistency_Common
                           canvas
		 	   art_Utilities
 			   ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
			   ${CETLIB}
                           ROOT::Eve
                           ${ROOT_EVE_LIB_LIST}
                           ROOT::EG
                           ${ROOT_RGL}
                           ROOT::X3d
                )

simple_plugin (EVD3D module
                  EventDisplay_EVD
                  EventDisplay_EVD3D
                  Geometry
                  nuevdb_EventDisplayBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_IO}
                  ${ART_ROOT_IO_Root_RootInput_source}
                  art_Utilities
                  canvas
				          GArRecoAlg
                  ${MF_MESSAGELOGGER}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${ROOT_GUI}
                  ${ROOT_GEOM}
                  ROOT::Eve
                  ROOT::EG
                  ${ROOT_RGL}
          )

simple_plugin (EVD3DToad module
                  EventDisplay_EVD
                  EventDisplay_EVD3D
                  Geometry
                  nuevdb_EventDisplayBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_IO}
                  ${ART_ROOT_IO_Root_RootInput_source}
                  art_Utilities
                  canvas
                                          GArRecoAlg
                  ${MF_MESSAGELOGGER}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${ROOT_GUI}
                  ${ROOT_GEOM}
                  ROOT::Eve
                  ROOT::EG
                  ${ROOT_RGL}
          )
simple_plugin (EventDisplay3DService service
                  EventDisplay_EVD3D
                  Geometry
                  nuevdb_EventDisplayBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_IO}
                  ART_ROOT_IO_RootInput_source
                  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                  art_Persistency_Common
                  art_Persistency_Provenance
                  art_Utilities
		  ART_ROOT_IO
                  canvas
                  ${MF_MESSAGELOGGER}
                  ${FHICLCPP}
                  ${CETLIB}
                    ROOT::Eve
                    ROOT::Core
                    ROOT::Hist
                    ROOT::Tree              
               )

install_headers()
install_source()
