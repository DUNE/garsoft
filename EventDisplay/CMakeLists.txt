
include(CetRootCint)

set( PACKAGE EventDisplay )
FILE( GLOB src_files *.cxx )
LIST( REMOVE_ITEM src_files "Display3DPad.cxx")
LIST( REMOVE_ITEM src_files "Display3DView.cxx")

cet_rootcint( ${PACKAGE} )

art_make_library(SOURCE    ${src_files} ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE}Cint.cc
                 LIBRARIES EventDisplay_ColorDrawingOptions_service
		 	                     EventDisplay_RawDrawingOptions_service
		                       EventDisplay_RecoDrawingOptions_service
			                     EventDisplay_SimulationDrawingOptions_service
 			                     SimulationDataProducts
 			                     RawDataProducts
	                         DetectorInfo_DetectorPropertiesServiceStandard_service
                           Geometry_Geometry_service
			                     ReconstructionDataProducts
                           nusimdata_SimulationBase
                           nutools_EventDisplayBase
			                     nutools_EventDisplayBase_Colors_service
                           ${ART_FRAMEWORK_CORE}
			                     ${ART_FRAMEWORK_PRINCIPAL}
                           ${ART_FRAMEWORK_IO}
			                     art_Persistency_Common
                           canvas
		 	                     art_Utilities
 			                     ${MF_MESSAGELOGGER}
                           ${FHICLCPP}
			                     ${CETLIB}
                           ${ROOT_EVE}
                           ${ROOT_EVE_LIB_LIST}
                           ${ROOT_EG}
                           ${ROOT_RGL}
                           ${ROOT_X3d}
                )

simple_plugin (EVD module
               EventDisplay
			         nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         art_Utilities
               canvas
         			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
     )

simple_plugin (AnalysisDrawingOptions service
			         Geometry_Geometry_service
			         Geometry
			         ReconstructionDataProducts
               nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
			         ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (ColorDrawingOptions service
			         nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
		         	 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (EvdLayoutOptions service
			         nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
        			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (RawDrawingOptions service
			         nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
			         ${ART_FRAMEWORK_PRINCIPAL}
			         ${ART_FRAMEWORK_SERVICES_REGISTRY}
			         ${ART_FRAMEWORK_SERVICES_BASIC}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			         ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			         art_Persistency_Common
			         art_Persistency_Provenance
			         art_Utilities
               canvas
        			 ${MF_MESSAGELOGGER}
	             ${FHICLCPP}
			         ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (RecoDrawingOptions service
			         nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
         			 ${ART_FRAMEWORK_PRINCIPAL}
         			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
         			 ${ART_FRAMEWORK_SERVICES_BASIC}
         			 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
         			 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
         			 art_Persistency_Common
	         		 art_Persistency_Provenance
	         		 art_Utilities
               canvas
	         		 ${MF_MESSAGELOGGER}
               ${FHICLCPP}
		         	 ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (SimulationDrawingOptions service
               nutools_EventDisplayBase
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
         			 ${ART_FRAMEWORK_SERVICES_REGISTRY}
         			 ${ART_FRAMEWORK_SERVICES_BASIC}
	         		 ${ART_FRAMEWORK_SERVICES_OPTIONAL}
		         	 ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	         		 art_Persistency_Common
	         		 art_Persistency_Provenance
		         	 art_Utilities
               canvas
	         		 ${MF_MESSAGELOGGER}
               ${FHICLCPP}
		         	 ${CETLIB}
               ${ROOT_BASIC_LIB_LIST}
     )

simple_plugin (EVD3D module
                  EventDisplay
                  Geometry
                  nutools_EventDisplayBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_IO}
                  ${ART_FRAMEWORK_IO_Root_RootInput_source}
                  art_Utilities
                  canvas
				  GArRecoAlg
                  ${MF_MESSAGELOGGER}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${ROOT_GUI}
                  ${ROOT_GEOM}
                  ${ROOT_EVE}
                  ${ROOT_EG}
                  ${ROOT_RGL}
          )

simple_plugin (EventDisplay3DService service
                  EventDisplay
                  Geometry
                  nutools_EventDisplayBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_IO}
                  ${ART_FRAMEWORK_IO_Root_RootInput_source}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                  art_Persistency_Common
                  art_Persistency_Provenance
                  art_Utilities
                  art_Framework_IO_Root_RootInput_source
                  art_Framework_IO_Root
                  canvas
                  ${MF_MESSAGELOGGER}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${ROOT_EVE}
                  ${ROOT_BASIC_LIB_LIST}
               )

install_headers()
install_fhicl()
install_source()
