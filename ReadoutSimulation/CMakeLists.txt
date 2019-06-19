include_directories( $ENV{GEANT4_FQ_DIR}/include )
include_directories( $ENV{NUTOOLS_DIR}/include   )
include_directories( $ENV{CLHEP_INCLUDE_DIR}     )

set(G4LIBS ${G4FR}
           ${G4GMOCREN}
           ${G4OPENGL}
           ${G4RAYTRACER}
           ${G4TREE}
           ${G4VRML}
           ${G4ANALYSIS}
           ${G4DIGITS_HITS}
           ${G4ERROR_PROPAGATION}
           ${G4EVENT}
           ${G4GEOMETRY}
           ${G4GL2PS}
           ${G4GLOBAL}
           ${G4GRAPHICS_REPS}
           ${G4INTERCOMS}
           ${G4INTERFACES}
           ${G4MATERIALS}
           ${G4MODELING}
           ${G4PARMODELS}
           ${G4PARTICLES}
           ${G4PERSISTENCY}
           ${G4PHYSICSLISTS}
           ${G4PROCESSES}
           ${G4READOUT}
           ${G4RUN}
           ${G4ERROR_PROPAGATION}
           ${G4EVENT}
           ${G4GEOMETRY}
           ${G4GL2PS}
           ${G4GLOBAL}
           ${G4GRAPHICS_REPS}
           ${G4TRACK}
           ${G4TRACKING}
           ${G4VISHEPREP}
           ${G4VISXXX}
           ${G4VIS_MANAGEMENT}
           ${G4ZLIB}
          )

art_make(LIB_LIBRARIES nusimdata_SimulationBase
                       Geometry
                       GArG4
                       SimulationDataProducts
                       RawDataProducts
                       GArRecoAlg
                       ${ART_FRAMEWORK_PRINCIPAL}
                       ${ART_FRAMEWORK_SERVICES_REGISTRY}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                       ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                       art_Persistency_Common
                       art_Persistency_Provenance
                       art_Utilities
                       canvas
                       ${MF_MESSAGELOGGER}
                       ${FHICLCPP}
                       ${CETLIB_EXCEPT}
                       ${ROOT_BASIC_LIB_LIST}
                       ${CLHEP}
                       ${G4LIBS}
                       ${Boost_SYSTEM_LIBRARY}
      MODULE_LIBRARIES Geometry
                       GArG4
                       ReadoutSimulation
                       SimulationDataProducts
                       RawDataProducts
                       GArRecoAlg
                       MCCheater_BackTracker_service
                       nusimdata_SimulationBase
                       nutools_RandomUtils_NuRandomService_service
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
