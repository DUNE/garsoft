simple_plugin(TOADInput "source"
  RawDataProducts
  Prototypes_TOAD_ChannelMap_TOADChannelMapService_service
  garsoft::dunedaqhdf5utils2
  art::Framework_Services_Registry
)

simple_plugin(TOADInput3 "source"
  RawDataProducts
  Prototypes_TOAD_ChannelMap_TOADChannelMapService_service
  garsoft::dunedaqhdf5utils3
  art::Framework_Services_Registry
)

add_subdirectory(dunedaqhdf5utils2)
add_subdirectory(dunedaqhdf5utils3)
install_headers()
install_fhicl()
install_source()
