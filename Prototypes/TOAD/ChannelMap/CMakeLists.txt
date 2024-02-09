set(PACKAGE Prototypes)
file(GLOB channel_map_file *ChannelMap*.txt)
file(GLOB src_files *.cxx)

art_make_library(LIBRARY_NAME Prototypes_TOADChannelMap
                 SOURCE ${src_files}
                 LIBRARIES art::Framework_Services_Registry
                           art::Framework_Principal
                           art_Framework_Core
                           art::Persistency_Provenance
                           ${MF_MESSAGELOGGER}
                           ROOT::Core
)

install_fhicl()
install_headers()
install_source()

install_fw( LIST ${channel_map_file} )


simple_plugin(TOADChannelMapService "service"
  art::Framework_Services_Registry
  Prototypes_TOADChannelMap
  ${MF_MESSAGELOGGER}
)

