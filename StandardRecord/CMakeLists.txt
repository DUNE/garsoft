# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

add_subdirectory(Proxy)

set( PACKAGE garsoft_StandardRecord )
FILE( GLOB src_files *.cxx )

cet_make_library( LIBRARY_NAME garsoft_StandardRecord
                  SOURCE       ${src_files}
                  LIBRARIES    ${ROOT_BASIC_LIB_LIST}
                )

build_dictionary( garsoft_StandardRecord
                  DICTIONARY_LIBRARIES garsoft_StandardRecord
                )

install_headers()
install_source()
