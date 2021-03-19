SET(INC_DIR $ENV{GARSOFT_DIR})

FILE(GLOB SR_DEPENDENCIES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${INC_DIR}/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT SRProxy.cxx SRProxy.h FwdDeclare.h
                   COMMAND gen_srproxy -i StandardRecord/StandardRecord.h -o SRProxy --target caf::StandardRecord --include-path ${INC_DIR}:$ENV{ROOT_INC} -op StandardRecord/Proxy/
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Proxy/
                   COMMAND cp SRProxy.h FwdDeclare.h $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Proxy/
  )

include_directories($ENV{SRPROXY_INC})

# This is a very picky error buried inside template instantiations
#add_definitions(-Wno-int-in-bool-context)

cet_make_library(LIBRARY_NAME StandardRecordProxy
                 SOURCE       SRProxy.cxx Instantiations.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

# install(FILES SRProxy.h FwdDeclare.h
#         DESTINATION $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Proxy/)
