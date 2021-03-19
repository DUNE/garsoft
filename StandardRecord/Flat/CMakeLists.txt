SET(INC_DIR $ENV{GARSOFT_DIR})

FILE(GLOB SR_DEPENDENCIES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${INC_DIR}/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT FlatRecord.cxx FlatRecord.h FwdDeclare.h
                   COMMAND gen_srproxy --flat -i StandardRecord/StandardRecord.h -o FlatRecord --target caf::StandardRecord --include-path ${INC_DIR}:$ENV{ROOT_INC} -op StandardRecord/Flat/
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Flat/
                   COMMAND echo cp FlatRecord.h FwdDeclare.h $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Flat/
                   COMMAND cp FlatRecord.h FwdDeclare.h $ENV{MRB_INSTALL}/${${product}_inc_dir}/StandardRecord/Flat/
  )

include_directories($ENV{SRPROXY_INC})

# This is a very picky error buried inside template instantiations
#add_definitions(-Wno-int-in-bool-context)

cet_make_library(LIBRARY_NAME FlatRecord
                 SOURCE       FlatRecord.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

cet_make_exec(flatten_caf
              SOURCE flatten_caf.cc
              LIBRARIES FlatRecord garsoft_StandardRecord_dict # CAFAnaCore
                        ${ROOT_BASIC_LIB_LIST})
