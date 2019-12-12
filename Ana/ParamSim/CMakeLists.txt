

cet_make_library(LIBRARY_NAME garsoft_paramsim_utils
                 SOURCE
                   CAF.cpp
                   Utils.cpp
		   MCP_Skimmer.cpp
		 LIBRARIES
                   ${ROOT_BASIC_LIB_LIST}
                   ${ROOT_EG}
                 )

cet_make_exec(cafanatree_module
  SOURCE 
    cafanatree_module.cxx
  LIBRARIES
    garsoft_paramsim_utils
    ${ROOT_BASIC_LIB_LIST}
)

cet_make_exec(mcp_skimmer
  SOURCE 
    mcp_skimmer.cxx
  LIBRARIES
    garsoft_paramsim_utils
    ${ROOT_BASIC_LIB_LIST}
)

install_source()
install_headers()
