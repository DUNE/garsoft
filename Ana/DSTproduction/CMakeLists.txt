

cet_make_exec(makeDST.exe
  SOURCE makeDST.cxx
  LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
    cetlib cetlib_except
)

install_source()
install_headers()
