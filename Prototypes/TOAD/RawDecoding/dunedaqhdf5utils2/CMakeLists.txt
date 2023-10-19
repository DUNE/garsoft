art_make(
	 LIBRARY_NAME dunedaqhdf5utils2
         LIB_LIBRARIES
	 cetlib::cetlib
         cetlib_except::cetlib_except
         messagefacility::MF_MessageLogger
         HDF5::HDF5
         pthread
         z
)

install_headers()
install_source()

add_subdirectory(hdf5filelayout)
add_subdirectory(hdf5sourceidmaps)
add_subdirectory(hdf5rawdatafile)

