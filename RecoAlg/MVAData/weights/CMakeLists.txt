# install files
file(GLOB xml_file *.xml)
install( FILES ${xml_file} DESTINATION ${${product}_MVAData_dir}/weights )