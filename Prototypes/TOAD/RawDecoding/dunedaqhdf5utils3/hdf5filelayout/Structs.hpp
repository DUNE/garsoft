/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5filelayout.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_STRUCTS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_STRUCTS_HPP

#include <cstdint>

#include <vector>
#include <string>

namespace dunedaq::hdf5libs::hdf5filelayout {

    // @brief A count of not too many things
    using Count = int32_t;


    // @brief A float number of 4 bytes
    using Factor = float;


    // @brief A string used in the hdf5 configuration
    using HDFString = std::string;

    // @brief Parameters for the HDF5 Group and DataSet names
    struct PathParams 
    {

        // @brief The special keyword that identifies this entry (e.g. "TPC", "PDS", "TPC_TP", etc.)
        HDFString detector_group_type = "unspecified";

        // @brief The detector name that should be used inside the HDF5 file
        HDFString detector_group_name = "unspecified";

        // @brief Prefix for the element name
        HDFString element_name_prefix = "Element";

        // @brief Number of digits to use for the element ID number inside the HDF5 file
        Count digits_for_element_number = 5;
    };

    // @brief List of subdetector path parameters
    using PathParamList = std::vector<dunedaq::hdf5libs::hdf5filelayout::PathParams>;

    // @brief Parameters for the layout of Groups and DataSets within the HDF5 file
    struct FileLayoutParams 
    {

        // @brief Prefix for the record name
        HDFString record_name_prefix = "TriggerRecord";

        // @brief Number of digits to use for the record number in the record name inside the HDF5 file
        Count digits_for_record_number = 6;

        // @brief Number of digits to use for the sequence number in the TriggerRecord name inside the HDF5 file
        Count digits_for_sequence_number = 4;

        // @brief Dataset name for the record header
        HDFString record_header_dataset_name = "TriggerRecordHeader";

        // @brief Group name to use for raw data
        HDFString raw_data_group_name = "RawData";

        // @brief Group name to use for views of the raw data
        HDFString view_group_name = "Views";

        // @brief 
        PathParamList path_param_list = {};
    };

    // @brief Parameter that can be used to enable or disable functionality
    using Flag = bool;

    // @brief A count of very many things
    using Size = uint64_t; // NOLINT


} // namespace dunedaq::hdf5libs::hdf5filelayout

#endif // DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_STRUCTS_HPP