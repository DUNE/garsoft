/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5filelayout to be serialized via nlohmann::json.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_NLJS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_NLJS_HPP

// My structs
#include "Structs.hpp"


#include <nlohmann/json.hpp>

namespace dunedaq::hdf5libs::hdf5filelayout {

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const PathParams& obj) {
        j["detector_group_type"] = obj.detector_group_type;
        j["detector_group_name"] = obj.detector_group_name;
        j["element_name_prefix"] = obj.element_name_prefix;
        j["digits_for_element_number"] = obj.digits_for_element_number;
    }
    
    inline void from_json(const data_t& j, PathParams& obj) {
        if (j.contains("detector_group_type"))
            j.at("detector_group_type").get_to(obj.detector_group_type);    
        if (j.contains("detector_group_name"))
            j.at("detector_group_name").get_to(obj.detector_group_name);    
        if (j.contains("element_name_prefix"))
            j.at("element_name_prefix").get_to(obj.element_name_prefix);    
        if (j.contains("digits_for_element_number"))
            j.at("digits_for_element_number").get_to(obj.digits_for_element_number);    
    }
    
    inline void to_json(data_t& j, const FileLayoutParams& obj) {
        j["record_name_prefix"] = obj.record_name_prefix;
        j["digits_for_record_number"] = obj.digits_for_record_number;
        j["digits_for_sequence_number"] = obj.digits_for_sequence_number;
        j["record_header_dataset_name"] = obj.record_header_dataset_name;
        j["raw_data_group_name"] = obj.raw_data_group_name;
        j["view_group_name"] = obj.view_group_name;
        j["path_param_list"] = obj.path_param_list;
    }
    
    inline void from_json(const data_t& j, FileLayoutParams& obj) {
        if (j.contains("record_name_prefix"))
            j.at("record_name_prefix").get_to(obj.record_name_prefix);    
        if (j.contains("digits_for_record_number"))
            j.at("digits_for_record_number").get_to(obj.digits_for_record_number);    
        if (j.contains("digits_for_sequence_number"))
            j.at("digits_for_sequence_number").get_to(obj.digits_for_sequence_number);    
        if (j.contains("record_header_dataset_name"))
            j.at("record_header_dataset_name").get_to(obj.record_header_dataset_name);    
        if (j.contains("raw_data_group_name"))
            j.at("raw_data_group_name").get_to(obj.raw_data_group_name);    
        if (j.contains("view_group_name"))
            j.at("view_group_name").get_to(obj.view_group_name);    
        if (j.contains("path_param_list"))
            j.at("path_param_list").get_to(obj.path_param_list);    
    }
    
} // namespace dunedaq::hdf5libs::hdf5filelayout

#endif // DUNEDAQ_HDF5LIBS_HDF5FILELAYOUT_NLJS_HPP
