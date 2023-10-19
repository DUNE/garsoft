/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5rawdatafile to be serialized via nlohmann::json.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_NLJS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_NLJS_HPP

// My structs
#include "hdf5libs/hdf5rawdatafile/Structs.hpp"


#include <nlohmann/json.hpp>

namespace dunedaq::hdf5libs::hdf5rawdatafile {

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const GeoID& obj) {
        j["det_id"] = obj.det_id;
        j["crate_id"] = obj.crate_id;
        j["slot_id"] = obj.slot_id;
        j["stream_id"] = obj.stream_id;
    }
    
    inline void from_json(const data_t& j, GeoID& obj) {
        if (j.contains("det_id"))
            j.at("det_id").get_to(obj.det_id);    
        if (j.contains("crate_id"))
            j.at("crate_id").get_to(obj.crate_id);    
        if (j.contains("slot_id"))
            j.at("slot_id").get_to(obj.slot_id);    
        if (j.contains("stream_id"))
            j.at("stream_id").get_to(obj.stream_id);    
    }
    
    inline void to_json(data_t& j, const SrcIDGeoIDEntry& obj) {
        j["src_id"] = obj.src_id;
        j["geo_id"] = obj.geo_id;
    }
    
    inline void from_json(const data_t& j, SrcIDGeoIDEntry& obj) {
        if (j.contains("src_id"))
            j.at("src_id").get_to(obj.src_id);    
        if (j.contains("geo_id"))
            j.at("geo_id").get_to(obj.geo_id);    
    }
    
} // namespace dunedaq::hdf5libs::hdf5rawdatafile

#endif // DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_NLJS_HPP