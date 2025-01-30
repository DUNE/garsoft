/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains functions struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5sourceidmaps to be serialized via nlohmann::json.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_NLJS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_NLJS_HPP

// My structs
#include "Structs.hpp"


#include <nlohmann/json.hpp>

namespace dunedaq::hdf5libs::hdf5sourceidmaps {

    using data_t = nlohmann::json;
    
    inline void to_json(data_t& j, const SourceID& obj) {
        j["subsys"] = obj.subsys;
        j["id"] = obj.id;
    }
    
    inline void from_json(const data_t& j, SourceID& obj) {
        if (j.contains("subsys"))
            j.at("subsys").get_to(obj.subsys);    
        if (j.contains("id"))
            j.at("id").get_to(obj.id);    
    }
    
    inline void to_json(data_t& j, const FragmentTypeSourceIDPair& obj) {
        j["fragment_type"] = obj.fragment_type;
        j["sourceids"] = obj.sourceids;
    }
    
    inline void from_json(const data_t& j, FragmentTypeSourceIDPair& obj) {
        if (j.contains("fragment_type"))
            j.at("fragment_type").get_to(obj.fragment_type);    
        if (j.contains("sourceids"))
            j.at("sourceids").get_to(obj.sourceids);    
    }
    
    inline void to_json(data_t& j, const FragmentTypeSourceIDMap& obj) {
        j["map_entries"] = obj.map_entries;
    }
    
    inline void from_json(const data_t& j, FragmentTypeSourceIDMap& obj) {
        if (j.contains("map_entries"))
            j.at("map_entries").get_to(obj.map_entries);    
    }
    
    inline void to_json(data_t& j, const SourceIDGeoIDPair& obj) {
        j["subsys"] = obj.subsys;
        j["id"] = obj.id;
        j["geoids"] = obj.geoids;
    }
    
    inline void from_json(const data_t& j, SourceIDGeoIDPair& obj) {
        if (j.contains("subsys"))
            j.at("subsys").get_to(obj.subsys);    
        if (j.contains("id"))
            j.at("id").get_to(obj.id);    
        if (j.contains("geoids"))
            j.at("geoids").get_to(obj.geoids);    
    }
    
    inline void to_json(data_t& j, const SourceIDPathPair& obj) {
        j["subsys"] = obj.subsys;
        j["id"] = obj.id;
        j["path"] = obj.path;
    }
    
    inline void from_json(const data_t& j, SourceIDPathPair& obj) {
        if (j.contains("subsys"))
            j.at("subsys").get_to(obj.subsys);    
        if (j.contains("id"))
            j.at("id").get_to(obj.id);    
        if (j.contains("path"))
            j.at("path").get_to(obj.path);    
    }
    
    inline void to_json(data_t& j, const SourceIDGeoIDMap& obj) {
        j["map_entries"] = obj.map_entries;
    }
    
    inline void from_json(const data_t& j, SourceIDGeoIDMap& obj) {
        if (j.contains("map_entries"))
            j.at("map_entries").get_to(obj.map_entries);    
    }
    
    inline void to_json(data_t& j, const SourceIDPathMap& obj) {
        j["map_entries"] = obj.map_entries;
    }
    
    inline void from_json(const data_t& j, SourceIDPathMap& obj) {
        if (j.contains("map_entries"))
            j.at("map_entries").get_to(obj.map_entries);    
    }
    
    inline void to_json(data_t& j, const SubdetectorSourceIDPair& obj) {
        j["subdetector"] = obj.subdetector;
        j["sourceids"] = obj.sourceids;
    }
    
    inline void from_json(const data_t& j, SubdetectorSourceIDPair& obj) {
        if (j.contains("subdetector"))
            j.at("subdetector").get_to(obj.subdetector);    
        if (j.contains("sourceids"))
            j.at("sourceids").get_to(obj.sourceids);    
    }
    
    inline void to_json(data_t& j, const SubdetectorSourceIDMap& obj) {
        j["map_entries"] = obj.map_entries;
    }
    
    inline void from_json(const data_t& j, SubdetectorSourceIDMap& obj) {
        if (j.contains("map_entries"))
            j.at("map_entries").get_to(obj.map_entries);    
    }
    
} // namespace dunedaq::hdf5libs::hdf5sourceidmaps

#endif // DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_NLJS_HPP
