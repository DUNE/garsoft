/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5sourceidmaps.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_STRUCTS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_STRUCTS_HPP

#include <cstdint>

#include <vector>
#include <string>

namespace dunedaq::hdf5libs::hdf5sourceidmaps {

    // @brief Reasonably sized number
    using NumericValue = uint32_t; // NOLINT


    // @brief A single SourceID
    struct SourceID 
    {

        // @brief SourceID subsystem
        NumericValue subsys = 0;

        // @brief SourceID ID
        NumericValue id = 0;
    };

    // @brief List of SourceIDs
    using SourceIDList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::SourceID>;

    // @brief A single FragmentType-to-SourceID map entry
    struct FragmentTypeSourceIDPair 
    {

        // @brief FragmentType
        NumericValue fragment_type = 0;

        // @brief List of SourceIDs for a given FragmentType
        SourceIDList sourceids = {};
    };

    // @brief List of FragmentType to SourceID map entries
    using FragmentTypeMapEntryList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::FragmentTypeSourceIDPair>;

    // @brief Information that is needed to build up the map of FragmentTypes to SourceIDs
    struct FragmentTypeSourceIDMap 
    {

        // @brief The list of entries in the map
        FragmentTypeMapEntryList map_entries = {};
    };

    // @brief 64-bit number
    using GeoIDValue = uint64_t; // NOLINT


    // @brief List of GeoIDs
    using GeoIDList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::GeoIDValue>;

    // @brief A single SourceID-to-GeoID map entry
    struct SourceIDGeoIDPair 
    {

        // @brief SourceID subsystem
        NumericValue subsys = 0;

        // @brief SourceID ID
        NumericValue id = 0;

        // @brief List of GeoIDs contained within the SourceID
        GeoIDList geoids = {};
    };

    // @brief List of SourceID to GeoID map entries
    using GeoIDMapEntryList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::SourceIDGeoIDPair>;

    // @brief A string representing an HDF5 path
    using PathString = std::string;

    // @brief A single SourceID-to-HDF5-path map entry
    struct SourceIDPathPair 
    {

        // @brief SourceID subsystem
        NumericValue subsys = 0;

        // @brief SourceID ID
        NumericValue id = 0;

        // @brief Path for the DataSet in the HDF5 file
        PathString path = "";
    };

    // @brief List of SourceID to HDF5 path map entries
    using PathMapEntryList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::SourceIDPathPair>;

    // @brief Information that is needed to build up the map of SourceIDs to GeoIDs
    struct SourceIDGeoIDMap 
    {

        // @brief The list of entries in the map
        GeoIDMapEntryList map_entries = {};
    };

    // @brief Information that is needed to build up the map of SourceIDs to HDF5 DataSet paths
    struct SourceIDPathMap 
    {

        // @brief The list of entries in the map
        PathMapEntryList map_entries = {};
    };

    // @brief A single Subdetector-to-SourceID map entry
    struct SubdetectorSourceIDPair 
    {

        // @brief DetID::Subdetector
        NumericValue subdetector = 0;

        // @brief List of SourceIDs for a given Subdetector
        SourceIDList sourceids = {};
    };

    // @brief List of Subdetector to SourceID map entries
    using SubdetectorMapEntryList = std::vector<dunedaq::hdf5libs::hdf5sourceidmaps::SubdetectorSourceIDPair>;

    // @brief Information that is needed to build up the map of Subdetectors to SourceIDs
    struct SubdetectorSourceIDMap 
    {

        // @brief The list of entries in the map
        SubdetectorMapEntryList map_entries = {};
    };

} // namespace dunedaq::hdf5libs::hdf5sourceidmaps

#endif // DUNEDAQ_HDF5LIBS_HDF5SOURCEIDMAPS_STRUCTS_HPP