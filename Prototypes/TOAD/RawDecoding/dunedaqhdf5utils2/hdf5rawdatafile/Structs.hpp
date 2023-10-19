/*
 * This file is 100% generated.  Any manual edits will likely be lost.
 *
 * This contains struct and other type definitions for shema in 
 * namespace dunedaq::hdf5libs::hdf5rawdatafile.
 */
#ifndef DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_STRUCTS_HPP
#define DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_STRUCTS_HPP

#include <cstdint>

#include <vector>

namespace dunedaq::hdf5libs::hdf5rawdatafile {

    // @brief A count of not too many things
    using Count = int32_t;


    // @brief Parameter that can be used to enable or disable functionality
    using Flag = bool;

    // @brief Geographic ID structure
    struct GeoID 
    {

        // @brief Detector ID
        Count det_id = 0;

        // @brief Crate ID
        Count crate_id = 0;

        // @brief Slot ID
        Count slot_id = 0;

        // @brief Stream ID
        Count stream_id = 0;
    };

    // @brief A count of very many things
    using Size = uint64_t; // NOLINT


    // @brief SourceID GeoID Map entry
    struct SrcIDGeoIDEntry 
    {

        // @brief Source ID
        Size src_id = 0;

        // @brief Geo ID
        GeoID geo_id = {0, 0, 0, 0};
    };

    // @brief SourceID to GeoID map
    using SrcIDGeoIDMap = std::vector<dunedaq::hdf5libs::hdf5rawdatafile::SrcIDGeoIDEntry>;

} // namespace dunedaq::hdf5libs::hdf5rawdatafile

#endif // DUNEDAQ_HDF5LIBS_HDF5RAWDATAFILE_STRUCTS_HPP