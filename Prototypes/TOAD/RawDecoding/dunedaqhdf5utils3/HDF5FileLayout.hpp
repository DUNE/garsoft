/**
 * @file HDF5FileLayout.hpp
 *
 * FileLayoutObject that is used to describe/provide instructions for
 * organizing DUNE DAQ HDF5 files.
 *
 * WK 02.02.2022 -- Future could make a FileLayout descriptor per SystemType
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 * 
 * Modified April 2024 Tom Junk -- update from dunedaqhdf5utils2
 */

#ifndef HDF5LIBS_INCLUDE_HDF5LIBS_HDF5FILELAYOUT_HPP_
#define HDF5LIBS_INCLUDE_HDF5LIBS_HDF5FILELAYOUT_HPP_

#include "hdf5filelayout/Structs.hpp"

#include "daqdataformats/v4_4_0/Fragment.hpp"
#include "daqdataformats/v4_4_0/SourceID.hpp"
#include "daqdataformats/v4_4_0/TimeSliceHeader.hpp"
#include "daqdataformats/v4_4_0/TriggerRecordHeader.hpp"

#include "nlohmann/json.hpp"

#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// Offline includes

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

namespace dunedaq {

namespace hdf5libs {

class HDF5FileLayout
{
public:
  /**
   * @brief Constructor from json conf, used in DataWriter. Version always most recent.
   */
  explicit HDF5FileLayout(hdf5filelayout::FileLayoutParams conf, uint32_t version = 5); // NOLINT(build/unsigned)

  uint32_t get_version() const noexcept // NOLINT(build/unsigned)
  {
    return m_version;
  }

  std::string get_record_name_prefix() const noexcept { return m_conf_params.record_name_prefix; }

  int get_digits_for_record_number() const noexcept { return m_conf_params.digits_for_record_number; }

  int get_digits_for_sequence_number() const noexcept { return m_conf_params.digits_for_sequence_number; }

  std::string get_record_header_dataset_name() const noexcept { return m_conf_params.record_header_dataset_name; }

  std::map<daqdataformats::SourceID::Subsystem, hdf5filelayout::PathParams> get_path_params_map() const
  {
    return m_path_params_map;
  }

  hdf5filelayout::PathParams get_path_params(daqdataformats::SourceID::Subsystem type) const;

  hdf5filelayout::FileLayoutParams get_file_layout_params() const { return m_conf_params; }

  /**
   * @brief get string for record number
   */
  std::string get_record_number_string(uint64_t record_number, // NOLINT(build/unsigned)
                                       daqdataformats::sequence_number_t seq_num = 0) const;

  /**
   * @brief get string for Trigger number
   */
  std::string get_trigger_number_string(daqdataformats::trigger_number_t trig_num,
                                        daqdataformats::sequence_number_t seq_num = 0) const;

  /**
   * @brief get string for TimeSlice number
   */
  std::string get_timeslice_number_string(daqdataformats::timeslice_number_t ts_num) const;

  /**
   * @brief get the correct path for the TriggerRecordHeader
   */
  std::vector<std::string> get_path_elements(const daqdataformats::TriggerRecordHeader& trh) const;

  /**
   * @brief get the correct path for the TimeSliceHeader
   */
  std::vector<std::string> get_path_elements(const daqdataformats::TimeSliceHeader& tsh) const;

  /**
   * @brief get the correct path for the Fragment
   */
  std::vector<std::string> get_path_elements(const daqdataformats::FragmentHeader& fh) const;

  /**
   * @brief extract Fragment GeoID given path elements
   */
  daqdataformats::SourceID get_source_id_from_path_elements(std::vector<std::string> const& path_elements) const;

  /**
   * @brief get the full path for a record header dataset based on trig/seq number
   */
  std::string get_record_header_path(uint64_t rec_num, // NOLINT (build/unsigned)
                                     daqdataformats::sequence_number_t seq_num = 0) const;

  /**
   * @brief get the full path for a TriggerRecordHeader dataset based on trig/seq number
   */
  std::string get_trigger_record_header_path(daqdataformats::trigger_number_t trig_num,
                                             daqdataformats::sequence_number_t seq_num = 0) const;

  /**
   * @brief get the full path for a TimesliceHeader dataset based on ts number
   */
  std::string get_timeslice_header_path(daqdataformats::timeslice_number_t ts_num) const;

  /**
   * @brief get the full path for a Fragment dataset based on trig/seq number and element ID
   */
  std::string get_fragment_path(uint64_t trig_num, // NOLINT(build/unsigned)
                                daqdataformats::sequence_number_t seq_num,
                                daqdataformats::SourceID element_id) const;
  /**
   * @brief get the full path for a Fragment dataset based on trig/seq number, give element_id pieces
   */
  std::string get_fragment_path(uint64_t trig_num, // NOLINT(build/unsigned)
                                daqdataformats::sequence_number_t seq_num,
                                daqdataformats::SourceID::Subsystem type,
                                uint32_t element_id) const; // NOLINT(build/unsigned)

  /**
   * @brief get the full path for a Fragment dataset based on trig/seq number, give element_id pieces
   */
  std::string get_fragment_path(uint64_t trig_num, // NOLINT(build/unsigned)
                                daqdataformats::sequence_number_t seq_num,
                                const std::string& typestring,
                                uint32_t element_id) const; // NOLINT(build/unsigned)

  /**
   * @brief get the path for a Fragment type group based on trig/seq number and type
   */
  std::string get_fragment_type_path(uint64_t trig_num, // NOLINT(build/unsigned)
                                     daqdataformats::sequence_number_t seq_num,
                                     daqdataformats::SourceID::Subsystem type) const;

  /**
   * @brief get the path for a Fragment type group based on trig/seq number and type
   */
  std::string get_fragment_type_path(uint64_t trig_num, // NOLINT(build/unsigned)
                                     daqdataformats::sequence_number_t seq_num,
                                     std::string typestring) const;

private:
  /**
   * @brief FileLayout configuration parameters
   */
  hdf5filelayout::FileLayoutParams m_conf_params;

  /**
   * @brief version number
   */
  uint32_t m_version; // NOLINT(build/unsigned)

  /**
   * @brief map translation for SourceID::Subsystem to dataset path parameters
   */
  std::map<daqdataformats::SourceID::Subsystem, hdf5filelayout::PathParams> m_path_params_map;

  /**
   * @brief map translation for the detector group name to SourceID::Subsystem
   */
  std::map<std::string, daqdataformats::SourceID::Subsystem> m_detector_group_name_to_type_map;

  // quick powers of ten lookup
  constexpr static uint64_t m_powers_ten[] // NOLINT(build/unsigned)
    = {
        1,                    // 1e0
        10,                   // 1e1
        100,                  // 1e2
        1000,                 // 1e3
        10000,                // 1e4
        100000,               // 1e5
        1000000,              // 1e6
        10000000,             // 1e7
        100000000,            // 1e8
        1000000000,           // 1e9
        10000000000,          // 1e10
        100000000000,         // 1e11
        1000000000000,        // 1e12
        10000000000000,       // 1e13
        100000000000000,      // 1e14
        1000000000000000,     // 1e15
        10000000000000000,    // 1e16
        100000000000000000,   // 1e17
        1000000000000000000,  // 1e18
        10000000000000000000u // 1e19
      };

  /**
   * @brief Fill path parameters maps from FileLayoutParams
   */
  void fill_path_params_maps(hdf5filelayout::FileLayoutParams const& flp);

  /**
   * @brief Version0 FileLayout parameters, for backward compatibility
   */
  hdf5filelayout::FileLayoutParams get_v0_file_layout_params();

  /**
   * @brief Check configuration for any errors.
   */
  void check_config();
};

} // namespace hdf5libs
} // namespace dunedaq

#endif // HDF5LIBS_INCLUDE_HDF5LIBS_HDF5FILELAYOUT_HPP_
