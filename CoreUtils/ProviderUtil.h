/**
 * @file   ProviderUtil.h
 * @brief  Simple utilities for service providers
 * @date   November 30, 2015
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * Currently this header is written so that it does not need a implementation
 * file.
 */

#ifndef DETECTORINFO_PROVIDERUTIL_H
#define DETECTORINFO_PROVIDERUTIL_H 1

// C/C++ standard libraries
#include <set>
#include <string>

namespace gar {

  /** **************************************************************************
   * @brief Returns a list of configuration keys that providers should ignore
   * @return a reference to a key list
   *
   * This function may be used for parameter validation, like in:
   *
   *     fhicl::Table<Config> cfg { pset, gar::IgnorableProviderConfigKeys() };
   *
   * where `pset` is a `fhicl::ParameterSet`. This will inform `cfg` that some
   * keys can be unexpectedly present, or missing.
   *
   * This implementation includes:
   * * art framework service keywords
   */
  inline std::set<std::string> const& IgnorableProviderConfigKeys()
  {
    static std::set<std::string> const ignorable{
      "service_type",    // added by art: service name (possibly interface)
      "service_provider" // art: service implementation name
    };
    return ignorable;
  } // IgnorableProviderConfigKeys()

} // namespace gar

#endif // DETECTORINFO_PROVIDERUTIL_H
