/**
 * @file   DetectorPropertiesStandardTestHelpers.h
 * @brief  Helpers for support of DetectorPropertiesService in GArSoft tests
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   May 6th, 2016
 * 
 * This library is a pure header.
 * It requires linkage with:
 * 
 * * `lardata_DetectorInfo`
 * * `mf_MessageLogger`
 * * `mf_Utilities`
 * * `fhiclcpp`
 * 
 */

#ifndef LARDATA_DETECTORINFO_DETECTORPROPERTIESSTANDARDTESTHELPERSTOAD_H
#define LARDATA_DETECTORINFO_DETECTORPROPERTIESSTANDARDTESTHELPERSTOAD_H 1

// GArSoft libraries
#include "DetectorInfo/DetectorPropertiesStandardToad.h"
#include "larcore/TestUtils/ProviderTestHelpers.h"

// framework and utility libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <string>
#include <memory> // std::unique_ptr<>


namespace testing {
  
  /**
   * @brief Set up a detinfo::DetectorPropertiesStandard from a parameter set
   * @param pset parameter set for the configuration
   * @param providers pack of providers needed by DetectorPropertiesService
   * @return pointer to a newly created and set up DetectorPropertiesService
   * 
   * This class specialisation enables the support of `SetupProvider()`
   * methods of `testing::TesterEnvironment`.
   * The resulting setup call will be something like:
   *     
   *     env.SetupProviderFor
   *       <detinfo::DetectorPropertiesStandard, detinfo::DetectorProperties>
   *       (pset, env.ProviderPackFor<detinfo::DetectorPropertiesStandard>());
   *     
   * (note the use of `SetupProviderFor()` to register the service also as an
   * implementation of DetectorProperties interface).
   * This implicitly assumes that whatever service providers
   * detinfo::DetectorPropertiesStandard needs, they have already been set up
   * in the environment.
   */
  template <>
  struct ProviderSetupClass<detinfo::DetectorPropertiesStandardToad> {
    
    static std::unique_ptr<detinfo::DetectorPropertiesStandardToad> setup
      (
        fhicl::ParameterSet const& pset,
        detinfo::DetectorPropertiesStandardToad::providers_type const& providers
      )
      {
        // some feedback about whether we are using the right configuration
        std::string ServiceProviderPath;
        if (pset.get_if_present("service_provider", ServiceProviderPath)) {
          std::string ServiceProviderName = ServiceProviderPath;
          size_t iSlash = ServiceProviderPath.rfind('/');
          if (iSlash != std::string::npos)
            ServiceProviderName.erase(0, iSlash + 1);
          
          if (ServiceProviderName == "DetectorPropertiesServiceStandardGAr") {
            LOG_TRACE("ProviderSetup") << "Verified service implementation for "
              "DetectorPropertiesService"
              ": '" << ServiceProviderPath << "'";
          }
          else {
            // this means you should not be using the simple set up...
            mf::LogWarning("setupProvider")
              << "This set up is for a DetectorPropertiesStandard provider.\n"
              "Your configuration specifies a '" << ServiceProviderPath
              << "' service implementation"
              " that is not known to use that provider.";
          }
        }
    
        /// @todo remove this debug stuff!
        // 
        // create the new DetectorPropertiesStandard service provider
        //
        // we choose to ignore a configuration parameter that is known to be used
        // by the art services but unknown to the provider; in this way we can use
        // for this test the same configuration as for art-based tests.
        //
        std::set<std::string> ignore_keys({ "InheritNumberTimeSamples" });
        {
           mf::LogInfo debug("setupProvider<DetectorPropertiesStandard>");
           
           debug << "Asked to ignore " << ignore_keys.size() << " keys:";
           for (auto const& key: ignore_keys) debug << " '" << key << "'";
        }
        
        return std::make_unique<detinfo::DetectorPropertiesStandardToad>
          (pset, providers, ignore_keys);
      } // setup()
      
  }; // ProviderSetupClass<DetectorPropertiesStandard>
  

  /**
   * @brief Environment setup helper for DetectorPropertiesStandard
   * @tparam TestEnv type of environment to set up
   * @see simpleEnvironmentSetup()
   * 
   * A service provider is set up in the environment, associated with the types
   * detinfo::DetectorPropertiesStandard and detinfo::DetectorProperties.
   * Its configuration is read from "services.DetectorPropertiesService".
   * 
   * The environment is expected to expose an interface equivalent to the one
   * of `testing::TesterEnvironment`.
   * 
   * This class specialisation enables the support of `SimpleProviderSetup()`
   * methods of `testing::TesterEnvironment`.
   * It should be possible to set up a testing environment by calling:
   *     
   *     env.SimpleProviderSetup<detinfo::DetectorPropertiesStandard>();
   *     
   * The provider will be available from any of these two calls:
   *     
   *     env.Provider<detinfo::DetectorPropertiesStandard>();
   *     env.Provider<detinfo::DetectorProperties>();
   *     
   * @note This function assumes there is a `DetectorPropertiesService` service
   * configuration available in the environment, and that the service providers
   * geo::GeometryCore, detinfo::GArProperties and detinfo::DetectorClocks are
   * configured and available.
   */
  template <typename TestEnv>
  struct SimpleEnvironmentSetupClass
    <detinfo::DetectorPropertiesStandardToad, TestEnv>
  {
    static detinfo::DetectorPropertiesStandardToad* setup(TestEnv& env)
      {
        return env.template SetupProviderFor
          <detinfo::DetectorProperties, detinfo::DetectorPropertiesStandardToad>
          (
            env.ServiceParameters("DetectorPropertiesService"),
            env.template ProviderPackFor<detinfo::DetectorPropertiesStandardToad>()
          );
      } // setup()
    
  }; // SimpleEnvironmentSetupClass<detinfo::DetectorPropertiesStandard>
  

} // namespace testing

  
#endif // LARDATA_DETECTORINFO_DETECTORPROPERTIESSTANDARDTESTHELPERS_H
