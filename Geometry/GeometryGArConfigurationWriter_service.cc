/**
 * @file   Geometry/GeometryGArConfigurationWriter_service.cc
 * @brief  Service writing geometry configuration information into _art_ runs.
 * @date   March 22, 2021
 * @author Eldwan Brianne (ebrianne@fnal.gov)
 *
 * Producing services have no header.
 */

// GArSoft libraries
#include "Geometry/GeometryGAr.h"
#include "SummaryDataProducts/GeometryConfigurationInfo.h"
#include "SummaryDataProducts/RunData.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ProducingService.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <memory> // std::make_unique()


// -----------------------------------------------------------------------------
namespace gar {
    namespace geo { 
        class GeometryGArConfigurationWriter; 
    }
}

/**
 * @brief Writes geometry configuration information into _art_ runs.
 * 
 * This service is part of the mandatory version check of `gar::geo::Geometry`
 * service.
 * It does not require any special configuration, but it must be listed
 * in the configuration in order for `Geometry` to work:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * services: {
 *   
 *   GeometryGArConfigurationWriter: {}
 *   
 *   Geometry: @local::experiment_geometry
 *   
 *   # ...
 *  
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * The configuration check is described in the documentation of `gar::geo::Geometry`
 * service.
 * 
 * 
 * Produced data products
 * -----------------------
 * 
 * The service guarantees that configuration information of type
 * `sumdata::GeometryConfigurationInfo` is present into the run, accessible with
 * an input tag `GeometryGArConfigurationWriter`:
 * 
 * * if such information is already available in the run, no further information
 *   is added
 * * legacy: if there is no such information, but there is a `sumdata::RunData`
 *   data product, a reduced version of the configuration information is created
 *   from the information in that data product (the first one, if multiple are
 *   present)
 * * finally, if no information is present neither in the full 
 *   `sumdata::GeometryConfigurationInfo` form nor in the legacy
 *   `sumdata::RunData` form, information is put together based on the current
 *   configuration of the `Geometry` service.
 * 
 * 
 * Service dependencies
 * ---------------------
 * 
 * * `Geometry` service
 *   (for obtaining the current configuration to put into the event)
 * 
 */
class gar::geo::GeometryGArConfigurationWriter: public art::ProducingService {
  
    public:
  
  /// Service configuration.
  struct Config {
    // no configuration parameters at this time
  };
  
  using Parameters = art::ServiceTable<Config>;
  
  /// Constructor: gets its configuration and does nothing with it.
  GeometryGArConfigurationWriter(Parameters const&);
  
    private:
  
  /// Writes the information from the service configuration into the `run`.
  virtual void postReadRun(art::Run& run) override;
  
    private:
  
  /// Alias for the pointer to the data product object to be put into the run.
  using InfoPtr_t = std::unique_ptr<gar::sumdata::GeometryConfigurationInfo>;
  
  
  /// Loads the geometry information from the `run` (either directly or legacy).
  InfoPtr_t loadInfo(art::Run& run) const;
  
  /// Creates configuration information based on the current `Geometry` service.
  static InfoPtr_t extractInfoFromGeometry();
  
  /// Reads geometry information from the run (returns null pointer if none).
  InfoPtr_t readGeometryInformation(art::Run& run) const;
  
  /// Upgrades legacy `sumdata::RunData` in `run` to geometry information
  /// (returns null pointer if no legacy information is present).
  InfoPtr_t makeInfoFromRunData(art::Run& run) const;
  
  /// Returns a pointer to the `sumdata::RunData` in `run` (nullptr if none).
  gar::sumdata::RunData const* readRunData(art::Run& run) const;
  
  /// Converts the legacy `data` into geometry configuration information.
  static InfoPtr_t convertRunDataToGeometryInformation(gar::sumdata::RunData const& data);
  
  /// Alias to `std::make_unique<sumdata::GeometryConfigurationInfo>`.
  static InfoPtr_t makeInfoPtr(gar::sumdata::GeometryConfigurationInfo const& info) { return std::make_unique<gar::sumdata::GeometryConfigurationInfo>(info); }
  
  
}; // gar::geo::GeometryGArConfigurationWriter


// -----------------------------------------------------------------------------
// --- implementation
// -----------------------------------------------------------------------------
gar::geo::GeometryGArConfigurationWriter::GeometryGArConfigurationWriter(Parameters const&)
{
  produces<gar::sumdata::GeometryConfigurationInfo, art::InRun>();
}


// -----------------------------------------------------------------------------
void gar::geo::GeometryGArConfigurationWriter::postReadRun(art::Run& run) {
  
  InfoPtr_t confInfo = gar::geo::GeometryGArConfigurationWriter::loadInfo(run);
  
  if (!confInfo) confInfo = extractInfoFromGeometry();
  
  run.put(std::move(confInfo), art::fullRun());
  
} // gar::geo::GeometryGArConfigurationWriter::postReadRun()


// -----------------------------------------------------------------------------
auto gar::geo::GeometryGArConfigurationWriter::loadInfo(art::Run& run) const -> InfoPtr_t
{
  
  /*
   * Read geometry configuration information from the run:
   * 
   * 1. first attempt to directly read information from past runs of this
   *    service
   * 2. if none is found, attempt reading legacy information and upgrade it
   * 3. if no legacy information is found either, return a null pointer
   * 
   */
  InfoPtr_t info = readGeometryInformation(run);
  
  return info ? std::move(info) : makeInfoFromRunData(run);
  
} // gar::geo::GeometryGArConfigurationWriter::loadInfo()


// -----------------------------------------------------------------------------
auto gar::geo::GeometryGArConfigurationWriter::extractInfoFromGeometry() -> InfoPtr_t {
  
  gar::sumdata::GeometryConfigurationInfo confInfo = art::ServiceHandle<gar::geo::GeometryGAr>()->configurationInfo();

  MF_LOG_DEBUG("GeometryGArConfigurationWriter") << "Geometry configuration information from service:\n"
                                                 << confInfo;

  return makeInfoPtr(std::move(confInfo));
  
} // gar::geo::GeometryGArConfigurationWriter::extractInfoFromGeometry()


// -----------------------------------------------------------------------------
auto gar::geo::GeometryGArConfigurationWriter::readGeometryInformation(art::Run& run) const -> InfoPtr_t
{
  auto infoHandle = run.getHandle<gar::sumdata::GeometryConfigurationInfo>(art::InputTag{"GeometryGArConfigurationWriter"});
  if (infoHandle)
    {
      return makeInfoPtr(*infoHandle);
    }
  return InfoPtr_t{};

} // gar::geo::GeometryGArConfigurationWriter::hasGeometryInformation()


// -----------------------------------------------------------------------------
auto gar::geo::GeometryGArConfigurationWriter::makeInfoFromRunData(art::Run& run) const -> InfoPtr_t
{
  
  gar::sumdata::RunData const* runData = readRunData(run);
  return runData ? convertRunDataToGeometryInformation(*runData): InfoPtr_t{};
  
} // gar::geo::GeometryGArConfigurationWriter::makeInfoFromRunData()


// -----------------------------------------------------------------------------
gar::sumdata::RunData const* gar::geo::GeometryGArConfigurationWriter::readRunData(art::Run& run) const
{
  auto allRunData = run.getMany<gar::sumdata::RunData>();

  return allRunData.empty() ? nullptr: allRunData.front().product();
} // gar::geo::GeometryGArConfigurationWriter::readRunData()


// -----------------------------------------------------------------------------
auto gar::geo::GeometryGArConfigurationWriter::convertRunDataToGeometryInformation(gar::sumdata::RunData const& data) -> InfoPtr_t
{
  
  gar::sumdata::GeometryConfigurationInfo confInfo;
  
  // we use the simplest version 1 data format (legacy format)
  confInfo.dataVersion = gar::sumdata::GeometryConfigurationInfo::DataVersion_t{1};
  confInfo.detectorName = data.DetName();

  MF_LOG_DEBUG("GeometryGArConfigurationWriter")
      << "Built geometry configuration information from run data:\n"
      << confInfo;

  return makeInfoPtr(std::move(confInfo));
  
} // gar::geo::GeometryGArConfigurationWriter::convertRunDataToGeometryInformation()


// -----------------------------------------------------------------------------
DEFINE_ART_PRODUCING_SERVICE(gar::geo::GeometryGArConfigurationWriter)


// -----------------------------------------------------------------------------
