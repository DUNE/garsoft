/**
 * @file   StandardGeometryHelper.h
 * @brief  Geometry helper service for detectors with strictly standard mapping
 * @author rs@fnal.gov
 *
 * Handles detector-specific information for the generic Geometry service
 * within GArSoft. Derived from the ExptGeoHelperInterface class. This version
 * provides strictly standard functionality
 */

#ifndef GEO_StandardGeometryHelper_h
#define GEO_StandardGeometryHelper_h

// GArSoft libraries
#include "Geometry/ExptGeoHelperInterface.h"

// C/C++ standard libraries
#include <memory> // std::shared_ptr<>

// Declaration
//
namespace gar
{
  namespace geo
  {
    /**
     * @brief Simple implementation of channel mapping
     *
     * This ExptGeoHelperInterface implementation serves a ChannelMapStandardAlg
     * for experiments that are known to work well with it.
     */
    class StandardGeometryHelper : public ExptGeoHelperInterface
    {
    public:

        /// Constructor; follows the standard art service signature
      StandardGeometryHelper
      ( fhicl::ParameterSet const & pset, ::art::ActivityRegistry &reg );

      /*
       Public interface for ExptGeoHelperInterface (for reference purposes)

       Configure, initialize and return the channel map:

       void ConfigureChannelMapAlg
       (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom);

       Returns null pointer if the initialization failed:

       ChannelMapAlgPtr_t GetChannelMapAlg() const;
       */

    private:

      virtual void doConfigureChannelMapAlg
      (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
      override;
      virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;

      virtual void doConfigureECALSegmentationAlg
      (fhicl::ParameterSet const& segParameters, geo::GeometryCore* geom)
      override;
      virtual ECALSegmentationAlgPtr_t doGetECALSegmentationAlg() const override;

      fhicl::ParameterSet fPset; ///< copy of configuration parameter set
      std::shared_ptr<geo::ChannelMapAlg> fChannelMap; ///< channel map algorithm
      std::shared_ptr<geo::ECALSegmentationAlg> fECALSegmentationAlg; ///< ECAL Segmentation Alg

    };

  }
} // end gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::geo::StandardGeometryHelper, gar::geo::ExptGeoHelperInterface, LEGACY)

#endif // GEO_StandardGeometryHelper_h
