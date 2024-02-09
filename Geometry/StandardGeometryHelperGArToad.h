/**
 * @file   StandardGeometryHelperGAr.h
 * @brief  Geometry helper service for detectors with strictly standard mapping
 * @author rs@fnal.gov
 *
 * Handles detector-specific information for the generic Geometry service
 * within GArSoft. Derived from the ExptGeoHelperInterface class. This version
 * provides strictly standard functionality
 */

#ifndef GEO_StandardGeometryHelperGArToad_h
#define GEO_StandardGeometryHelperGArToad_h

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
    class StandardGeometryHelperGArToad : public ExptGeoHelperInterface
    {
    public:

        /// Constructor; follows the standard art service signature
      StandardGeometryHelperGArToad
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
      virtual SegmentationAlgPtr_t doGetECALSegmentationAlg() const override;

      virtual void doConfigureMinervaSegmentationAlg
      (fhicl::ParameterSet const& segParameters, geo::GeometryCore* geom)
      override;
      virtual SegmentationAlgPtr_t doGetMinervaSegmentationAlg() const override;

      virtual void doConfigureMuIDSegmentationAlg
      (fhicl::ParameterSet const& segParameters, geo::GeometryCore* geom)
      override;
      virtual SegmentationAlgPtr_t doGetMuIDSegmentationAlg() const override;

      fhicl::ParameterSet fPset; ///< copy of configuration parameter set
      std::shared_ptr<geo::seg::ChannelMapAlg> fChannelMap; ///< channel map algorithm
      std::shared_ptr<geo::seg::SegmentationAlg> fECALSegmentationAlg; ///< ECAL Segmentation Alg
      std::shared_ptr<geo::seg::SegmentationAlg> fMinervaSegmentationAlg; ///< Tracker Sc Segmentation Alg
      std::shared_ptr<geo::seg::SegmentationAlg> fMuIDSegmentationAlg; ///< MuID Segmentation Alg

    };

  }
} // end gar

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::geo::StandardGeometryHelperGArToad, gar::geo::ExptGeoHelperInterface, LEGACY)

#endif // GEO_StandardGeometryHelperGAr_h
