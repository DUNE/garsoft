////////////////////////////////////////////////////////////////////////////////
/// \file  ExptGeoHelperInterface.h
/// \brief Interface to a service that handles any experiment-specific knowledge
///        that is needed by the Geometry service.
///
///  This is an interface to a service that virtualizes detector or experiment-specific
///  knowledge that is required by the Geometry service. Experiments implement the
///  private virtual functions within a concrete service provider class to perform
///  the specified actions as appropriate for the particular experiment. It is
///  expected that such requests will occur infrequently within a job. Calculations
///  that occur frequently should be handled via interfaces that are passed
///  back to the Geometry service.
///
///  Note that the public interface for this service cannot be overriden. The
///  experiment-specific sub-classes should implement only the private methods
///  without promoting their visibility.
///
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////


#ifndef GEO_ExptGeoHelperInterface_h
#define GEO_ExptGeoHelperInterface_h


// framework libraries
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"

// C/C++ standard libraries
#include <memory> // std::shared_ptr<>
#include <vector>


// prototypes of geometry classes
namespace gar {

    namespace geo
    {

        class GeometryCore;
        namespace seg {
            class ChannelMapAlg;
            class SegmentationAlg;
        }

        /**
        * @brief Interface to a service with detector-specific geometry knowledge
        *
        * This is an interface to a service that virtualizes detector or
        * experiment-specific knowledge that is required by the Geometry service.
        * Experiments implement the private virtual functions within a concrete
        * service provider class to perform the specified actions as appropriate for
        * the particular experiment.
        * It is expected that such requests will occur infrequently within a job.
        * Calculations that occur frequently should be handled via interfaces that
        * are passed back to the Geometry service.
        *
        * @note The public interface for this service cannot be overriden.
        * The experiment-specific sub-classes should implement only the private
        * methods without promoting their visibility.
        */
        class ExptGeoHelperInterface
        {
        public:
            using ChannelMapAlgPtr_t = std::shared_ptr<const seg::ChannelMapAlg>;
            using SegmentationAlgPtr_t = std::shared_ptr<const seg::SegmentationAlg>;

            /// Virtual destructor; does nothing
            virtual ~ExptGeoHelperInterface() = default;

            /**
            * @brief Configure and initialize the channel map
            * @param sortingParameters parameters for the channel map algorithm
            * @param geom pointer to a geometry description object
            * @return a (shared) pointer to the channel mapping algorithm
            *
            * This method creates a new ChannelMapAlg according to the geometry and
            * specified configuration, then it configures the geometry itself
            * according to the channel map (usually, it resorts the data).
            */
            void ConfigureChannelMapAlg
            (fhicl::ParameterSet const & sortingParameters, geo::GeometryCore* geom);

            /// Returns null pointer if the initialization failed
            /// NOTE:  the sub-class owns the ChannelMapAlg object
            ///
            ChannelMapAlgPtr_t GetChannelMapAlg() const;

            void ConfigureECALSegmentationAlg
            (fhicl::ParameterSet const & segParameters, geo::GeometryCore* geom);

            SegmentationAlgPtr_t GetECALSegmentationAlg() const;

            void ConfigureMinervaSegmentationAlg
            (fhicl::ParameterSet const & segParameters, geo::GeometryCore* geom);

            SegmentationAlgPtr_t GetMinervaSegmentationAlg() const;

            void ConfigureMuIDSegmentationAlg
            (fhicl::ParameterSet const & segParameters, geo::GeometryCore* geom);

            SegmentationAlgPtr_t GetMuIDSegmentationAlg() const;

        private:

            /// Implementation of ConfigureChannelMapAlg (pure virtual)
            virtual void doConfigureChannelMapAlg(fhicl::ParameterSet const & sortingParameters,
            gar::geo::GeometryCore* geom) = 0;

            /// Returns the ChannelMapAlg
            virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const    = 0;

            /// Implementation of ConfigureECALSegmentationAlg (pure virtual)
            virtual void doConfigureECALSegmentationAlg(fhicl::ParameterSet const & segParameters,
            gar::geo::GeometryCore* geom) = 0;

            /// Returns the ECAL SegmentationAlg
            virtual SegmentationAlgPtr_t doGetECALSegmentationAlg() const    = 0;

            /// Implementation of ConfigureECALSegmentationAlg (pure virtual)
            virtual void doConfigureMinervaSegmentationAlg(fhicl::ParameterSet const & segParameters,
            gar::geo::GeometryCore* geom) = 0;

            /// Returns the Tracker Sc SegmentationAlg
            virtual SegmentationAlgPtr_t doGetMinervaSegmentationAlg() const    = 0;

            /// Implementation of ConfigureECALSegmentationAlg (pure virtual)
            virtual void doConfigureMuIDSegmentationAlg(fhicl::ParameterSet const & segParameters,
            gar::geo::GeometryCore* geom) = 0;

            /// Returns the MuID SegmentationAlg
            virtual SegmentationAlgPtr_t doGetMuIDSegmentationAlg() const    = 0;

        }; // end ExptGeoHelperInterface class declaration



        //-------------------------------------------------------------------------------------------
        inline void ExptGeoHelperInterface::ConfigureChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
        gar::geo::GeometryCore* geom)
        {
            doConfigureChannelMapAlg(sortingParameters, geom);
        }

        //-------------------------------------------------------------------------------------------
        inline ExptGeoHelperInterface::ChannelMapAlgPtr_t ExptGeoHelperInterface::GetChannelMapAlg() const
        {
            return doGetChannelMapAlg();
        }

        //-------------------------------------------------------------------------------------------
        inline void ExptGeoHelperInterface::ConfigureECALSegmentationAlg(fhicl::ParameterSet const& segParameters,
        gar::geo::GeometryCore* geom)
        {
            doConfigureECALSegmentationAlg(segParameters, geom);
        }

        //-------------------------------------------------------------------------------------------
        inline ExptGeoHelperInterface::SegmentationAlgPtr_t ExptGeoHelperInterface::GetECALSegmentationAlg() const
        {
            return doGetECALSegmentationAlg();
        }

        //-------------------------------------------------------------------------------------------
        inline void ExptGeoHelperInterface::ConfigureMinervaSegmentationAlg(fhicl::ParameterSet const& segParameters,
        gar::geo::GeometryCore* geom)
        {
            doConfigureMinervaSegmentationAlg(segParameters, geom);
        }

        //-------------------------------------------------------------------------------------------
        inline ExptGeoHelperInterface::SegmentationAlgPtr_t ExptGeoHelperInterface::GetMinervaSegmentationAlg() const
        {
            return doGetMinervaSegmentationAlg();
        }

        //-------------------------------------------------------------------------------------------
        inline void ExptGeoHelperInterface::ConfigureMuIDSegmentationAlg(fhicl::ParameterSet const& segParameters,
        gar::geo::GeometryCore* geom)
        {
            doConfigureMuIDSegmentationAlg(segParameters, geom);
        }

        //-------------------------------------------------------------------------------------------
        inline ExptGeoHelperInterface::SegmentationAlgPtr_t ExptGeoHelperInterface::GetMuIDSegmentationAlg() const
        {
            return doGetMuIDSegmentationAlg();
        }
    }
} // gar

DECLARE_ART_SERVICE_INTERFACE(gar::geo::ExptGeoHelperInterface, LEGACY)

#endif // GEO_ExptGeoHelperInterface_h
