////////////////////////////////////////////////////////////////////////////////
/// \file StadardAuxDetGeometryHelper.h
/// \brief Example of the AuxDetGeoHelper service for Auxiliary detector geometries.
///        Do not attempt to use this in practice, need to define experiment
///        specific helpers
///
/// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef AUXDETExptGeoHelperExample_h
#define AUXDETExptGeoHelperExample_h

#include "Geometry/AuxDetExptGeoHelperInterface.h"

#include <memory> // std::shared_ptr<>

// Declaration
//
namespace gar {
  namespace geo {
    class AuxDetGeometryHelperExample : public geo::AuxDetExptGeoHelperInterface {
    public:
      AuxDetGeometryHelperExample(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

      /*
       Public interface for AuxDetExptGeoHelperInterface (for reference purposes)

       Configure, initialize and return the channel map:

       void ConfigureChannelMapAlg
       (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom);

       Returns null pointer if the initialization failed:

       ChannelMapAlgPtr_t GetChannelMapAlg() const;
       */

    private:
      virtual void doConfigureAuxDetChannelMapAlg(fhicl::ParameterSet const& sortingParameters,
                                                  geo::AuxDetGeometryCore* geom) override;
      virtual AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;

      fhicl::ParameterSet fPset; ///< copy of configuration parameter set
      std::shared_ptr<geo::seg::AuxDetChannelMapAlg> fChannelMap; ///< channel map
    };

  }
}

DECLARE_ART_SERVICE_INTERFACE_IMPL(gar::geo::AuxDetGeometryHelperExample,
                                   gar::geo::AuxDetExptGeoHelperInterface,
                                   LEGACY)

#endif // StandardAUXDETExptGeoHelperExample_h
