////////////////////////////////////////////////////////////////////////////////
/// \file AuxDetGeometryHelperExample_service.cc
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geometry/AuxDetGeometryCore.h"
#include "Geometry/AuxDetGeometryHelperExample.h"

#include <memory> // std::make_shared()

namespace gar {
  namespace geo {

    //------------------------------------------------------------------------
    AuxDetGeometryHelperExample::AuxDetGeometryHelperExample(fhicl::ParameterSet const& pset,
                                                             art::ActivityRegistry&)
      : fPset(pset), fChannelMap()
    {}

    //------------------------------------------------------------------------
    void AuxDetGeometryHelperExample::doConfigureAuxDetChannelMapAlg(
      fhicl::ParameterSet const&,    //sortingParameters,
      gar::geo::AuxDetGeometryCore*) //  geom)
    {

      // This is where your experiment specific helper has to instantiate its
      // sorting algorithm
      // fChannelMap = std::make_shared<gar::geo::AuxDetChannelMapStandardAlg>(sortingParameters);
      // if(fChannelMap) geom->ApplyChannelMap(fChannelMap);

      return;
    }

    //------------------------------------------------------------------------
    AuxDetGeometryHelperExample::AuxDetChannelMapAlgPtr_t
    AuxDetGeometryHelperExample::doGetAuxDetChannelMapAlg() const
    {
      return fChannelMap;
    }

  }
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::geo::AuxDetGeometryHelperExample,
                                  gar::geo::AuxDetExptGeoHelperInterface)
