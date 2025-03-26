////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTER_H
#define GEO_AUXDETGEOOBJECTSORTER_H

#include <vector>

#include "fhiclcpp/ParameterSet.h"

namespace gar {
  namespace geo {

    class AuxDetGeo;
    class AuxDetSensitiveGeo;

    class AuxDetGeoObjectSorter {

    public:
      explicit AuxDetGeoObjectSorter(fhicl::ParameterSet const& pset);
      ~AuxDetGeoObjectSorter();

      virtual void SortAuxDets(std::vector<geo::AuxDetGeo*>& adgeo) const = 0;
      virtual void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*>& adsgeo) const = 0;

    private:
    };

  }
} // gar
#endif // GEO_GEOOBJECTSORTER_H
