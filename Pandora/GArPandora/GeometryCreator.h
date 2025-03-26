#ifndef GEOMETRY_CREATOR_H
#define GEOMETRY_CREATOR_H 1

#include "Api/PandoraApi.h"

#include "Geometry/GeometryGAr.h"

namespace gar {
  namespace gar_pandora {

    class GeometryCreator {
    public:
      class Settings {
      public:
        Settings();
      };

      GeometryCreator(const Settings& settings, const pandora::Pandora* const pPandora);

      ~GeometryCreator();

      pandora::StatusCode CreateGeometry() const;

    private:
      typedef std::map<pandora::SubDetectorType, PandoraApi::Geometry::SubDetector::Parameters>
        SubDetectorTypeMap;
      typedef std::map<std::string, PandoraApi::Geometry::SubDetector::Parameters>
        SubDetectorNameMap;

      void SetMandatorySubDetectorParameters(SubDetectorTypeMap& subDetectorTypeMap) const;

      void SetAdditionalSubDetectorParameters(SubDetectorNameMap& subDetectorNameMap) const;

      void SetDefaultSubDetectorParameters(
        const gar::geo::LayeredCalorimeterData& inputParameters,
        const std::string& subDetectorName,
        const pandora::SubDetectorType subDetectorType,
        PandoraApi::Geometry::SubDetector::Parameters& parameters) const;

      const Settings m_settings;          ///< The geometry creator settings
      const pandora::Pandora& m_pPandora; ///< Address of the pandora object to create the geometry

      const geo::GeometryCore* fGeo; ///< Geometry provider
    };
  }
}

#endif // #ifndef GEOMETRY_CREATOR_H
