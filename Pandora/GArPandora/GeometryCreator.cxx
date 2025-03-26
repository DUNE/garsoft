#include "GeometryCreator.h"

#include "CoreUtils/ServiceUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <exception>

namespace gar {
  namespace gar_pandora {

    GeometryCreator::GeometryCreator(const Settings& settings,
                                     const pandora::Pandora* const pPandora)
      : m_settings(settings), m_pPandora(*pPandora)
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    GeometryCreator::~GeometryCreator() {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode GeometryCreator::CreateGeometry() const
    {
      try {
        SubDetectorTypeMap subDetectorTypeMap;
        this->SetMandatorySubDetectorParameters(subDetectorTypeMap);

        SubDetectorNameMap subDetectorNameMap;
        this->SetAdditionalSubDetectorParameters(subDetectorNameMap);

        std::string detectorName = "NDGAr";

        MF_LOG_INFO("GeometryCreator::CreateGeometry()")
          << "Creating geometry for detector " << detectorName << std::endl;

        for (SubDetectorTypeMap::const_iterator iter = subDetectorTypeMap.begin(),
                                                iterEnd = subDetectorTypeMap.end();
             iter != iterEnd;
             ++iter)
          PANDORA_RETURN_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));

        for (SubDetectorNameMap::const_iterator iter = subDetectorNameMap.begin(),
                                                iterEnd = subDetectorNameMap.end();
             iter != iterEnd;
             ++iter) {
          MF_LOG_INFO("GeometryCreator::CreateGeometry()")
            << "Creating geometry for additional subdetector " << iter->first << std::endl;
          PANDORA_RETURN_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            PandoraApi::Geometry::SubDetector::Create(m_pPandora, iter->second));
        }
      }
      catch (std::exception& exception) {
        MF_LOG_ERROR("GeometryCreator::CreateGeometry()")
          << "Failure in pandora geometry creator, exception: " << exception.what() << std::endl;
        throw exception;
      }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void GeometryCreator::SetMandatorySubDetectorParameters(
      SubDetectorTypeMap& subDetectorTypeMap) const
    {
      PandoraApi::Geometry::SubDetector::Parameters eCalBarrelParameters, eCalEndCapParameters,
        muonBarrelParameters, muonEndCapParameters;

      this->SetDefaultSubDetectorParameters(
        *const_cast<gar::geo::LayeredCalorimeterData*>(
          fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::BarrelLayout]
            .get()),
        "ECalBarrel",
        pandora::ECAL_BARREL,
        eCalBarrelParameters);
      this->SetDefaultSubDetectorParameters(
        *const_cast<gar::geo::LayeredCalorimeterData*>(
          fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::EndcapLayout]
            .get()),
        "ECalEndcap",
        pandora::ECAL_ENDCAP,
        eCalEndCapParameters);

      subDetectorTypeMap[pandora::ECAL_BARREL] = eCalBarrelParameters;
      subDetectorTypeMap[pandora::ECAL_ENDCAP] = eCalEndCapParameters;

      PandoraApi::Geometry::SubDetector::Parameters trackerParameters;

      trackerParameters.m_subDetectorName = "Tracker";
      trackerParameters.m_subDetectorType = pandora::INNER_TRACKER;
      trackerParameters.m_innerRCoordinate = 0.f; // inner R of TPC
      trackerParameters.m_innerZCoordinate = 0.f;
      trackerParameters.m_innerPhiCoordinate = 0.f;
      trackerParameters.m_innerSymmetryOrder = 0;
      trackerParameters.m_outerRCoordinate = fGeo->TPCRadius() * CLHEP::cm;      // outer R of TPC
      trackerParameters.m_outerZCoordinate = fGeo->TPCLength() / 2. * CLHEP::cm; // outer X of TPC
      trackerParameters.m_outerPhiCoordinate = 0.f;
      trackerParameters.m_outerSymmetryOrder = 0;
      trackerParameters.m_isMirroredInZ = true;
      trackerParameters.m_nLayers = 0;
      subDetectorTypeMap[pandora::INNER_TRACKER] = trackerParameters;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void GeometryCreator::SetAdditionalSubDetectorParameters(
      SubDetectorNameMap& /* subDetectorNameMap */) const
    {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    void GeometryCreator::SetDefaultSubDetectorParameters(
      const gar::geo::LayeredCalorimeterData& inputParameters,
      const std::string& subDetectorName,
      const pandora::SubDetectorType subDetectorType,
      PandoraApi::Geometry::SubDetector::Parameters& parameters) const
    {
      const std::vector<gar::geo::LayeredCalorimeterStruct::Layer>& layers = inputParameters.layers;

      parameters.m_subDetectorName = subDetectorName;
      parameters.m_subDetectorType = subDetectorType;
      parameters.m_innerRCoordinate = inputParameters.extent[0] * CLHEP::cm;
      parameters.m_innerZCoordinate = inputParameters.extent[2] * CLHEP::cm;
      parameters.m_innerSymmetryOrder = inputParameters.inner_symmetry;
      parameters.m_outerRCoordinate = inputParameters.extent[1] * CLHEP::cm;
      parameters.m_outerZCoordinate = inputParameters.extent[3] * CLHEP::cm;
      parameters.m_innerPhiCoordinate = inputParameters.inner_phi0;
      parameters.m_outerPhiCoordinate = inputParameters.outer_phi0;
      parameters.m_outerSymmetryOrder = inputParameters.outer_symmetry;
      parameters.m_isMirroredInZ = true;
      parameters.m_nLayers = layers.size();

      MF_LOG_DEBUG("GeometryCreator::SetDefaultSubDetectorParameters")
        << " parameters.m_subDetectorName = " << parameters.m_subDetectorName.Get()
        << " parameters.m_subDetectorType = " << parameters.m_subDetectorType.Get()
        << " parameters.m_innerRCoordinate = " << parameters.m_innerRCoordinate.Get()
        << " parameters.m_innerZCoordinate = " << parameters.m_innerZCoordinate.Get()
        << " parameters.m_innerSymmetryOrder = " << parameters.m_innerSymmetryOrder.Get()
        << " parameters.m_outerRCoordinate = " << parameters.m_outerRCoordinate.Get()
        << " parameters.m_outerZCoordinate = " << parameters.m_outerZCoordinate.Get()
        << " parameters.m_innerPhiCoordinate = " << parameters.m_innerPhiCoordinate.Get()
        << " parameters.m_outerPhiCoordinate = " << parameters.m_outerPhiCoordinate.Get()
        << " parameters.m_nLayers = " << parameters.m_nLayers.Get();

      for (size_t i = 0; i < layers.size(); i++) {
        const gar::geo::LayeredCalorimeterStruct::Layer& theLayer = layers.at(i);
        PandoraApi::Geometry::LayerParameters layerParameters;

        double totalNumberOfRadLengths = theLayer.inner_nRadiationLengths;
        double totalNumberOfIntLengths = theLayer.inner_nInteractionLengths;

        if (i > 0) {
          //Add the numbers from previous layer's outer side
          totalNumberOfRadLengths += layers.at(i - 1).outer_nRadiationLengths;
          totalNumberOfIntLengths += layers.at(i - 1).outer_nInteractionLengths;
        }

        layerParameters.m_closestDistanceToIp = (theLayer.distance + theLayer.inner_thickness) *
                                                CLHEP::cm; //Distance to center of sensitive element
        layerParameters.m_nRadiationLengths = totalNumberOfRadLengths;
        layerParameters.m_nInteractionLengths = totalNumberOfIntLengths;

        MF_LOG_DEBUG("GeometryCreator::SetDefaultSubDetectorParameters")
          << " layer " << i << " layerParameters.m_closestDistanceToIp = "
          << layerParameters.m_closestDistanceToIp.Get()
          << " layerParameters.m_nRadiationLengths = " << layerParameters.m_nRadiationLengths.Get()
          << " layerParameters.m_nInteractionLengths = "
          << layerParameters.m_nInteractionLengths.Get();

        parameters.m_layerParametersVector.push_back(layerParameters);
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    GeometryCreator::Settings::Settings() {}
  }
}
