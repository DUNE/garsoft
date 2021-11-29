#include "CaloHitCreator.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <exception>

#include "CoreUtils/ServiceUtil.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace gar_pandora {

    CaloHitCreator::CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora, const RotationTransformation *const pRotation)
      : m_settings(settings),
        m_pandora(*pPandora),
        m_rotation(*pRotation),
        m_eCalBarrelLayerThickness(0.f),
        m_eCalEndCapLayerThickness(0.f),
        artCalorimeterHitVector(0)
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();
      std::string fEncoding = fGeo->GetECALCellIDEncoding();
      m_fieldDecoder = new gar::geo::BitFieldCoder( fEncoding );
      m_origin[0] = fGeo->GetOriginX();
      m_origin[1] = fGeo->GetOriginY();
      m_origin[2] = fGeo->GetOriginZ();

      const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &barrelLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::BarrelLayout].get())->layers;

      const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &endcapLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::EndcapLayout].get())->layers;

      ///Take thicknesses from last layer (was like that before with gear)
      m_eCalEndCapLayerThickness = (endcapLayers.back().inner_thickness + endcapLayers.back().outer_thickness) * CLHEP::cm;
      m_eCalBarrelLayerThickness = (barrelLayers.back().inner_thickness + barrelLayers.back().outer_thickness) * CLHEP::cm;

      if ( (m_eCalEndCapLayerThickness < std::numeric_limits<float>::epsilon()) || (m_eCalBarrelLayerThickness < std::numeric_limits<float>::epsilon()) )
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    CaloHitCreator::~CaloHitCreator()
    {
      if ( m_fieldDecoder ) delete m_fieldDecoder;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode CaloHitCreator::CreateCaloHits(const art::Event &pEvent)
    {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectECALCaloHits(pEvent, m_settings.m_CaloHitCollection, m_settings.m_CaloHitInstanceName, artCalorimeterHitVector));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateECalCaloHits());

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode CaloHitCreator::CreateECalCaloHits() const
    {
      const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &barrelLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::BarrelLayout].get())->layers;

      const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &endcapLayers = (fGeo->GetECALLayeredCalorimeterData()[gar::geo::LayeredCalorimeterData::EndcapLayout].get())->layers;

      for (CalorimeterHitVector::const_iterator iter = artCalorimeterHitVector.begin(), iterEnd = artCalorimeterHitVector.end(); iter != iterEnd; ++iter)
        {
          try
            {
              art::Ptr<gar::rec::CaloHit> artPtrCaloHit = *iter;
              const gar::rec::CaloHit *pCaloHit = artPtrCaloHit.get();

              if (nullptr == pCaloHit)
                throw;

              PandoraApi::CaloHit::Parameters caloHitParameters;
              caloHitParameters.m_hitType = pandora::ECAL;
              caloHitParameters.m_isDigital = false;
              caloHitParameters.m_layer = m_fieldDecoder->get(pCaloHit->CellID(), "layer");
              caloHitParameters.m_isInOuterSamplingLayer = (this->GetNLayersFromEdge(pCaloHit) <= m_settings.m_nOuterSamplingLayers);
              this->GetCommonCaloHitProperties(pCaloHit, caloHitParameters);

              float absorberCorrection(1.);

              //Need to use X instead of Z
              if (std::fabs( (pCaloHit->Position()[0] - m_origin[0]) * CLHEP::cm) < m_settings.m_eCalBarrelOuterZ)
                {
                  this->GetBarrelCaloHitProperties(pCaloHit, barrelLayers, m_settings.m_eCalBarrelInnerSymmetry, caloHitParameters, m_settings.m_eCalBarrelNormalVector, absorberCorrection);
                  caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVBarrel * pCaloHit->Energy();
                }
              else
                {
                  this->GetEndCapCaloHitProperties(pCaloHit, endcapLayers, caloHitParameters, absorberCorrection);
                  caloHitParameters.m_hadronicEnergy = m_settings.m_eCalToHadGeVEndCap * pCaloHit->Energy();
                }

              caloHitParameters.m_mipEquivalentEnergy = pCaloHit->Energy() * m_settings.m_eCalToMip * absorberCorrection;

              if (caloHitParameters.m_mipEquivalentEnergy.Get() < m_settings.m_eCalMipThreshold)
                continue;

              caloHitParameters.m_electromagneticEnergy = m_settings.m_eCalToEMGeV * pCaloHit->Energy();

              //Check if the hit is inside the calo... could be SSA reco problems....
              const float rCoordinate( std::sqrt( (pCaloHit->Position()[1] - m_origin[1]) * (pCaloHit->Position()[1] - m_origin[1]) + (pCaloHit->Position()[2] - m_origin[2]) * (pCaloHit->Position()[2] - m_origin[2]) ) * CLHEP::cm );
              if( rCoordinate > m_settings.m_eCalBarrelOuterR )
                {
                  MF_LOG_WARNING("CaloHitCreator::CreateECalCaloHits")
                    << " Position x: " << pCaloHit->Position()[0] * CLHEP::cm
                    << " y: " << pCaloHit->Position()[1] * CLHEP::cm
                    << " z: " << pCaloHit->Position()[2] * CLHEP::cm
                    << " rCoordinate " << rCoordinate << " > " << " m_settings.m_eCalBarrelOuterR " << m_settings.m_eCalBarrelOuterR
                    << " caloHitParameters.m_hitType = " << caloHitParameters.m_hitType.Get()
                    << " caloHitParameters.m_isDigital = " << caloHitParameters.m_isDigital.Get()
                    << " caloHitParameters.m_layer = " << caloHitParameters.m_layer.Get()
                    << " caloHitParameters.m_isInOuterSamplingLayer = " << caloHitParameters.m_isInOuterSamplingLayer.Get()
                    << " caloHitParameters.m_cellGeometry = " << caloHitParameters.m_cellGeometry.Get()
                    << " caloHitParameters.m_expectedDirection = " << caloHitParameters.m_expectedDirection.Get()
                    << " caloHitParameters.m_inputEnergy = " << caloHitParameters.m_inputEnergy.Get()
                    << " caloHitParameters.m_time = " << caloHitParameters.m_time.Get()
                    << " caloHitParameters.m_cellSize0 = " << caloHitParameters.m_cellSize0.Get()
                    << " caloHitParameters.m_cellSize1 = " << caloHitParameters.m_cellSize1.Get()
                    << " caloHitParameters.m_cellThickness = " << caloHitParameters.m_cellThickness.Get()
                    << " caloHitParameters.m_nCellRadiationLengths = " << caloHitParameters.m_nCellRadiationLengths.Get()
                    << " caloHitParameters.m_nCellInteractionLengths = " << caloHitParameters.m_nCellInteractionLengths.Get()
                    << " caloHitParameters.m_hadronicEnergy = " << caloHitParameters.m_hadronicEnergy.Get()
                    << " caloHitParameters.m_mipEquivalentEnergy = " << caloHitParameters.m_mipEquivalentEnergy.Get()
                    << " caloHitParameters.m_electromagneticEnergy = " << caloHitParameters.m_electromagneticEnergy.Get();
                  continue;
                }

              PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(m_pandora, caloHitParameters));
            }
          catch (pandora::StatusCodeException &statusCodeException)
            {
              MF_LOG_ERROR("CaloHitCreator::CreateECalCaloHits()")
                << "Failed to extract ecal calo hit: " << statusCodeException.ToString();
            }
          catch (std::exception &exception)
            {
              MF_LOG_WARNING("CaloHitCreator::CreateECalCaloHits()")
                << "Failed to extract ecal calo hit: " << exception.what();
            }
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void CaloHitCreator::GetCommonCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const
    {
      const float *pCaloHitPosition(pCaloHit->Position());
      //Inverse X and Z for pandora to cope with the change in beam axis
      const pandora::CartesianVector positionVector( (pCaloHitPosition[0] - m_origin[0]) * CLHEP::cm, (pCaloHitPosition[1] - m_origin[1]) * CLHEP::cm, (pCaloHitPosition[2] - m_origin[2]) * CLHEP::cm);
      const pandora::CartesianVector newPositionVector = m_rotation.MakeRotation(positionVector);

      caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
      caloHitParameters.m_positionVector = newPositionVector;
      caloHitParameters.m_expectedDirection = newPositionVector.GetUnitVector();
      caloHitParameters.m_pParentAddress = (void*)pCaloHit;
      caloHitParameters.m_inputEnergy = pCaloHit->Energy();
      caloHitParameters.m_time = pCaloHit->Time().first;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void CaloHitCreator::GetEndCapCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &layers, PandoraApi::CaloHit::Parameters &caloHitParameters, float &absorberCorrection) const
    {
      caloHitParameters.m_hitRegion = pandora::ENDCAP;

      const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
      caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0 * CLHEP::cm;
      caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1 * CLHEP::cm;

      double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
      double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
      double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;
      double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;

      if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
      }

      caloHitParameters.m_cellThickness = thickness;
      caloHitParameters.m_nCellRadiationLengths = nRadLengths;
      caloHitParameters.m_nCellInteractionLengths = nIntLengths;

      if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
        {
          MF_LOG_ERROR("CaloHitCreator::GetEndcapCaloHitProperties")
            << "Calo hit has 0 radiation length or interaction length: not creating a Pandora calo hit.";
          throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }

      absorberCorrection = 1.;
      for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
        {
          float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 ) * CLHEP::cm);

          if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0) * CLHEP::cm;

          if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

          if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

          break;
        }

      caloHitParameters.m_cellNormalVector = (pCaloHit->Position()[0] * CLHEP::cm > 0) ? pandora::CartesianVector(1, 0, 0) : pandora::CartesianVector(-1, 0, 0);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void CaloHitCreator::GetBarrelCaloHitProperties( const gar::rec::CaloHit *const pCaloHit,
                                                     const std::vector<gar::geo::LayeredCalorimeterStruct::Layer> &layers, unsigned int barrelSymmetryOrder, PandoraApi::CaloHit::Parameters &caloHitParameters, FloatVector const& /* normalVector */, float &absorberCorrection ) const
    {
      caloHitParameters.m_hitRegion = pandora::BARREL;

      const int physicalLayer(std::min(static_cast<int>(caloHitParameters.m_layer.Get()), static_cast<int>(layers.size()-1)));
      caloHitParameters.m_cellSize0 = layers[physicalLayer].cellSize0 * CLHEP::cm;
      caloHitParameters.m_cellSize1 = layers[physicalLayer].cellSize1 * CLHEP::cm;

      double thickness = (layers[physicalLayer].inner_thickness+layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
      double nRadLengths = layers[physicalLayer].inner_nRadiationLengths;
      double nIntLengths = layers[physicalLayer].inner_nInteractionLengths;

      double layerAbsorberThickness = (layers[physicalLayer].inner_thickness-layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
      if(physicalLayer>0){
        thickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
        nRadLengths += layers[physicalLayer-1].outer_nRadiationLengths;
        nIntLengths += layers[physicalLayer-1].outer_nInteractionLengths;
        layerAbsorberThickness += (layers[physicalLayer-1].outer_thickness -layers[physicalLayer].sensitive_thickness/2.0) * CLHEP::cm;
      }

      caloHitParameters.m_cellThickness = thickness;
      caloHitParameters.m_nCellRadiationLengths = nRadLengths;
      caloHitParameters.m_nCellInteractionLengths = nIntLengths;

      if (caloHitParameters.m_nCellRadiationLengths.Get() < std::numeric_limits<float>::epsilon() || caloHitParameters.m_nCellInteractionLengths.Get() < std::numeric_limits<float>::epsilon())
        {
          MF_LOG_ERROR("CaloHitCreator::GetBarrelCaloHitProperties")
            << "Calo hit has 0 radiation length or interaction length: not creating a Pandora calo hit.";
          throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }

      absorberCorrection = 1.;
      for (unsigned int i = 0, iMax = layers.size(); i < iMax; ++i)
        {
          float absorberThickness((layers[i].inner_thickness - layers[i].sensitive_thickness/2.0 ) * CLHEP::cm);

          if (i>0)
            absorberThickness += (layers[i-1].outer_thickness - layers[i-1].sensitive_thickness/2.0) * CLHEP::cm;

          if (absorberThickness < std::numeric_limits<float>::epsilon())
            continue;

          if (layerAbsorberThickness > std::numeric_limits<float>::epsilon())
            absorberCorrection = absorberThickness / layerAbsorberThickness;

          break;
        }

      if (barrelSymmetryOrder > 2)
        {
          const double phi = atan2( pCaloHit->Position()[2], pCaloHit->Position()[1] );

          MF_LOG_DEBUG( "CaloHitCreator::GetBarrelCaloHitProperties" )
                                                            << "This hit does not have any cellIDs set, will use phi-direction for normal vector "
                                                            << " phi:" << std::setw(15) << phi*180/M_PI;

          caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0, std::cos(phi), std::sin(phi));
        }
      else
        {
          const float *pCaloHitPosition( pCaloHit->Position() );
          const float phi = std::atan2( pCaloHitPosition[2], pCaloHitPosition[1] );
          caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0, std::cos(phi), std::sin(phi));
        }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    int CaloHitCreator::GetNLayersFromEdge(const gar::rec::CaloHit *const pCaloHit) const
    {
      // Calo hit coordinate calculations
      const float barrelMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_eCalBarrelOuterSymmetry, m_settings.m_eCalBarrelOuterPhi0));
      const float endCapMaximumRadius(this->GetMaximumRadius(pCaloHit, m_settings.m_eCalEndCapInnerSymmetryOrder, m_settings.m_eCalEndCapInnerPhiCoordinate));
      const float caloHitAbsZ(std::fabs(pCaloHit->Position()[0]) * CLHEP::cm);

      // Distance from radial outer
      float radialDistanceToEdge(std::numeric_limits<float>::max());

      if (caloHitAbsZ < m_settings.m_eCalBarrelOuterZ)
        {
          radialDistanceToEdge = (m_settings.m_eCalBarrelOuterR - barrelMaximumRadius) / m_eCalBarrelLayerThickness;
        }
      else
        {
          radialDistanceToEdge = (m_settings.m_eCalEndCapOuterR - endCapMaximumRadius) / m_eCalEndCapLayerThickness;
        }

      // Distance from rear of endcap outer
      float rearDistanceToEdge(std::numeric_limits<float>::max());

      if (caloHitAbsZ >= m_settings.m_eCalBarrelOuterZ)
        {
          rearDistanceToEdge = (m_settings.m_eCalEndCapOuterZ - caloHitAbsZ) / m_eCalEndCapLayerThickness;
        }
      else
        {
          const float rearDistance((m_settings.m_eCalBarrelOuterZ - caloHitAbsZ) / m_eCalBarrelLayerThickness);

          if (rearDistance < m_settings.m_layersFromEdgeMaxRearDistance)
            {
              const float overlapDistance((m_settings.m_eCalEndCapOuterR - endCapMaximumRadius) / m_eCalEndCapLayerThickness);
              rearDistanceToEdge = std::max(rearDistance, overlapDistance);
            }
        }

      return static_cast<int>(std::min(radialDistanceToEdge, rearDistanceToEdge));
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    float CaloHitCreator::GetMaximumRadius(const gar::rec::CaloHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const
    {
      const float *pCaloHitPosition(pCaloHit->Position());

      if (symmetryOrder <= 2)
        return std::sqrt((pCaloHitPosition[1] * pCaloHitPosition[1]) + (pCaloHitPosition[2] * pCaloHitPosition[2]) * CLHEP::cm);

      float maximumRadius(0.f);
      const float twoPi(2.f * M_PI);

      for (unsigned int i = 0; i < symmetryOrder; ++i)
        {
          const float phi = phi0 + i * twoPi / static_cast<float>(symmetryOrder);
          float radius = pCaloHitPosition[2] * std::cos(phi) + pCaloHitPosition[1] * std::sin(phi);

          if (radius > maximumRadius)
            maximumRadius = radius;
        }

      return maximumRadius * CLHEP::cm;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode CaloHitCreator::CollectECALCaloHits(const art::Event &pEvent, const std::string &label, const std::string &instanceName, CalorimeterHitVector &ecalCaloHitVector)
    {
      art::InputTag itag(label, instanceName);
      auto theHits = pEvent.getHandle< RawCalorimeterHitVector >(itag);
      if (!theHits)
        {
          MF_LOG_WARNING("CaloHitCreator::CreateECalCaloHits") << "  Failed to find ecal hits for label " << label << std::endl;
          return pandora::STATUS_CODE_NOT_FOUND;
        }

      MF_LOG_DEBUG("CaloHitCreator::CreateECalCaloHits") << "  Found: " << theHits->size() << " Hits " << std::endl;

      for (unsigned int i = 0; i < theHits->size(); ++i)
        {
          const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
          ecalCaloHitVector.push_back(hit);
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    CaloHitCreator::Settings::Settings()
      : m_CaloHitCollection( "" ),
        m_CaloHitInstanceName( "" ),
        m_eCalToMip(1.f),
        m_eCalMipThreshold(0.f),
        m_eCalToEMGeV(1.f),
        m_eCalToHadGeVBarrel(1.f),
        m_eCalToHadGeVEndCap(1.f),
        m_maxECalHitHadronicEnergy(10000.f),
        m_nOuterSamplingLayers(3),
        m_layersFromEdgeMaxRearDistance(250.f),
        m_eCalBarrelOuterR(0.f),
        m_eCalBarrelOuterZ(0.f),
        m_eCalBarrelInnerPhi0(0.f),
        m_eCalBarrelOuterPhi0(0.f),
        m_eCalBarrelInnerSymmetry(0.f),
        m_eCalBarrelOuterSymmetry(8),
        m_eCalBarrelNormalVector({0.0, 0.0, 1.0}),
        m_eCalEndCapOuterR(0.f),
        m_eCalEndCapOuterZ(0.f),
        m_eCalEndCapInnerSymmetryOrder(0),
        m_eCalEndCapInnerPhiCoordinate(0.f)
    {
    }
  }
}
