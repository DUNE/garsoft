#include "PseudoLayerPlugin.h"

//Pandora
#include "Helpers/XmlHelper.h"
#include "Pandora/AlgorithmHeaders.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace pandora;

namespace gar {
  namespace gar_pandora {

    PseudoLayerPlugin::PseudoLayerPlugin()
      : m_barrelInnerR(0.f)
      , m_endCapInnerZ(0.f)
      , m_rCorrection(0.f)
      , m_zCorrection(0.f)
      , m_barrelEdgeR(0.f)
      , m_endCapEdgeZ(0.f)
    {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode PseudoLayerPlugin::Initialize()
    {
      try {
        this->StoreLayerPositions();
        this->StoreDetectorOuterEdge();
        this->StorePolygonAngles();
        this->StoreOverlapCorrectionDetails();
      }
      catch (StatusCodeException& statusCodeException) {
        std::cout << "PseudoLayerPlugin: Incomplete geometry - consider using a different "
                     "PseudoLayerCalculator."
                  << std::endl;
        return statusCodeException.GetStatusCode();
      }

      return STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    unsigned int PseudoLayerPlugin::GetPseudoLayer(
      const pandora::CartesianVector& positionVector) const
    {
      const float zCoordinate(std::fabs(positionVector.GetZ()));

      if (zCoordinate > m_endCapEdgeZ) {
        MF_LOG_WARNING("PseudoLayerPlugin::GetPseudoLayer")
          << "zCoordinate " << zCoordinate << " > m_endCapEdgeZ " << m_endCapEdgeZ << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
      }

      const float rCoordinate(this->GetMaximumRadius(
        m_eCalBarrelAngleVector, positionVector.GetX(), positionVector.GetY()));

      if (rCoordinate > m_barrelEdgeR) {
        MF_LOG_WARNING("PseudoLayerPlugin::GetPseudoLayer")
          << "rCoordinate " << rCoordinate << " > m_barrelEdgeR " << m_barrelEdgeR << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
      }

      unsigned int pseudoLayer;

      const StatusCode statusCode(this->GetPseudoLayer(rCoordinate,
                                                       zCoordinate,
                                                       m_rCorrection,
                                                       m_zCorrection,
                                                       m_barrelInnerR,
                                                       m_endCapInnerZ,
                                                       pseudoLayer));

      if (STATUS_CODE_SUCCESS != statusCode) throw StatusCodeException(statusCode);

      // Reserve a pseudo layer for track projections, etc.
      return (1 + pseudoLayer);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode PseudoLayerPlugin::GetPseudoLayer(const float rCoordinate,
                                                 const float zCoordinate,
                                                 const float rCorrection,
                                                 const float zCorrection,
                                                 const float barrelInnerR,
                                                 const float endCapInnerZ,
                                                 unsigned int& pseudoLayer) const
    {
      if (zCoordinate < endCapInnerZ) {
        return this->FindMatchingLayer(rCoordinate, m_barrelLayerPositions, pseudoLayer);
      }
      else if (rCoordinate < barrelInnerR) {
        return this->FindMatchingLayer(zCoordinate, m_endCapLayerPositions, pseudoLayer);
      }
      else {
        unsigned int bestBarrelLayer(0);
        const StatusCode barrelStatusCode(this->FindMatchingLayer(
          rCoordinate - rCorrection, m_barrelLayerPositions, bestBarrelLayer));

        unsigned int bestEndCapLayer(0);
        const StatusCode endCapStatusCode(this->FindMatchingLayer(
          zCoordinate - zCorrection, m_endCapLayerPositions, bestEndCapLayer));

        if ((STATUS_CODE_SUCCESS != barrelStatusCode) && (STATUS_CODE_SUCCESS != endCapStatusCode))
          return STATUS_CODE_NOT_FOUND;

        pseudoLayer = std::max(bestBarrelLayer, bestEndCapLayer);
        return STATUS_CODE_SUCCESS;
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode PseudoLayerPlugin::FindMatchingLayer(const float position,
                                                    const LayerPositionList& layerPositionList,
                                                    unsigned int& layer) const
    {
      LayerPositionList::const_iterator upperIter =
        std::upper_bound(layerPositionList.begin(), layerPositionList.end(), position);

      if (layerPositionList.end() == upperIter) { return STATUS_CODE_NOT_FOUND; }

      if (layerPositionList.begin() == upperIter) {
        layer = 0;
        return STATUS_CODE_SUCCESS;
      }

      LayerPositionList::const_iterator lowerIter = upperIter - 1;

      if (std::fabs(position - *lowerIter) < std::fabs(position - *upperIter)) {
        layer = std::distance(layerPositionList.begin(), lowerIter);
      }
      else {
        layer = std::distance(layerPositionList.begin(), upperIter);
      }

      return STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::StoreLayerPositions()
    {
      const GeometryManager* const pGeometryManager(this->GetPandora().GetGeometry());
      this->StoreLayerPositions(pGeometryManager->GetSubDetector(ECAL_BARREL),
                                m_barrelLayerPositions);
      this->StoreLayerPositions(pGeometryManager->GetSubDetector(ECAL_ENDCAP),
                                m_endCapLayerPositions);

      if (m_barrelLayerPositions.empty() || m_endCapLayerPositions.empty()) {
        std::cout << "PseudoLayerPlugin: No layer positions specified." << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
      }

      std::sort(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end());
      std::sort(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end());

      LayerPositionList::const_iterator barrelIter =
        std::unique(m_barrelLayerPositions.begin(), m_barrelLayerPositions.end());
      LayerPositionList::const_iterator endcapIter =
        std::unique(m_endCapLayerPositions.begin(), m_endCapLayerPositions.end());

      if ((m_barrelLayerPositions.end() != barrelIter) ||
          (m_endCapLayerPositions.end() != endcapIter)) {
        std::cout << "PseudoLayerPlugin: Duplicate layer position detected." << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::StoreLayerPositions(const SubDetector& subDetector,
                                                LayerPositionList& layerPositionList)
    {
      if (!subDetector.IsMirroredInZ()) {
        std::cout << "PseudoLayerPlugin: Error, detector must be symmetrical about z=0 plane."
                  << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
      }

      const SubDetector::SubDetectorLayerVector& subDetectorLayerVector(
        subDetector.GetSubDetectorLayerVector());

      for (SubDetector::SubDetectorLayerVector::const_iterator
             iter = subDetectorLayerVector.begin(),
             iterEnd = subDetectorLayerVector.end();
           iter != iterEnd;
           ++iter) {
        layerPositionList.push_back(iter->GetClosestDistanceToIp());
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::StoreDetectorOuterEdge()
    {
      const GeometryManager* const pGeometryManager(this->GetPandora().GetGeometry());

      m_barrelEdgeR = pGeometryManager->GetSubDetector(ECAL_BARREL).GetOuterRCoordinate();

      m_endCapEdgeZ =
        std::fabs(pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetOuterZCoordinate());

      if ((m_barrelLayerPositions.end() != std::upper_bound(m_barrelLayerPositions.begin(),
                                                            m_barrelLayerPositions.end(),
                                                            m_barrelEdgeR)) ||
          (m_endCapLayerPositions.end() != std::upper_bound(m_endCapLayerPositions.begin(),
                                                            m_endCapLayerPositions.end(),
                                                            m_endCapEdgeZ))) {
        std::cout << "PseudoLayerPlugin: Layers specified outside detector edge." << std::endl;
        throw StatusCodeException(STATUS_CODE_FAILURE);
      }

      m_barrelLayerPositions.push_back(m_barrelEdgeR);
      m_endCapLayerPositions.push_back(m_endCapEdgeZ);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::StorePolygonAngles()
    {
      const GeometryManager* const pGeometryManager(this->GetPandora().GetGeometry());

      this->FillAngleVector(pGeometryManager->GetSubDetector(ECAL_BARREL).GetInnerSymmetryOrder(),
                            pGeometryManager->GetSubDetector(ECAL_BARREL).GetInnerPhiCoordinate(),
                            m_eCalBarrelAngleVector);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::StoreOverlapCorrectionDetails()
    {
      const GeometryManager* const pGeometryManager(this->GetPandora().GetGeometry());

      m_barrelInnerR = pGeometryManager->GetSubDetector(ECAL_BARREL).GetInnerRCoordinate();
      m_endCapInnerZ =
        std::fabs(pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetInnerZCoordinate());

      const float barrelOuterZ =
        std::fabs(pGeometryManager->GetSubDetector(ECAL_BARREL).GetOuterZCoordinate());
      const float endCapOuterR =
        pGeometryManager->GetSubDetector(ECAL_ENDCAP).GetOuterRCoordinate();

      const bool IsEnclosingEndCap(endCapOuterR > m_barrelInnerR);
      m_rCorrection =
        ((!IsEnclosingEndCap) ? 0.f : m_barrelInnerR * ((m_endCapInnerZ / barrelOuterZ) - 1.f));
      m_zCorrection =
        ((IsEnclosingEndCap) ? 0.f : m_endCapInnerZ * ((m_barrelInnerR / endCapOuterR) - 1.f));
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    float PseudoLayerPlugin::GetMaximumRadius(const AngleVector& angleVector,
                                              const float x,
                                              const float y) const
    {
      if (angleVector.size() <= 2) return std::sqrt((x * x) + (y * y));

      float maxRadius(0.);
      for (AngleVector::const_iterator iter = angleVector.begin(), iterEnd = angleVector.end();
           iter != iterEnd;
           ++iter) {
        const float radius((x * iter->first) + (y * iter->second));

        if (radius > maxRadius) maxRadius = radius;
      }

      return maxRadius;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PseudoLayerPlugin::FillAngleVector(const unsigned int symmetryOrder,
                                            const float phi0,
                                            AngleVector& angleVector) const
    {
      const float twoPi = static_cast<float>(2. * std::acos(-1.));
      angleVector.clear();

      for (unsigned int iSymmetry = 0; iSymmetry < symmetryOrder; ++iSymmetry) {
        const float phi =
          phi0 + ((twoPi * static_cast<float>(iSymmetry)) / static_cast<float>(symmetryOrder));
        angleVector.push_back(std::make_pair(std::cos(phi), std::sin(phi)));
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode PseudoLayerPlugin::ReadSettings(const pandora::TiXmlHandle /* xmlHandle */)
    {
      return STATUS_CODE_SUCCESS;
    }

  }
}
