#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Geometry/BitFieldCoder.h"
#include "Geometry/GeometryGAr.h"
#include "ReconstructionDataProducts/CaloHit.h"

#include "RotationTransformation.h"

#include "Api/PandoraApi.h"

namespace gar {
  namespace gar_pandora {

    typedef std::vector<art::Ptr<gar::rec::CaloHit>> CalorimeterHitVector;
    typedef std::vector<gar::rec::CaloHit> RawCalorimeterHitVector;

    class CaloHitCreator {
    public:
      typedef std::vector<float> FloatVector;

      class Settings {
      public:
        Settings();

        std::string m_CaloHitCollection; ///< The calorimeter hit collection

        std::string m_CaloHitInstanceName; ///< The calorimeter hit instance name

        float m_eCalToMip;        ///< The calibration from deposited ECal energy to mip
        float m_eCalMipThreshold; ///< Threshold for creating calo hits in the ECal, units mip
        float m_eCalToEMGeV;      ///< The calibration from deposited ECal energy to EM energy
        float
          m_eCalToHadGeVBarrel; ///< The calibration from deposited ECal barrel energy to hadronic energy
        float
          m_eCalToHadGeVEndCap; ///< The calibration from deposited ECal endcap energy to hadronic energy
        float
          m_maxECalHitHadronicEnergy; ///< The maximum hadronic energy allowed for a single hcal hit
        int
          m_nOuterSamplingLayers; ///< Number of layers from edge for hit to be flagged as an outer layer hit
        float
          m_layersFromEdgeMaxRearDistance; ///< Maximum number of layers from candidate outer layer hit to rear of detector
        float m_eCalBarrelOuterR;                    ///< ECal barrel outer r coordinate
        float m_eCalBarrelOuterZ;                    ///< ECal barrel outer z coordinate
        float m_eCalBarrelInnerPhi0;                 ///< ECal barrel inner phi0 coordinate
        float m_eCalBarrelOuterPhi0;                 ///< ECal barrel outer phi0 coordinate
        unsigned int m_eCalBarrelInnerSymmetry;      ///< ECal barrel inner symmetry order
        unsigned int m_eCalBarrelOuterSymmetry;      ///< ECal barrel outer symmetry order
        FloatVector m_eCalBarrelNormalVector;        ///< ECal barrel normal vector
        float m_eCalEndCapOuterR;                    ///< ECal endcap outer r coordinate
        float m_eCalEndCapOuterZ;                    ///< ECal endcap outer z coordinate
        unsigned int m_eCalEndCapInnerSymmetryOrder; ///< ECal endcap inner symmetry
        float m_eCalEndCapInnerPhiCoordinate;        ///< ECal endcap inner phi
      };

      CaloHitCreator(const Settings& settings,
                     const pandora::Pandora* const pPandora,
                     const RotationTransformation* const pRotation);

      ~CaloHitCreator();

      pandora::StatusCode CreateCaloHits(const art::Event& pEvent);
      const CalorimeterHitVector& GetCalorimeterHitVector() const;

      void Reset();

    private:
      pandora::StatusCode CreateECalCaloHits() const;

      void GetCommonCaloHitProperties(const gar::rec::CaloHit* const pCaloHit,
                                      PandoraApi::CaloHit::Parameters& caloHitParameters) const;
      void GetEndCapCaloHitProperties(
        const gar::rec::CaloHit* const pCaloHit,
        const std::vector<gar::geo::LayeredCalorimeterStruct::Layer>& layers,
        PandoraApi::CaloHit::Parameters& caloHitParameters,
        float& absorberCorrection) const;
      void GetBarrelCaloHitProperties(
        const gar::rec::CaloHit* const pCaloHit,
        const std::vector<gar::geo::LayeredCalorimeterStruct::Layer>& layers,
        unsigned int barrelSymmetryOrder,
        PandoraApi::CaloHit::Parameters& caloHitParameters,
        FloatVector const& normalVector,
        float& absorberCorrection) const;
      int GetNLayersFromEdge(const gar::rec::CaloHit* const pCaloHit) const;
      float GetMaximumRadius(const gar::rec::CaloHit* const pCaloHit,
                             const unsigned int symmetryOrder,
                             const float phi0) const;

      pandora::StatusCode CollectECALCaloHits(const art::Event& pEvent,
                                              const std::string& label,
                                              const std::string& instanceName,
                                              CalorimeterHitVector& ecalCaloHitVector);

      const Settings m_settings;         ///< The calo hit creator settings
      const pandora::Pandora& m_pandora; ///< Reference to the pandora object to create calo hits
      const geo::GeometryCore* fGeo;     //Geometry Manager
      gar::geo::BitFieldCoder const* m_fieldDecoder;
      float m_origin[3] = {0, 0, 0};

      const RotationTransformation& m_rotation;

      float m_eCalBarrelLayerThickness; ///< ECal barrel layer thickness
      float m_eCalEndCapLayerThickness; ///< ECal endcap layer thickness

      CalorimeterHitVector artCalorimeterHitVector;
    };

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline const CalorimeterHitVector& CaloHitCreator::GetCalorimeterHitVector() const
    {
      return artCalorimeterHitVector;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    inline void CaloHitCreator::Reset()
    {
      artCalorimeterHitVector.clear();
    }

  }
}

#endif // #ifndef CALO_HIT_CREATOR_H
