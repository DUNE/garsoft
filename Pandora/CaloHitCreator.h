#ifndef CALO_HIT_CREATOR_H
#define CALO_HIT_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Geometry/Geometry.h"
#include "Geometry/BitFieldCoder.h"
#include "ReconstructionDataProducts/CaloHit.h"

#include "Api/PandoraApi.h"

namespace gar {
    namespace gar_pandora {

        typedef std::vector<const gar::rec::CaloHit *> CalorimeterHitVector;

        class CaloHitCreator
        {
        public:
            typedef std::vector<float> FloatVector;

            class Settings
            {
            public:
                Settings();

                std::string     m_CaloHitCollection;               ///< The calorimeter hit collection
                
                float           m_eCalToMip;                            ///< The calibration from deposited ECal energy to mip
                float           m_eCalMipThreshold;                     ///< Threshold for creating calo hits in the ECal, units mip

                float           m_eCalToEMGeV;                          ///< The calibration from deposited ECal energy to EM energy
                float           m_eCalToHadGeVBarrel;                   ///< The calibration from deposited ECal barrel energy to hadronic energy
                float           m_eCalToHadGeVEndCap;                   ///< The calibration from deposited ECal endcap energy to hadronic energy

                float           m_maxECalHitHadronicEnergy;             ///< The maximum hadronic energy allowed for a single hcal hit
                int             m_nOuterSamplingLayers;                 ///< Number of layers from edge for hit to be flagged as an outer layer hit
                float           m_layersFromEdgeMaxRearDistance;        ///< Maximum number of layers from candidate outer layer hit to rear of detector

                float           m_eCalBarrelInnerPhi0;              ///< ECal barrel inner phi0 coordinate
                unsigned int    m_eCalBarrelInnerSymmetry;          ///< ECal barrel inner symmetry order

                float           m_eCalEndCapOuterR;                 ///< ECal endcap outer r coordinate
                float           m_eCalEndCapInnerX;                 ///< ECal endcap inner x coordinate
                float           m_eCalEndCapOuterX;                 ///< ECal endcap outer x coordinate
                float           m_eCalBarrelOuterR;                 ///< ECal barrel outer r coordinate

                FloatVector m_eCalBarrelNormalVector;
            };

            CaloHitCreator(const Settings &settings, const pandora::Pandora *const pPandora);

            ~CaloHitCreator();

            pandora::StatusCode CreateCaloHits(const art::Event *const pEvent);

            const CalorimeterHitVector &GetCalorimeterHitVector() const;

            void Reset();

        private:

            pandora::StatusCode CreateECalCaloHits(const art::Event *const pEvent);

            void GetCommonCaloHitProperties(const gar::rec::CaloHit *const pCaloHit, PandoraApi::CaloHit::Parameters &caloHitParameters) const;

            int GetNLayersFromEdge(const  gar::rec::CaloHit *const pCaloHit) const;

            float GetMaximumRadius(const  gar::rec::CaloHit *const pCaloHit, const unsigned int symmetryOrder, const float phi0) const;

            const Settings                      m_settings;                         ///< The calo hit creator settings

            const pandora::Pandora &            m_pandora;                          ///< Reference to the pandora object to create calo hits

            float                               m_eCalBarrelLayerThickness;         ///< ECal barrel layer thickness
            float                               m_eCalEndCapLayerThickness;         ///< ECal endcap layer thickness

            CalorimeterHitVector                m_calorimeterHitVector;             ///< The calorimeter hit vector

            const geo::GeometryCore*            fGeo; //Geometry Manager
            gar::geo::BitFieldCoder const*      m_fieldDecoder;
        };

        //------------------------------------------------------------------------------------------------------------------------------------------

        inline const CalorimeterHitVector &CaloHitCreator::GetCalorimeterHitVector() const
        {
            return m_calorimeterHitVector;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------

        inline void CaloHitCreator::Reset()
        {
            m_calorimeterHitVector.clear();
        }

    }
}

#endif // #ifndef CALO_HIT_CREATOR_H