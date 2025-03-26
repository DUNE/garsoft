#ifndef PFO_CREATOR_H
#define PFO_CREATOR_H 1

#include "art/Framework/Principal/Event.h"

#include "Api/PandoraApi.h"

#include "RotationTransformation.h"

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/PFParticle.h"
#include "ReconstructionDataProducts/Track.h"

#include "RecoAlg/ClusterShapes.h"

namespace gar {
  namespace gar_pandora {

    class PfoCreator {
    public:
      typedef std::unique_ptr<std::vector<rec::PFParticle>> PFParticleCollection;
      typedef std::unique_ptr<std::vector<rec::Cluster>> ClusterCollection;

      class Settings {
      public:
        Settings();

        float m_emStochasticTerm = 0;  ///< The stochastic term for EM shower energy resolution
        float m_emConstantTerm = 0;    ///< The constant term for EM shower energy resolution
        float m_hadStochasticTerm = 0; ///< The stochastic term for HAD shower energy resolution
        float m_hadConstantTerm = 0;   ///< The constant term for HAD shower energy resolution
      };

      PfoCreator(const Settings& settings,
                 const pandora::Pandora* const pPandora,
                 const RotationTransformation* const pRotation);

      ~PfoCreator();

      pandora::StatusCode CreateParticleFlowObjects(art::Event& pEvent);

    private:
      /**
            *  @brief  initialise sub detector name strings
            *
            *  @param  subDetectorNames to receive the list of sub detector names
            */
      void InitialiseSubDetectorNames(pandora::StringVector& subDetectorNames) const;

      /**
            *  @brief  Set sub detector energies for a cluster
            *
            *  @param  subDetectorNames the list of sub detector names
            *  @param  pLcioCluster the address of the lcio cluster to be set sub detector energies
            *  @param  pandoraCaloHitList the pandora calorimeter hit list
            *  @param  hitE the vector to receive the energy of hits
            *  @param  hitX the vector to receive the x position of hits
            *  @param  hitY the vector to receive the y position of hits
            *  @param  hitZ the vector to receive the z position of hits
            */
      void SetClusterSubDetectorEnergies(const pandora::StringVector& subDetectorNames,
                                         const pandora::CaloHitList& pandoraCaloHitList,
                                         pandora::FloatVector& hitE,
                                         pandora::FloatVector& hitX,
                                         pandora::FloatVector& hitY,
                                         pandora::FloatVector& hitZ) const;

      /**
            *  @brief  Set cluster energies and errors
            *
            *  @param  pPandoraPfo the address of the pandora pfo
            *  @param  pPandoraCluster the address of the pandora cluster
            *  @param  pLcioCluster the address of the lcio cluster to be set energies and erros
            *  @param  clusterCorrectEnergy a number to receive the cluster correct energy
            */
      void SetClusterEnergyAndError(const pandora::ParticleFlowObject* const pPandoraPfo,
                                    const pandora::Cluster* const pPandoraCluster,
                                    gar::rec::Cluster& pCluster,
                                    float& clusterCorrectEnergy) const;

      /**
            *  @brief  Set cluster position, errors and other shape info, by calculating culster shape first
            *
            *  @param  nHitsInCluster number of hits in cluster
            *  @param  hitE the vector of the energy of hits
            *  @param  hitX the vector of the x position of hits
            *  @param  hitY the vector of the y position of hits
            *  @param  hitZ the vector of the z position of hits
            *  @param  pLcioCluster the lcio cluster to be set positions and errors
            *  @param  clusterPosition a CartesianVector to receive the cluster position
            */
      void SetClusterPositionAndError(const unsigned int nHitsInCluster,
                                      pandora::FloatVector& hitE,
                                      pandora::FloatVector& hitX,
                                      pandora::FloatVector& hitY,
                                      pandora::FloatVector& hitZ,
                                      gar::rec::Cluster& pCluster,
                                      pandora::CartesianVector& clusterPositionVec) const;

      /**
            *  @brief  Calculate reference point for pfo with tracks
            *
            *  @param  pPandoraPfo the address of the pandora pfo
            *  @param  referencePoint a CartesianVector to receive the reference point
            */
      pandora::StatusCode CalculateTrackBasedReferencePoint(
        const pandora::ParticleFlowObject* const pPandoraPfo,
        pandora::CartesianVector& referencePoint) const;

      /**
            *  @brief  Set reference point of the reconstructed particle
            *
            *  @param  referencePoint a CartesianVector of the reference point
            *  @param  pReconstructedParticle the address of the reconstructed particle to be reference point
            */
      void SetRecoParticleReferencePoint(const pandora::CartesianVector& referencePoint,
                                         gar::rec::PFParticle& pReconstructedParticle) const;

      /**
            *  @brief  Add tracks to reconstructed particle
            *
            *  @param  pPandoraPfo the address of the pandora pfo
            *  @param  pReconstructedParticle the address of the reconstructed particle to be added tracks
            */
      void AddTracksToRecoParticle(const pandora::ParticleFlowObject* const pPandoraPfo,
                                   gar::rec::PFParticle& pReconstructedParticle) const;

      /**
            *  @brief  Set properties of reconstructed particle from pandora pfo
            *
            *  @param  pPandoraPfo the address of the pandora pfo
            *  @param  pReconstructedParticle the address of the reconstructed particle to be set properties
            */
      void SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject* const pPandoraPfo,
                                            gar::rec::PFParticle& pReconstructedParticle) const;

      /**
            *  @brief  Whether parent and daughter tracks are associated with the same pfo
            *
            *  @param  pPandoraTrack the address of the pandora track
            *  @param  allTrackList list of all tracks associated with reconstructed particle
            *
            *  @return boolean
            */
      bool IsValidParentTrack(const pandora::Track* const pPandoraTrack,
                              const pandora::TrackList& allTrackList) const;

      /**
            *  @brief  Whether sibling tracks are associated with the same pfo
            *
            *  @param  pPandoraTrack the address of the pandora track
            *  @param  allTrackList list of all tracks associated with reconstructed particle
            *
            *  @return boolean
            */
      bool HasValidSiblingTrack(const pandora::Track* const pPandoraTrack,
                                const pandora::TrackList& allTrackList) const;

      /**
            *  @brief  Whether the track is the closest (of those associated with the same pfo) to the interaction point
            *
            *  @param  pPandoraTrack the address of the pandora track
            *  @param  allTrackList list of all tracks associated to reconstructed particle
            *
            *  @return boolean
            */
      bool IsClosestTrackToIP(const pandora::Track* const pPandoraTrack,
                              const pandora::TrackList& allTrackList) const;

      /**
            *  @brief  Whether at least one track sibling track is associated to the reconstructed particle
            *
            *  @param  pPandoraTrack the address of the pandora track
            *  @param  allTrackList list of all tracks associated to reconstructed particle
            *
            *  @return boolean
            */
      bool AreAnyOtherSiblingsInList(const pandora::Track* const pPandoraTrack,
                                     const pandora::TrackList& allTrackList) const;

      const Settings m_settings; ///< The pfo creator settings
      const pandora::Pandora&
        m_pandora; ///< Reference to the pandora object from which to extract the pfos
      const RotationTransformation& m_rotation;
    };
  }
}

#endif // #ifndef PFO_CREATOR_H
