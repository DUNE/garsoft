#include "PfoCreator.h"

#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>

#include "art/Framework/Core/EDProducer.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace gar_pandora {

    PfoCreator::PfoCreator(const Settings& settings,
                           const pandora::Pandora* const pPandora,
                           const RotationTransformation* const pRotation)
      : m_settings(settings), m_pandora(*pPandora), m_rotation(*pRotation)
    {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    PfoCreator::~PfoCreator() {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode PfoCreator::CreateParticleFlowObjects(art::Event& pEvent)
    {
      PFParticleCollection outputParticles(new std::vector<rec::PFParticle>);
      ClusterCollection outputClusters(new std::vector<rec::Cluster>);

      const pandora::PfoList* pPandoraPfoList = NULL;
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS,
                               !=,
                               PandoraApi::GetCurrentPfoList(m_pandora, pPandoraPfoList));

      pandora::StringVector subDetectorNames;
      this->InitialiseSubDetectorNames(subDetectorNames);

      for (pandora::PfoList::const_iterator pIter = pPandoraPfoList->begin(),
                                            pIterEnd = pPandoraPfoList->end();
           pIter != pIterEnd;
           ++pIter) {
        const pandora::ParticleFlowObject* const pPandoraPfo(*pIter);
        gar::rec::PFParticle pReconstructedParticle;

        const bool hasTrack(!pPandoraPfo->GetTrackList().empty());
        const pandora::ClusterList& clusterList(pPandoraPfo->GetClusterList());

        float clustersTotalEnergy(0.f);
        pandora::CartesianVector referencePoint(0.f, 0.f, 0.f),
          clustersWeightedPosition(0.f, 0.f, 0.f);

        for (pandora::ClusterList::const_iterator cIter = clusterList.begin(),
                                                  cIterEnd = clusterList.end();
             cIter != cIterEnd;
             ++cIter) {
          const pandora::Cluster* const pPandoraCluster(*cIter);
          pandora::CaloHitList pandoraCaloHitList;
          pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
          pandoraCaloHitList.insert(pandoraCaloHitList.end(),
                                    pPandoraCluster->GetIsolatedCaloHitList().begin(),
                                    pPandoraCluster->GetIsolatedCaloHitList().end());

          pandora::FloatVector hitE, hitX, hitY, hitZ;
          gar::rec::Cluster pCluster;
          this->SetClusterSubDetectorEnergies(
            subDetectorNames, pandoraCaloHitList, hitE, hitX, hitY, hitZ);

          float clusterCorrectEnergy(0.f);
          this->SetClusterEnergyAndError(
            pPandoraPfo, pPandoraCluster, pCluster, clusterCorrectEnergy);

          pandora::CartesianVector clusterPosition(0.f, 0.f, 0.f);
          const unsigned int nHitsInCluster(pandoraCaloHitList.size());
          this->SetClusterPositionAndError(
            nHitsInCluster, hitE, hitX, hitY, hitZ, pCluster, clusterPosition);

          if (!hasTrack) {
            clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
            clustersTotalEnergy += clusterCorrectEnergy;
          }

          MF_LOG_DEBUG("PfoCreator::CreateParticleFlowObjects")
            << "Adding cluster " << &pCluster << " with energy " << pCluster.Energy();

          outputClusters->emplace_back(pCluster);
        }

        if (!hasTrack) {
          if (clustersTotalEnergy < std::numeric_limits<float>::epsilon()) {
            MF_LOG_DEBUG("PfoCreator::CreateParticleFlowObjects")
              << "invalid cluster energy " << clustersTotalEnergy;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
          }
          else {
            referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
          }
        }
        else {
          PANDORA_RETURN_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS,
            !=,
            this->CalculateTrackBasedReferencePoint(pPandoraPfo, referencePoint));
        }

        this->SetRecoParticleReferencePoint(referencePoint, pReconstructedParticle);
        this->AddTracksToRecoParticle(pPandoraPfo, pReconstructedParticle);
        this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, pReconstructedParticle);

        MF_LOG_DEBUG("PfoCreator::CreateParticleFlowObjects")
          << "Adding PFO " << &pReconstructedParticle << " with energy "
          << pReconstructedParticle.Energy() << " , mass " << pReconstructedParticle.Mass();

        // a hack to make clang happy -- make sure m_rotation is used in the code
        // this should never print

        if (pReconstructedParticle.Energy() < -1000)
          MF_LOG_DEBUG("PfoCreator::CreateParticleFlowObjects")
            << " Negative energy.  Message should not be printed. " << m_rotation.kAxisX;

        outputParticles->emplace_back(pReconstructedParticle);
      }

      pEvent.put(std::move(outputParticles));
      pEvent.put(std::move(outputClusters));

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::InitialiseSubDetectorNames(pandora::StringVector& subDetectorNames) const
    {
      subDetectorNames.push_back("ecal");
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::SetClusterSubDetectorEnergies(
      const pandora::StringVector& /* subDetectorNames */,
      const pandora::CaloHitList& pandoraCaloHitList,
      pandora::FloatVector& hitE,
      pandora::FloatVector& hitX,
      pandora::FloatVector& hitY,
      pandora::FloatVector& hitZ) const
    {
      for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(),
                                                hIterEnd = pandoraCaloHitList.end();
           hIter != hIterEnd;
           ++hIter) {
        const pandora::CaloHit* const pPandoraCaloHit(*hIter);
        gar::rec::CaloHit* const pCalorimeterHit =
          (gar::rec::CaloHit*)(pPandoraCaloHit->GetParentAddress());

        MF_LOG_DEBUG("PfoCreator::SetClusterSubDetectorEnergies")
          << " hit " << pCalorimeterHit << " energy " << pCalorimeterHit->Energy() << " position x "
          << pCalorimeterHit->Position()[0] << " y " << pCalorimeterHit->Position()[1] << " z "
          << pCalorimeterHit->Position()[2];

        const float caloHitEnergy(pCalorimeterHit->Energy());
        hitE.push_back(caloHitEnergy);
        hitX.push_back(pCalorimeterHit->Position()[0]);
        hitY.push_back(pCalorimeterHit->Position()[1]);
        hitZ.push_back(pCalorimeterHit->Position()[2]);
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject* const pPandoraPfo,
                                              const pandora::Cluster* const pPandoraCluster,
                                              gar::rec::Cluster& pCluster,
                                              float& clusterCorrectEnergy) const
    {
      const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) ||
                            (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
      clusterCorrectEnergy =
        (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(m_pandora) :
                      pPandoraCluster->GetCorrectedHadronicEnergy(m_pandora));

      if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

      const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm :
                                              m_settings.m_hadStochasticTerm);
      const float constantTerm(isEmShower ? m_settings.m_emConstantTerm :
                                            m_settings.m_hadConstantTerm);
      const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy +
                                        constantTerm * constantTerm) *
                              clusterCorrectEnergy);

      pCluster.setEnergy(clusterCorrectEnergy);
      pCluster.setEnergyError(energyError);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::SetClusterPositionAndError(const unsigned int nHitsInCluster,
                                                pandora::FloatVector& hitE,
                                                pandora::FloatVector& hitX,
                                                pandora::FloatVector& hitY,
                                                pandora::FloatVector& hitZ,
                                                gar::rec::Cluster& pCluster,
                                                pandora::CartesianVector& clusterPositionVec) const
    {
      util::ClusterShapes pClusterShapes(
        nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data());

      try {
        pCluster.setIPhi(std::atan2(pClusterShapes.getEigenVecInertia()[1],
                                    pClusterShapes.getEigenVecInertia()[0]));
        pCluster.setITheta(std::acos(pClusterShapes.getEigenVecInertia()[2]));
        pCluster.setPosition(pClusterShapes.getCenterOfGravity());
        //pCluster->setPositionError(pClusterShapes->getCenterOfGravityErrors());
        //pCluster->setDirectionError(pClusterShapes->getEigenVecInertiaErrors());
        clusterPositionVec.SetValues(pClusterShapes.getCenterOfGravity()[0],
                                     pClusterShapes.getCenterOfGravity()[1],
                                     pClusterShapes.getCenterOfGravity()[2]);
      }
      catch (...) {
        MF_LOG_WARNING("PfoCreator::SetClusterPositionAndError")
          << "unidentified exception caught.";
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode PfoCreator::CalculateTrackBasedReferencePoint(
      const pandora::ParticleFlowObject* const pPandoraPfo,
      pandora::CartesianVector& referencePoint) const
    {
      const pandora::TrackList& trackList(pPandoraPfo->GetTrackList());

      float totalTrackMomentumAtDca(0.f), totalTrackMomentumAtStart(0.f);
      pandora::CartesianVector referencePointAtDCAWeighted(0.f, 0.f, 0.f),
        referencePointAtStartWeighted(0.f, 0.f, 0.f);

      bool hasSiblings(false);
      for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end();
           tIter != tIterEnd;
           ++tIter) {
        const pandora::Track* const pPandoraTrack(*tIter);

        if (!this->IsValidParentTrack(pPandoraTrack, trackList)) continue;

        if (this->HasValidSiblingTrack(pPandoraTrack, trackList)) {
          // Presence of sibling tracks typically represents a conversion
          const pandora::CartesianVector& trackStartPoint(
            (pPandoraTrack->GetTrackStateAtStart()).GetPosition());
          const float trackStartMomentum(
            ((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
          referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
          totalTrackMomentumAtStart += trackStartMomentum;
          hasSiblings = true;
        }
        // else
        // {
        //     const gar::rec::Track *const pTrack = (gar::rec::Track*)(pPandoraTrack->GetParentAddress());
        //     const float z0(pPandoraTrack->GetZ0());
        //     pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);
        //
        //     intersectionPoint.SetValues(pTrack->getD0() * std::cos(pTrack->getPhi()), pTrack->getD0() * std::sin(pTrack->getPhi()), z0);
        //     const float trackMomentumAtDca((pPandoraTrack->GetMomentumAtDca()).GetMagnitude());
        //     referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
        //     totalTrackMomentumAtDca += trackMomentumAtDca;
        // }
      }

      if (hasSiblings) {
        if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon()) {
          MF_LOG_WARNING("PfoCreator::CalculateTrackBasedReferencePoint")
            << "invalid track momentum " << totalTrackMomentumAtStart;
          throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else {
          referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
        }
      }
      else {
        if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon()) {
          MF_LOG_WARNING("PfoCreator::CalculateTrackBasedReferencePoint")
            << "invalid track momentum " << totalTrackMomentumAtDca;
          throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else {
          referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
        }
      }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool PfoCreator::IsValidParentTrack(const pandora::Track* const pPandoraTrack,
                                        const pandora::TrackList& allTrackList) const
    {
      const pandora::TrackList& parentTrackList(pPandoraTrack->GetParentList());

      for (pandora::TrackList::const_iterator iter = parentTrackList.begin(),
                                              iterEnd = parentTrackList.end();
           iter != iterEnd;
           ++iter) {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
          continue;

        // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
        MF_LOG_WARNING("PfoCreator::IsValidParentTrack")
          << "mismatch in track relationship information, use information as available ";

        if (this->IsClosestTrackToIP(pPandoraTrack, allTrackList)) return true;

        return false;
      }

      // Ideal case: All parents are associated to same pfo
      return true;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool PfoCreator::HasValidSiblingTrack(const pandora::Track* const pPandoraTrack,
                                          const pandora::TrackList& allTrackList) const
    {
      const pandora::TrackList& siblingTrackList(pPandoraTrack->GetSiblingList());
      for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(),
                                              iterEnd = siblingTrackList.end();
           iter != iterEnd;
           ++iter) {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
          continue;

        // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in the list
        MF_LOG_WARNING("PfoCreator::HasValidSiblingTrack")
          << "mismatch in track relationship information, use information as available ";

        if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList)) return true;

        return false;
      }

      return true;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool PfoCreator::IsClosestTrackToIP(const pandora::Track* const pPandoraTrack,
                                        const pandora::TrackList& allTrackList) const
    {
      const pandora::Track* pClosestTrack(NULL);
      float closestTrackDisplacement(std::numeric_limits<float>::max());

      for (pandora::TrackList::const_iterator iter = allTrackList.begin(),
                                              iterEnd = allTrackList.end();
           iter != iterEnd;
           ++iter) {
        const pandora::Track* const pTrack(*iter);
        const float trialTrackDisplacement(
          pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

        if (trialTrackDisplacement < closestTrackDisplacement) {
          closestTrackDisplacement = trialTrackDisplacement;
          pClosestTrack = pTrack;
        }
      }

      return (pPandoraTrack == pClosestTrack);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool PfoCreator::AreAnyOtherSiblingsInList(const pandora::Track* const pPandoraTrack,
                                               const pandora::TrackList& allTrackList) const
    {
      const pandora::TrackList& siblingTrackList(pPandoraTrack->GetSiblingList());

      for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(),
                                              iterEnd = siblingTrackList.end();
           iter != iterEnd;
           ++iter) {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
          return true;
      }

      return false;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::SetRecoParticleReferencePoint(
      const pandora::CartesianVector& referencePoint,
      gar::rec::PFParticle& pReconstructedParticle) const
    {
      const float referencePointArray[3] = {
        referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
      pReconstructedParticle.setPosition(referencePointArray);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::AddTracksToRecoParticle(
      const pandora::ParticleFlowObject* const pPandoraPfo,
      gar::rec::PFParticle& /* pReconstructedParticle */) const
    {
      const pandora::TrackList& trackList(pPandoraPfo->GetTrackList());

      for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end();
           tIter != tIterEnd;
           ++tIter) {
        // const pandora::Track *const pTrack(*tIter);
        // pReconstructedParticle->addTrack((gar::rec::Track*)(pTrack->GetParentAddress()));
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void PfoCreator::SetRecoParticlePropertiesFromPFO(
      const pandora::ParticleFlowObject* const pPandoraPfo,
      gar::rec::PFParticle& pReconstructedParticle) const
    {
      const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(),
                                 pPandoraPfo->GetMomentum().GetY(),
                                 pPandoraPfo->GetMomentum().GetZ()};
      pReconstructedParticle.setMomentum(momentum);
      pReconstructedParticle.setEnergy(pPandoraPfo->GetEnergy());
      pReconstructedParticle.setMass(pPandoraPfo->GetMass());
      pReconstructedParticle.setCharge(pPandoraPfo->GetCharge());
      pReconstructedParticle.setType(pPandoraPfo->GetParticleId());
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    PfoCreator::Settings::Settings()
      : m_emStochasticTerm(0.17f)
      , m_emConstantTerm(0.01f)
      , m_hadStochasticTerm(0.30f)
      , m_hadConstantTerm(0.01f)
    {}
  }
}
