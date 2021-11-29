#include "TrackCreator.h"

#include "Pandora/PdgTable.h"
#include "GArObjects/Helix.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "MCCheater/BackTracker.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "RecoAlg/TrackPropagator.h"

namespace gar {
  namespace gar_pandora {

    TrackCreator::TrackCreator(const Settings &settings , const pandora::Pandora *const pPandora, const RotationTransformation *const pRotation)
      : m_settings(settings),
        m_pandora(*pPandora),
        m_rotation(*pRotation),
        m_trackVector(0),
        m_v0TrackList( TrackList() ),
        m_parentTrackList( TrackList() ),
        m_daughterTrackList( TrackList() )
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    TrackCreator::~TrackCreator()
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::CreateTracks(const art::Event &pEvent)
    {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectTracks(pEvent, m_settings.m_trackCollection, artTrkVector));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateTracks());

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::CreateTrackAssociations(const art::Event &pEvent)
    {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(pEvent));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(pEvent));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(pEvent));

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::ExtractKinks(const art::Event & /* pEvent */)
    {
      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::ExtractProngsAndSplits(const art::Event & /* pEvent */)
    {
      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::ExtractV0s(const art::Event & /* pEvent */)
    {
      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::CreateTracks()
    {
      for (TrackVector::const_iterator iter = artTrkVector.begin(), iterEnd = artTrkVector.end(); iter != iterEnd; ++iter)
        {
          art::Ptr<gar::rec::Track> artPtrTrack = *iter;
          const gar::rec::Track *pTrack = artPtrTrack.get();

          const unsigned int nTrackHits(static_cast<int>(pTrack->NHits()));
          if ((nTrackHits < m_settings.m_minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
            continue;

          PandoraApi::Track::Parameters trackParameters;

          const float *trackParams = pTrack->TrackParBeg(); //y, z, omega, phi, lambda
          const float omega = trackParams[2] / CLHEP::cm;
          const float d0 = std::sqrt(trackParams[0]*trackParams[0] + trackParams[1]*trackParams[1]) * CLHEP::cm;
          const float z0 = pTrack->Vertex()[0] * CLHEP::cm;

          trackParameters.m_pParentAddress = (void*)pTrack;
          trackParameters.m_d0 = d0;
          trackParameters.m_z0 = z0;

          trackParameters.m_particleId = (omega > 0) ? pandora::MU_PLUS : pandora::MU_MINUS;
          trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::MU_PLUS);

          if (std::numeric_limits<float>::epsilon() < std::fabs(omega))
            trackParameters.m_charge = static_cast<int>(omega / std::fabs(omega));

          try
            {
              this->GetTrackStates(pTrack, trackParameters);
              this->TrackReachesECAL(pTrack, trackParameters);
              this->DefineTrackPfoUsage(pTrack, trackParameters);

              MF_LOG_DEBUG("TrackCreator::CreateTracks()")
                << "Creating Track " << pTrack << "\n"
                << " starting at " << trackParameters.m_trackStateAtStart.Get() << "\n"
                << " ending at " << trackParameters.m_trackStateAtEnd.Get() << "\n"
                << " state at calo " << trackParameters.m_trackStateAtCalorimeter.Get() << "\n"
                << " with momentum " << trackParameters.m_momentumAtDca.Get() << "\n"
                << " magnitude " << trackParameters.m_momentumAtDca.Get().GetMagnitude() << " GeV\n"
                << " can form PFO " << trackParameters.m_canFormPfo.Get() << "\n"
                << " can form clusterless PFO " << trackParameters.m_canFormClusterlessPfo.Get() << "\n"
                << " time at calo " << trackParameters.m_timeAtCalorimeter.Get() << " ns";

              PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(m_pandora, trackParameters));
              m_trackVector.push_back(artPtrTrack);
            }
          catch (pandora::StatusCodeException &statusCodeException)
            {
              MF_LOG_DEBUG("TrackCreator::CreateTracks")
                << "Failed to extract a track: " << statusCodeException.ToString();
            }

        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void TrackCreator::GetTrackStates(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
    {
      const float *trackParams = pTrack->TrackParEnd(); //y, z, omega, phi, lambda
      const float omega = trackParams[2] / CLHEP::cm;

      //Track momentum
      const double pt(m_settings.m_bField * 2.99792e-4 / std::fabs(omega));
      const pandora::CartesianVector momentum = pandora::CartesianVector(std::cos(trackParams[3]), std::sin(trackParams[3]), std::tan(trackParams[4])) * pt;
      const pandora::CartesianVector newmomentum = m_rotation.MakeRotation(momentum);
      trackParameters.m_momentumAtDca = newmomentum;

      //Track state at begin and end
      const pandora::CartesianVector startPosition(pTrack->Vertex()[0] * CLHEP::cm, pTrack->Vertex()[1] * CLHEP::cm, pTrack->Vertex()[2] * CLHEP::cm);
      const pandora::CartesianVector endPosition(pTrack->End()[0] * CLHEP::cm, pTrack->End()[1] * CLHEP::cm, pTrack->End()[2] * CLHEP::cm);
      const pandora::CartesianVector newstartPosition = m_rotation.MakeRotation(startPosition);
      const pandora::CartesianVector newendPosition = m_rotation.MakeRotation(endPosition);

      trackParameters.m_trackStateAtStart = pandora::TrackState(newstartPosition, momentum);
      trackParameters.m_trackStateAtEnd = pandora::TrackState(newendPosition, momentum);

      pandora::CartesianVector calorimeterPosition(0., 0., 0.);
      float minGenericTime(0.);
      this->CalculateTrackStateAtCalo(pTrack, calorimeterPosition);
      this->CalculateTimeAtCalo(pTrack, minGenericTime);

      const pandora::CartesianVector newcalorimeterPosition = m_rotation.MakeRotation(calorimeterPosition);
      trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(newcalorimeterPosition, momentum);
      trackParameters.m_isProjectedToEndCap = ((std::fabs(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetX()) < m_settings.m_eCalEndCapInnerZ) ? false : true);

      // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
      const float particleMass(trackParameters.m_mass.Get());
      const float particleEnergy(std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
      trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 299.792f;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool TrackCreator::PassesQualityCuts(const gar::rec::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
    {
      if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp)
        return false;

      if (std::fabs(pTrack->TrackParEnd()[2] / CLHEP::cm) < std::numeric_limits<float>::epsilon())
        {
          MF_LOG_ERROR("TrackCreator::PassesQualityCuts")
            << "Track has Omega = 0 ";
          return false;
        }

      // Check momentum uncertainty is reasonable to use track
      const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
      const float sigmaPOverP(std::sqrt(pTrack->CovMatEndPacked()[5]) / std::fabs( pTrack->TrackParEnd()[2] / CLHEP::cm ));

      if (sigmaPOverP > m_settings.m_maxTrackSigmaPOverP)
        {
          MF_LOG_WARNING("TrackCreator::PassesQualityCuts")
            << " Dropping track : " << momentumAtDca.GetMagnitude()
            << " +- " << sigmaPOverP * (momentumAtDca.GetMagnitude())
            << " chi2 = " <<  pTrack->ChisqForward();
          return false;
        }

      return true;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void TrackCreator::TrackReachesECAL(const gar::rec::Track *const /* pTrack */, PandoraApi::Track::Parameters &trackParameters) const
    {
      trackParameters.m_reachesCalorimeter = true;
      return;

      //Do a check later on with the track state at the calo!
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void TrackCreator::DefineTrackPfoUsage(const gar::rec::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
    {
      bool canFormPfo(false);
      bool canFormClusterlessPfo(false);

      if (trackParameters.m_reachesCalorimeter.Get() && !this->IsParent(pTrack))
        {
          const float *trackParams = pTrack->TrackParEnd(); //y, z, omega, phi, lambda
          const float d0(std::sqrt(trackParams[0]*trackParams[0] + trackParams[1]+trackParams[1]) * CLHEP::cm);
          const float z0(std::fabs(pTrack->Vertex()[0]) * CLHEP::cm);

          if (this->PassesQualityCuts(pTrack, trackParameters))
            {
              const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
              const bool isV0(this->IsV0(pTrack));
              const bool isDaughter(this->IsDaughter(pTrack));

              // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
              if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut))
                {
                  canFormPfo = true;
                }
              else if (isV0 || isDaughter)
                {
                  canFormPfo = true;
                }

              // Decide whether track can be used to form a charged PFO, even if track fails to be associated with a pandora cluster
              const float particleMass(trackParameters.m_mass.Get());
              const float trackEnergy(std::sqrt(momentumAtDca.GetMagnitudeSquared() + particleMass * particleMass));

              if (trackEnergy < m_settings.m_unmatchedVertexTrackMaxEnergy)
                {
                  if ((d0 < m_settings.m_d0UnmatchedVertexTrackCut) && (z0 < m_settings.m_z0UnmatchedVertexTrackCut))
                    {
                      canFormClusterlessPfo = true;
                    }
                  else if (isV0 || isDaughter)
                    {
                      canFormClusterlessPfo = true;
                    }
                }
            }
          else if (this->IsDaughter(pTrack) || this->IsV0(pTrack))
            {
              MF_LOG_DEBUG("TrackCreator::DefineTrackPfoUsage")
                << "Recovering daughter or v0 track "
                << trackParameters.m_momentumAtDca.Get().GetMagnitude();
              canFormPfo = true;
            }

          MF_LOG_DEBUG("TrackCreator::DefineTrackPfoUsage")
            << " -- track canFormPfo = " << canFormPfo
            << " -  canFormClusterlessPfo = " << canFormClusterlessPfo;

        }

      trackParameters.m_canFormPfo = canFormPfo;
      trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;

      return;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------
    void TrackCreator::CalculateTrackStateAtCalo(const gar::rec::Track *const pTrack, pandora::CartesianVector &posAtCalo) const
    {
      pandora::CartesianVector bestECalProjection(0., 0., 0);
      float fbestECalProjection[3] = {0., 0., 0.};    float otherECALProjection[3];
      //Use our TrackPropagator
      //First propagate to the barrel; 2nd extrapolation isn't used here.  Guess it could be.
      int result = util::TrackPropagator::PropagateToCylinder(pTrack->TrackParEnd(), pTrack->End(), m_settings.m_eCalBarrelInnerR / CLHEP::cm, m_settings.m_GArCenterY / CLHEP::cm, m_settings.m_GArCenterZ / CLHEP::cm, fbestECalProjection,otherECALProjection, m_settings.m_eCalEndCapInnerZ / CLHEP::cm );
      bestECalProjection.SetValues(fbestECalProjection[0] * CLHEP::cm, fbestECalProjection[1] * CLHEP::cm, fbestECalProjection[2] * CLHEP::cm);

      if( result != 0 ) {
        //Propagate to the Endcap
        float endcapProjection[3] = {0., 0., 0.};
        int result_local = util::TrackPropagator::PropagateToX( pTrack->TrackParEnd(), pTrack->End(), (pTrack->End()[0] > 0) ? m_settings.m_eCalEndCapInnerZ / CLHEP::cm : -m_settings.m_eCalEndCapInnerZ / CLHEP::cm, endcapProjection, m_settings.m_eCalBarrelInnerR / CLHEP::cm );

        if(result_local == 0) {
          bestECalProjection.SetValues(endcapProjection[0] * CLHEP::cm, endcapProjection[1] * CLHEP::cm, endcapProjection[2] * CLHEP::cm);
        }
      }

      if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

      posAtCalo = bestECalProjection;

      return;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------
    void TrackCreator::CalculateTimeAtCalo(const gar::rec::Track *const pTrack, float &timeAtCalo) const
    {
      const float *trackParams = pTrack->TrackParEnd(); //y, z, omega, phi, lambda
      const float omega = trackParams[2] / CLHEP::cm;
      const float tanl = std::tan(trackParams[4]);
      const float d0(std::sqrt(trackParams[0]*trackParams[0] + trackParams[1]+trackParams[1]) * CLHEP::cm);
      const float z0(std::fabs(pTrack->Vertex()[0]) * CLHEP::cm);

      const gar_content::Helix helix(trackParams[3], d0, z0, omega, tanl, m_settings.m_bField);
      const pandora::CartesianVector &referencePoint(helix.GetReferencePoint());

      // First project to endcap
      float minGenericTime(std::numeric_limits<float>::max());

      pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
      const int signPx((helix.GetMomentum().GetX() > 0.f) ? 1 : -1);
      (void) helix.GetPointInX(static_cast<float>(signPx) * m_settings.m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);

      // Then project to barrel surface(s)
      pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
      if (m_settings.m_eCalBarrelInnerSymmetry > 0)
        {
          // Polygon
          float twopi_n = 2. * M_PI / (static_cast<float>(m_settings.m_eCalBarrelInnerSymmetry));

          for (int i = 0; i < m_settings.m_eCalBarrelInnerSymmetry; ++i)
            {
              float genericTime(std::numeric_limits<float>::max());
              const float phi(twopi_n * static_cast<float>(i) + m_settings.m_eCalBarrelInnerPhi0);

              const pandora::StatusCode statusCode(helix.GetPointInZY(m_settings.m_eCalBarrelInnerR * std::cos(phi), m_settings.m_eCalBarrelInnerR * std::sin(phi), std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

              if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
                {
                  minGenericTime = genericTime;
                  bestECalProjection = barrelProjection;
                }
            }
        }
      else
        {
          // Cylinder
          float genericTime(std::numeric_limits<float>::max());
          const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_settings.m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

          if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
              minGenericTime = genericTime;
              bestECalProjection = barrelProjection;
            }
        }

      if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

      timeAtCalo = minGenericTime;

      return;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode TrackCreator::CollectTracks(const art::Event &pEvent, const std::string &label, TrackVector &trkVector)
    {
      auto theTrk = pEvent.getHandle< RawTrackVector >(label);
      if (!theTrk)
        {
          MF_LOG_DEBUG("TrackCreator::CreateTracks") << "  Failed to find tracks... " << std::endl;
          return pandora::STATUS_CODE_NOT_FOUND;
        }

      MF_LOG_DEBUG("TrackCreator::CreateTracks") << "  Found: " << theTrk->size() << " tracks " << std::endl;

      for (unsigned int i = 0; i < theTrk->size(); ++i)
        {
          const art::Ptr<gar::rec::Track> trk(theTrk, i);
          trkVector.push_back(trk);
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    TrackCreator::Settings::Settings()
      : m_trackCollection( "" ),
        m_V0Collection( "" ),
        m_minTrackHits(3),
        m_maxTrackHits(5000.f),
        m_d0TrackCut(2500.f),
        m_z0TrackCut(3000.f),
        m_unmatchedVertexTrackMaxEnergy(0.1f),
        m_d0UnmatchedVertexTrackCut(2500.f),
        m_z0UnmatchedVertexTrackCut(3000.f),
        m_minTrackECalDistanceFromIp(100.f),
        m_maxTrackSigmaPOverP(0.15f),
        m_bField(0.f),
        m_eCalBarrelInnerSymmetry(0),
        m_eCalBarrelInnerPhi0(0.f),
        m_eCalBarrelInnerR(0.f),
        m_eCalEndCapInnerZ(0.f),
        m_GArCenterY(0.f),
        m_GArCenterZ(0.f)
    {
    }
  }
}
