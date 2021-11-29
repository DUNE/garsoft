#include "art/Framework/Principal/Handle.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "CoreUtils/ServiceUtil.h"

#include "CaloHitCreator.h"
#include "MCParticleCreator.h"
#include "TrackCreator.h"

#include "MCCheater/BackTracker.h"

#include "Pandora/PdgTable.h"

#include <cmath>
#include <limits>

namespace gar {
  namespace gar_pandora {

    MCParticleCreator::MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora, const RotationTransformation *const pRotation)
      : m_settings(settings),
        m_pandora(*pPandora),
        m_rotation(*pRotation)
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();
      auto const *magFieldService = gar::providerFrom<mag::MagneticFieldService>();
      G4ThreeVector zerovec(0, 0, 0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);
      m_bField = magfield[0]; //x component at (0, 0, 0)
      m_origin[0] = fGeo->GetOriginX();
      m_origin[1] = fGeo->GetOriginY();
      m_origin[2] = fGeo->GetOriginZ();
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    MCParticleCreator::~MCParticleCreator()
    {
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CreateMCParticles(const art::Event &pEvent)
    {
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectMCParticles(pEvent, m_settings.m_geantModuleLabel, artMCParticleVector));
      if (!m_settings.m_generatorModuleLabel.empty())
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectGeneratorMCParticles(pEvent, m_settings.m_generatorModuleLabel, generatorArtMCParticleVector));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CollectMCParticles(pEvent, m_settings.m_geantModuleLabel, artMCTruthToMCParticles, artMCParticlesToMCTruth));
      PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CreateMCParticles());

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CreateMCParticles() const
    {
      MCParticleMap particleMap;
      for (MCParticlesToMCTruth::const_iterator iter = artMCParticlesToMCTruth.begin(), iterEnd = artMCParticlesToMCTruth.end(); iter != iterEnd; ++iter)
        {
          const art::Ptr<simb::MCParticle> particle = iter->first;
          particleMap[particle->TrackId()] = particle;
        }

      int neutrinoCounter(0);
      for (MCTruthToMCParticles::const_iterator iter1 = artMCTruthToMCParticles.begin(), iterEnd1 = artMCTruthToMCParticles.end(); iter1 != iterEnd1; ++iter1)
        {
          const art::Ptr<simb::MCTruth> truth = iter1->first;
          if (truth->NeutrinoSet())
            {
              const simb::MCNeutrino neutrino(truth->GetNeutrino());
              ++neutrinoCounter;

              const pandora::CartesianVector momentum(neutrino.Nu().Px(), neutrino.Nu().Py(), neutrino.Nu().Pz());
              const pandora::CartesianVector vertex( (neutrino.Nu().Vx() - m_origin[0]) * CLHEP::cm, (neutrino.Nu().Vy()  - m_origin[1]) * CLHEP::cm, (neutrino.Nu().Vz()  - m_origin[2]) * CLHEP::cm);
              const pandora::CartesianVector endpoint( (neutrino.Nu().EndX() - m_origin[0]) * CLHEP::cm, (neutrino.Nu().EndY() - m_origin[1]) * CLHEP::cm, (neutrino.Nu().EndZ() - m_origin[2]) * CLHEP::cm);

              const pandora::CartesianVector newmomentum = m_rotation.MakeRotation(momentum);
              const pandora::CartesianVector newvertex = m_rotation.MakeRotation(vertex);
              const pandora::CartesianVector newendpoint = m_rotation.MakeRotation(endpoint);

              PandoraApi::MCParticle::Parameters mcParticleParameters;
              mcParticleParameters.m_energy = neutrino.Nu().E();
              mcParticleParameters.m_momentum = newmomentum;
              mcParticleParameters.m_vertex = newvertex;
              mcParticleParameters.m_endpoint = newendpoint;
              mcParticleParameters.m_particleId = neutrino.Nu().PdgCode();
              mcParticleParameters.m_mcParticleType = pandora::MC_3D;
              mcParticleParameters.m_pParentAddress = &neutrino;

              PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

              // Loop over associated particles
              const MCParticleVector &particleVector = iter1->second;
              for (MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
                {
                  const art::Ptr<simb::MCParticle> particle = *iter2;
                  // Mother/Daughter Links
                  if (particle->Mother() == 0)
                    {
                      try
                        {
                          PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora, &neutrino, particle.get()));
                        }
                      catch (const pandora::StatusCodeException &)
                        {
                          MF_LOG_WARNING("MCParticleCreator") << "CreatePandoraMCParticles - unable to create mc particle relationship, invalid information supplied " << std::endl;
                          continue;
                        }
                    }
                }
            }
        }

      MF_LOG_DEBUG("MCParticleCreator") << " Number of Pandora neutrinos: " << neutrinoCounter << std::endl;

      for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI)
        {
          const art::Ptr<simb::MCParticle> pMcParticle = iterI->second;

          // Lookup position and kinematics at start and end points
          const pandora::CartesianVector momentum(pMcParticle->Px(), pMcParticle->Py(), pMcParticle->Pz());
          const pandora::CartesianVector vertex( (pMcParticle->Vx() - m_origin[0]) * CLHEP::cm, (pMcParticle->Vy() - m_origin[1]) * CLHEP::cm, (pMcParticle->Vz() - m_origin[2]) * CLHEP::cm);
          const pandora::CartesianVector endpoint( (pMcParticle->EndX() - m_origin[0]) * CLHEP::cm, (pMcParticle->EndY() - m_origin[1]) * CLHEP::cm, (pMcParticle->EndZ() - m_origin[2]) * CLHEP::cm);

          const pandora::CartesianVector newmomentum = m_rotation.MakeRotation(momentum);
          const pandora::CartesianVector newvertex = m_rotation.MakeRotation(vertex);
          const pandora::CartesianVector newendpoint = m_rotation.MakeRotation(endpoint);

          PandoraApi::MCParticle::Parameters mcParticleParameters;
          mcParticleParameters.m_energy = pMcParticle->E();
          mcParticleParameters.m_particleId = pMcParticle->PdgCode();
          mcParticleParameters.m_mcParticleType = pandora::MC_3D;
          mcParticleParameters.m_pParentAddress = pMcParticle.get();
          mcParticleParameters.m_momentum = newmomentum;
          mcParticleParameters.m_vertex = newvertex;
          mcParticleParameters.m_endpoint = newendpoint;

          MF_LOG_DEBUG("MCParticleCreator") << " Adding MC Particle with parameters "
                                            << " mcParticleParameters.m_energy = " << mcParticleParameters.m_energy.Get()
                                            << " mcParticleParameters.m_particleId = " << mcParticleParameters.m_particleId.Get()
                                            << " mcParticleParameters.m_mcParticleType = " << mcParticleParameters.m_mcParticleType.Get()
                                            << " mcParticleParameters.m_pParentAddress = " << mcParticleParameters.m_pParentAddress.Get()
                                            << " mcParticleParameters.m_momentum = " << mcParticleParameters.m_momentum.Get()
                                            << " mcParticleParameters.m_vertex = " << mcParticleParameters.m_vertex.Get()
                                            << " mcParticleParameters.m_endpoint = " << mcParticleParameters.m_endpoint.Get();

          try
            {
              MF_LOG_DEBUG("MCParticleCreator::CreateMCParticles")
                << " Creating mc particle " << pMcParticle.get()
                << " of pdg " << pMcParticle->PdgCode()
                << " with TrackID " << pMcParticle->TrackId()
                << " and energy " << pMcParticle->E() << " GeV";

              PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(m_pandora, mcParticleParameters));

              // Create parent-daughter relationships
              const int id_mother(pMcParticle->Mother());
              MCParticleMap::const_iterator iterJ = particleMap.find(id_mother);
              if (iterJ != particleMap.end())
                {
                  try
                    {
                      MF_LOG_DEBUG("MCParticleCreator::CreateMCParticles")
                        << " Adding daughter relation " << iterJ->second.get()
                        << " to mc particle " << pMcParticle.get();

                      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(m_pandora, iterJ->second.get(), pMcParticle.get()));
                    }
                  catch (const pandora::StatusCodeException &)
                    {
                      MF_LOG_WARNING("MCParticleCreator") << "CreatePandoraMCParticles - Unable to create mc particle relationship, invalid information supplied " << std::endl;
                      continue;
                    }
                }
            }
          catch (const pandora::StatusCodeException &)
            {
              MF_LOG_WARNING("MCParticleCreator") << "CreatePandoraMCParticles - Unable to create MCParticle " << std::endl;
              continue;
            }
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const TrackVector &trackVector) const
    {
      MCParticleMap particleMap;
      for (MCParticlesToMCTruth::const_iterator iter = artMCParticlesToMCTruth.begin(), iterEnd = artMCParticlesToMCTruth.end(); iter != iterEnd; ++iter)
        {
          const art::Ptr<simb::MCParticle> particle = iter->first;
          particleMap[particle->TrackId()] = particle;
        }

      for( auto const& itr : trackVector )
        {
          try
            {
              const gar::rec::Track *pTrack = itr.get();

              const float *trackParams = pTrack->TrackParEnd(); //y, z, omega, phi, lambda
              const float omega = trackParams[2] / CLHEP::cm;
              const float d0 = std::sqrt(trackParams[0]*trackParams[0] + trackParams[1]*trackParams[1]) * CLHEP::cm;
              const float z0 = pTrack->End()[0] * CLHEP::cm;

              const pandora::Helix helixFit(trackParams[3], d0, z0, omega, std::tan(trackParams[4]), m_bField);
              const float recoMomentum(helixFit.GetMomentum().GetMagnitude());

              // Use momentum magnitude to identify best mc particle
              simb::MCParticle *pBestMCParticle = nullptr;
              float bestDeltaMomentum(std::numeric_limits<float>::max());

              //Loop over the MCParticles
              for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI)
                {
                  simb::MCParticle *pMCParticle = const_cast<simb::MCParticle*>(iterI->second.get());

                  if (nullptr == pMCParticle)
                    continue;

                  const pandora::CartesianVector momentum(pMCParticle->Px(), pMCParticle->Py(), pMCParticle->Pz());
                  const pandora::CartesianVector newmomentum = m_rotation.MakeRotation(momentum);
                  const float trueMomentum(newmomentum.GetMagnitude());

                  const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));

                  if (deltaMomentum < bestDeltaMomentum)
                    {
                      pBestMCParticle = pMCParticle;
                      bestDeltaMomentum = deltaMomentum;
                    }
                }

              if (nullptr == pBestMCParticle)
                continue;

              MF_LOG_DEBUG("MCParticleCreator::CreateTrackToMCParticleRelationships")
                << "Found MCParticle " << pBestMCParticle
                << " associated to track " << pTrack
                << " with best delta momentum " << bestDeltaMomentum;

              PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(m_pandora, pTrack, pBestMCParticle));
            }
          catch (pandora::StatusCodeException &statusCodeException)
            {
              MF_LOG_ERROR("MCParticleCreator::CreateTrackToMCParticleRelationships")
                << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString();
            }
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const CalorimeterHitVector &calorimeterHitVector) const
    {
      cheat::BackTrackerCore const* bt = gar::providerFrom<cheat::BackTracker>();
      std::map< eveLoc, std::vector< art::Ptr<gar::rec::CaloHit> > > eveCaloHitMap;

      MCParticleMap particleMap;
      for (MCParticlesToMCTruth::const_iterator iter = artMCParticlesToMCTruth.begin(), iterEnd = artMCParticlesToMCTruth.end(); iter != iterEnd; ++iter)
        {
          const art::Ptr<simb::MCParticle> particle = iter->first;
          particleMap[particle->TrackId()] = particle;
        }

      // loop over all hits and fill in the map
      for( auto const& itr : calorimeterHitVector )
        {
          std::vector<gar::cheat::CalIDE> eveides = bt->CaloHitToCalIDEs(itr);

          // loop over all eveides for this hit
          for(size_t ieve = 0; ieve < eveides.size(); ieve++) {
            MF_LOG_DEBUG("MCParticleCreator::CreateCaloHitToMCParticleRelationships")
              << " Found eveID " << eveides[ieve].trackID
              << " associated to art hit " << itr;

            if(eveides[ieve].energyFrac < 0.1) continue;

            eveLoc el(eveides[ieve].trackID);
            eveCaloHitMap[el].push_back(itr);
          } // end loop over eve IDs for this hit
        }// end loop over hits

      for(auto const &hitMapItr : eveCaloHitMap)
        {
          MF_LOG_DEBUG("MCParticleCreator::CreateCaloHitToMCParticleRelationships")
            << " Trying to find mcp associated to eveID " << hitMapItr.first.GetEveID();

          const simb::MCParticle *part = bt->TrackIDToParticle(hitMapItr.first.GetEveID());

          if( nullptr == part ) {
            MF_LOG_WARNING("MCParticleCreator::CreateCaloHitToMCParticleRelationships")
              << "Cannot find MCParticle for eveid: " << hitMapItr.first.GetEveID();
            continue;
          }

          // const int eveID = hitMapItr.first.GetEveID();
          const int trackID = part->TrackId();
          // const float partE = part->E();

          for (MCParticleMap::const_iterator iterI = particleMap.begin(), iterEndI = particleMap.end(); iterI != iterEndI; ++iterI)
            {
              simb::MCParticle *pMCParticleCalo = const_cast<simb::MCParticle*>(iterI->second.get());

              if (nullptr == pMCParticleCalo)
                continue;
              if (trackID != pMCParticleCalo->TrackId())
                continue;

              //Need to get the original monte carlo particle that created this one (mother <- daughters)
              const simb::MCParticle *pMCPrimary = MCParticleCreator::GetFinalStateMCParticle(particleMap, pMCParticleCalo);
              if (nullptr == pMCPrimary)
                continue;

              for(auto const& itr : hitMapItr.second)
                {
                  const gar::rec::CaloHit *hit = itr.get();
                  try
                    {
                      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(m_pandora, hit, pMCPrimary, hit->Energy()));
                    }
                  catch (const pandora::StatusCodeException &)
                    {
                      MF_LOG_WARNING("MCParticleCreator") << "CreateCaloHitToMCParticleRelationships - unable to create calo hit to mc particle relationship, invalid information supplied " << std::endl;
                      continue;
                    }
                }
            }
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CollectMCParticles(const art::Event &pEvent, const std::string &label, MCParticleVector &particleVector)
    {
      auto theParticles = pEvent.getHandle< RawMCParticleVector >(label);
      if (!theParticles)
        {
          MF_LOG_WARNING("MCParticleCreator") << "  Failed to find MC particles for label " << label << std::endl;
          return pandora::STATUS_CODE_NOT_FOUND;
        }

      MF_LOG_DEBUG("MCParticleCreator") << "  Found: " << theParticles->size() << " MC particles " << std::endl;

      for (unsigned int i = 0; i < theParticles->size(); ++i)
        {
          const art::Ptr<simb::MCParticle> particle(theParticles, i);
          particleVector.push_back(particle);
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CollectGeneratorMCParticles(const art::Event &pEvent, const std::string &label, RawMCParticleVector &particleVector)
    {
      auto mcTruthBlocks = pEvent.getHandle< std::vector<simb::MCTruth> >(label);
      if (!mcTruthBlocks)
        {
          MF_LOG_WARNING("MCParticleCreator") << "  Failed to find MC Truth for generator " << label << std::endl;
          return pandora::STATUS_CODE_NOT_FOUND;
        }

      MF_LOG_DEBUG("MCParticleCreator") << "  Found: " << mcTruthBlocks->size() << " MC truth blocks " << std::endl;

      if (mcTruthBlocks->size() != 1)
        throw cet::exception("MCParticleCreator") << " PandoraCollector::CollectGeneratorMCParticles --- Unexpected number of MC truth blocks ";

      const art::Ptr<simb::MCTruth> mcTruth(mcTruthBlocks, 0);
      for (int i = 0; i < mcTruth->NParticles(); ++i)
        {
          particleVector.push_back(mcTruth->GetParticle(i));
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode MCParticleCreator::CollectMCParticles(const art::Event &pEvent, const std::string &label, MCTruthToMCParticles &truthToParticles, MCParticlesToMCTruth &particlesToTruth)
    {
      auto theParticles = pEvent.getHandle< RawMCParticleVector >(label);
      if (!theParticles)
        {
          MF_LOG_WARNING("MCParticleCreator") << "  Failed to find MC particles for label " << label << std::endl;
          return pandora::STATUS_CODE_NOT_FOUND;
        }

      MF_LOG_DEBUG("MCParticleCreator") << "  Found: " << theParticles->size() << " MC particles " << std::endl;

      art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, pEvent, label);

      for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
        {
          const art::Ptr<simb::MCParticle> particle(theParticles, i);
          const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
          truthToParticles[truth].push_back(particle);
          particlesToTruth[particle] = truth;
        }

      return pandora::STATUS_CODE_SUCCESS;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    const simb::MCParticle* MCParticleCreator::GetFinalStateMCParticle(const MCParticleMap &particleMap, const simb::MCParticle *inputParticle)
    {
      // Navigate upward through MC daughter/parent links - collect this particle and all its parents
      MCParticleVector mcVector;

      int trackID(inputParticle->TrackId());

      while(1)
        {
          MCParticleMap::const_iterator pIter = particleMap.find(trackID);
          if (particleMap.end() == pIter)
            break; // Can't find MC Particle for this track ID [break]

          const art::Ptr<simb::MCParticle> particle = pIter->second;
          mcVector.push_back(particle);

          trackID = particle->Mother();
        }

      // Navigate downward through MC parent/daughter links - return the first long-lived charged particle
      for (MCParticleVector::const_reverse_iterator iter = mcVector.rbegin(), iterEnd = mcVector.rend(); iter != iterEnd; ++iter)
        {
          const art::Ptr<simb::MCParticle> nextParticle = *iter;

          if (MCParticleCreator::IsVisible(nextParticle))
            return nextParticle.get();
        }

      throw cet::exception("LArPandora"); // need to catch this exception
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    bool MCParticleCreator::IsVisible(const art::Ptr<simb::MCParticle> particle)
    {
      // Include long-lived charged particles
      const int pdg(particle->PdgCode());

      if ((pandora::E_MINUS == std::abs(pdg)) || (pandora::MU_MINUS == std::abs(pdg)) || (pandora::PROTON == std::abs(pdg)) ||
          (pandora::PI_PLUS == std::abs(pdg)) || (pandora::K_PLUS == std::abs(pdg)) ||
          (pandora::SIGMA_MINUS == std::abs(pdg)) || (pandora::SIGMA_PLUS == std::abs(pdg)) || (pandora::HYPERON_MINUS == std::abs(pdg)) ||
          (pandora::PHOTON == std::abs(pdg)) || (pandora::NEUTRON == std::abs(pdg)))
        return true;

      // TODO: What about ions, neutrons, photons? (Have included neutrons and photons for now)

      return false;
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    MCParticleCreator::Settings::Settings()
      : m_geantModuleLabel( "" ),
        m_generatorModuleLabel( "" )
    {
    }
  }
}
