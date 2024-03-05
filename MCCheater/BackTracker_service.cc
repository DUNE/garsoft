////////////////////////////////////////////////////////////////////////
//
//
// \file: BackTracker_service.cc
//
// brebel@fnal.gov
// modified Feb - Mar 2020 Leo Bellantoni, bellanto@fnal.gov
//
//
////////////////////////////////////////////////////////////////////////

// GArSoft includes
#include "MCCheater/BackTracker.h"
#include "Geometry/GeometryGAr.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "nug4/ParticleNavigation/EmEveIdCalculator.h"

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"

namespace gar {
  namespace cheat {

    //----------------------------------------------------------------------
    BackTracker::BackTracker(fhicl::ParameterSet const& pset,
                             art::ActivityRegistry    & reg)
      : BackTrackerCore(pset) {
      reg.sPreProcessEvent.watch(this, &BackTracker::Rebuild);
      fHasMC = fHasHits = fHasCalHits = fHasTracks = fHasClusters = false;
      fSTFU = 0;
    }

    //----------------------------------------------------------------------
    BackTracker::~BackTracker() {}

    //----------------------------------------------------------------------
    void BackTracker::beginJob() {}



    //----------------------------------------------------------------------
    void BackTracker::Rebuild(art::Event const& evt, art::ScheduleContext) {

      // do nothing if this is data
      if (evt.isRealData()) return;

      // do nothing if we are asked to do nothing
      if (fDisableRebuild) return;

      // We have a new event, so clear the collections from the previous events
      // Just the MC collections now, others as needed
      fParticleList.clear();
      fMCTruthList .clear();
      fTrackIDToMCTruthIndex.clear();





      // get the MCParticles from the event
      // getValidHandle will throw an exception if product not found.
      auto partCol = evt.getValidHandle<std::vector<simb::MCParticle> >(fG4ModuleLabel);
      if (partCol.failedToGet()) {
          ++fSTFU;
          if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::Rebuild")
                << "failed to get handle to simb::MCParticle from " << fG4ModuleLabel
                << ", return";
          }
          return;
      }
 
      // get mapping from MCParticles to MCTruths; loop to fill fTrackIDToMCTruthIndex
      art::FindOneP<simb::MCTruth> fo(partCol, evt, fG4ModuleLabel);
      if ( !fo.isValid() ) {
          ++fSTFU;
          if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::Rebuild")
                  << "failed to associate MCTruth to MCParticle with " << fG4ModuleLabel
                  << "; all attempts to backtrack will fail\n";
          }
          return;

      } else {
          for (size_t iPart = 0; iPart < partCol->size(); ++iPart) {

              simb::MCParticle *part = new simb::MCParticle(partCol->at(iPart));
              fParticleList.Add(part);

              // get the simb::MCTruth associated to this MCParticle
              try {
                  art::Ptr<simb::MCTruth> mct = fo.at(iPart);
                  if (fMCTruthList.size() < 1) {
                      fMCTruthList.push_back(mct);
                  } else {
                      // check that we are not adding a simb::MCTruth twice to
                      // the collection we know that all the particles for a
                      // given simb::MCTruth are put into the collection of
                      // particles at the same time, so we can just check that the
                      // current ::art::Ptr has a different id than the last one put
                      if (!(mct == fMCTruthList.back())) fMCTruthList.push_back(mct);
                  }

                 // Fill the track id to mctruth index map.  partCol will not contain
                 // negative track ids.
                 fTrackIDToMCTruthIndex[partCol->at(iPart).TrackId()] = fMCTruthList.size() - 1;
              } catch (cet::exception &ex) {
                  ++fSTFU;
                  if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                    MF_LOG_WARNING("BackTracker_service::Rebuild")
                      << "unable to find MCTruth from ParticleList created in "
                      << fG4ModuleLabel
                      << "; all attempts to backtrack will fail\n"
                      << "message from caught exception:\n" << ex;
                  }
              }
          } // end loop over MCParticles
      }

      // Register the electromagnetic Eve calculator
      fParticleList.AdoptEveIdCalculator(new sim::EmEveIdCalculator);

      fHasMC = true;





      // Inverse of drift velocity and longitudinal diffusion constant will
      // be needed by the backtracker when fHasHits, and also the geometry
      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      fInverseVelocity = 1.0/detProp->DriftVelocity(detProp->Efield(),
                                                    detProp->Temperature());

      auto garProp = gar::providerFrom<detinfo::GArPropertiesService>();
      fLongDiffConst = std::sqrt(2.0 * garProp->LongitudinalDiffusion());

      fGeo = gar::providerFrom<geo::GeometryGAr>();





      // Create the TPC channel to energy deposit collection.  Start by
      // looking for the RawDigit collection from the event and create
      // a FindMany mapping between the digits and energy deposits if it exists.
      for (auto vec : fChannelToEDepCol) vec.clear();
      fChannelToEDepCol.clear();
      auto digCol = evt.getHandle<std::vector<raw::RawDigit>>(fRawTPCDataLabel);
      if (!digCol) {
          ++fSTFU;
          if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::Rebuild")
                  << "Unable to find RawDigits in " << fRawTPCDataLabel <<
                  "; no backtracking in TPC will be possible" <<
                  "\n Well did you EXPECT RawDigits in this file?";
          }

      } else {
          art::FindMany<sdp::EnergyDeposit, float> fmEnergyDep(digCol, evt, fRawTPCDataLabel);
          if (!fmEnergyDep.isValid()) {
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
            	MF_LOG_WARNING("BackTracker_service::Rebuild")
            		<< "Unable to find valid association between RawDigits and "
            		<< "energy deposits in " << fRawTPCDataLabel <<
            		"; no backtracking in TPC will be possible" <<
            		"\n Well did you EXPECT those Assns in this file?";
              }

          } else {
              fChannelToEDepCol.resize(fGeo->NChannels());
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_DEBUG("BackTracker_service::Rebuild")
                	  << "There are " << fChannelToEDepCol.size()
                	  << " channels in the geometry";
              }

              std::vector<float const*> depweights;
              std::vector<const sdp::EnergyDeposit*> eDeps;
              for (size_t iDig = 0; iDig < digCol->size(); ++iDig) {
            	  eDeps.clear();
            	  depweights.clear();
            	  fmEnergyDep.get(iDig, eDeps, depweights);    // uses vector::insert
            	  if (eDeps.size() < 1) continue;

            	  for (size_t iDep=0; iDep<eDeps.size(); iDep++) {
            		  float w = 1.0;
            		  if (fSplitEDeps) w = *depweights[iDep];
            		  fChannelToEDepCol[ (*digCol)[iDig].Channel() ].emplace_back(std::make_pair(eDeps[iDep],w));

            		  MF_LOG_DEBUG("BackTracker_service::Rebuild") << "eDep\t" << iDep << " from RawDigit \t" << iDig
            			  << "\t TrkID, energy, dl, position:\t " << eDeps[iDep]->TrackID()
            			  << "\t " << eDeps[iDep]->Energy() << "\t " << eDeps[iDep]->dX()
            			  << "\t " << eDeps[iDep]->X() << "\t " << eDeps[iDep]->Y()
            			  << "\t " << eDeps[iDep]->Z() << std::endl;
 
            	  }
              }
              fHasHits = true;
          }
      }





      // Create the CellID to energy deposit collection for both ECAL & MuID.  Start by
      // looking for the CaloRawDigit collection from the event and create
      // a FindMany mapping between these digits and CaloDeposits if it exists.
      fOdetTrackToIdetTrack->clear();
      fCellIDToEDepCol.clear();
      art::InputTag* rawtag[2];
      rawtag[0] = new art::InputTag(fRawECALDataLabel, fRawECALDataInstance);
      rawtag[1] = new art::InputTag(fRawMuIDDataLabel, fRawMuIDDataInstance);

      for (int i=0; i<2; ++i) {
          auto caloDigCol = evt.getHandle<std::vector<raw::CaloRawDigit>>(*(rawtag[i]));
          if (!caloDigCol) {
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_WARNING("BackTracker_service::Rebuild")
                      << "Unable to find CaloRawDigits in " << 
                      (i==0 ? fRawECALDataLabel : fRawMuIDDataLabel) <<
                      "; no backtracking in ECAL/MuID will be possible" <<
                      "\n Well did you EXPECT CaloRawDigits in this file?";
              }

          } else {
              art::FindMany<sdp::CaloDeposit> fmCaloDep(caloDigCol, evt, *(rawtag[i]));
              if (!fmCaloDep.isValid()) {
                  ++fSTFU;
                  if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                      MF_LOG_WARNING("BackTracker_service::Rebuild")
                          << "Unable to find valid association between CaloRawDigits and "
                          << "calorimeter energy deposits in " << 
                          (i==0 ? fRawECALDataLabel : fRawMuIDDataLabel) <<
                          "; no backtracking in ECAL/MuID will be possible" <<
                          "\n Well did you EXPECT those Assns in this file?";
                  }

              } else {
                  std::vector<const sdp::CaloDeposit*> CeDeps;
                  for (size_t iCalDig = 0; iCalDig < caloDigCol->size(); ++iCalDig) {
                      CeDeps.clear();
                      fmCaloDep.get(iCalDig, CeDeps);    // uses vector::insert
                      if (CeDeps.size() < 1) continue;
                      fCellIDToEDepCol[ (*caloDigCol)[iCalDig].CellID() ] = CeDeps;
                  }
                  fHasCalHits = true;				  
              }
          }
      }



      // Create an unordered map of IDNumbers of reco'd Tracks to the hits in those tracks
      // Just to keep the code comprehensible, first map Tracks to TPCClusters, then
      // TPCClusters to Hits, then build the combined map.
      fTrackIDToTPCClusters.clear();
      auto recoTrackCol = evt.getHandle<std::vector<rec::Track>>(fTrackLabel);
      if (!recoTrackCol) {
          ++fSTFU;
          if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::Rebuild")
                  << "Unable to find rec::Tracks in " << fTrackLabel <<
                  "; no backtracking of reconstructed tracks will be possible" <<
                  "\n Well did you EXPECT tracks in this file?";
          }

      } else {
          art::FindMany<rec::TPCCluster> fmTPCClusts(recoTrackCol, evt, fTrackLabel);
          if (!fmTPCClusts.isValid()) {
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_WARNING("BackTracker_service::Rebuild")
                    << "Unable to find valid association between Tracks and "
                    << "TPCClusters in " << fTrackLabel <<
                    "; no backtracking of reconstructed tracks will be possible" <<
                    "\n Well did you EXPECT those Assns in this file?";
              }

          } else {
              std::vector<const rec::TPCCluster*> tpclusThisTrack;
              for (size_t iTrack = 0; iTrack < recoTrackCol->size(); ++iTrack) {
                  tpclusThisTrack.clear();
                  fmTPCClusts.get(iTrack, tpclusThisTrack);    // uses vector::insert
                  if (tpclusThisTrack.size() < 1) continue;
                  fTrackIDToTPCClusters[ (*recoTrackCol)[iTrack].getIDNumber() ] = tpclusThisTrack;
              }
          }
      }

      // Construct map of TPCCluster to Hits
      fTPCClusterIDToHits.clear();
      auto tpclusCol = evt.getHandle<std::vector<rec::TPCCluster>>(fTPCClusterLabel);
      if (!tpclusCol) {
          ++fSTFU;
          if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
              MF_LOG_WARNING("BackTracker_service::Rebuild")
                  << "Unable to find rec::TPCClusters in " << fTPCClusterLabel <<
                  "; no backtracking of reconstructed tracks will be possible" <<
                  "\n Well did you EXPECT TPCClusters in this file?";
          }

      } else {
          art::FindManyP<rec::Hit> fmHits(tpclusCol, evt, fTPCClusterLabel);
          if (!fmHits.isValid()) {
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_WARNING("BackTracker_service::Rebuild")
                      << "Unable to find valid association between TPCClusters and "
                      << "Hits in " << fTPCClusterLabel <<
                      "; no backtracking of reconstructed tracks will be possible" <<
                      "\n Well did you EXPECT those Assns in this file?";
              }

          } else {
              std::vector<art::Ptr<rec::Hit>> hitsThisTPCCluster;
              for (size_t iTPClus = 0; iTPClus < tpclusCol->size(); ++iTPClus) {
                  hitsThisTPCCluster.clear();
                  fmHits.get(iTPClus, hitsThisTPCCluster);    // uses vector::insert
                  if (hitsThisTPCCluster.size() < 1) continue;
                  fTPCClusterIDToHits[(*tpclusCol)[iTPClus].getIDNumber()] = hitsThisTPCCluster;
              }
          }
      }

      // Now compound the lists to make fTrackIDToHits
      fTrackIDToHits.clear();
      if ( !fTrackIDToTPCClusters.empty() && !fTPCClusterIDToHits.empty()) {
          for (auto aTrack : *recoTrackCol) {
              std::vector<const rec::TPCCluster*> tpclusThisTrack =
                  fTrackIDToTPCClusters[ aTrack.getIDNumber() ];
              std::vector<art::Ptr<rec::Hit>> hitsThisTrack;
              for (const rec::TPCCluster* aTPClus : tpclusThisTrack) {
                  std::vector<art::Ptr<rec::Hit>> hitsThisClus =
                      fTPCClusterIDToHits[aTPClus->getIDNumber()];
                  hitsThisTrack.insert( hitsThisTrack.end(),
                                        hitsThisClus.begin(),hitsThisClus.end() );
              }
              // Explicit assumption is that each hit is in at most 1 cluster so
              // there are no duplicates in hitsThisTrack
              fTrackIDToHits[aTrack.getIDNumber()] = hitsThisTrack;
          }
      }
      fHasTracks = true;





      // Create an unordered map of IDNumbers of reco'd Clusters to the CaloHits in those tracks
      art::InputTag* clustertag[2];
      clustertag[0] = new art::InputTag(fECALClusterLabel, fClusterECALInstance);
      clustertag[1] = new art::InputTag(fMuIDClusterLabel, fClusterMuIDInstance);

      for (int i=0; i<2; ++i) {
           auto recoClusterCol = evt.getHandle<std::vector<rec::Cluster>>(*(clustertag[i]));
          if (!recoClusterCol) {
              ++fSTFU;
              if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                  MF_LOG_WARNING("BackTracker_service::Rebuild")
                      << "Unable to find rec::Clusters in " << 
                      (i==0 ? fECALClusterLabel : fMuIDClusterLabel ) << ", " <<
                      (i==0 ? fClusterECALInstance : fClusterMuIDInstance ) << 
                      " instance; no backtracking of ECAL/MuID clusters will be possible" <<
                      "\n Well did you EXPECT ECAL/MuID Clusters in this file?";
              }

          } else {
              art::FindManyP<rec::CaloHit> fmCaloHits(recoClusterCol, evt, *(clustertag[i]));
              if (!fmCaloHits.isValid()) {
                  ++fSTFU;
                  if (fSTFU<=12) {    // Ye who comprehend messagelogger, doeth ye better.
                      MF_LOG_WARNING("BackTracker_service::Rebuild")
                          << "Unable to find valid association between Clusters and "
                          << "CaloHits in " << fECALClusterLabel << ", " << fClusterECALInstance 
                          << " instance; no backtracking of reconstructed tracks will be possible" <<
                          "\n Well did you EXPECT those Assns in this file?";
                  }

              } else {
                  std::vector<art::Ptr<rec::CaloHit>> caloHitsThisCluster;
                  for (size_t iCluster = 0; iCluster < recoClusterCol->size(); ++iCluster) {
                      caloHitsThisCluster.clear();
                      fmCaloHits.get(iCluster, caloHitsThisCluster);    // uses vector::insert
                      if (caloHitsThisCluster.size() < 1) continue;
                      fClusterIDToCaloHits[ (*recoClusterCol)[iCluster].getIDNumber() ] = caloHitsThisCluster;
                  }
                  fHasClusters = true;
              }
          }
      }
      return;
    }





    DEFINE_ART_SERVICE(BackTracker)
  } // namespace cheat
} // namespace gar
