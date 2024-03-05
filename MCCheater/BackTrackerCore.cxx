//
//  BackTrackerCore.cxx
//  Created by Brian Rebel on 3/13/17.
//  Major rewrite, Leo Bellantoni Mar 2020
//

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nug4/ParticleNavigation/EmEveIdCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "SimulationDataProducts/sim.h"
#include "CoreUtils/ServiceUtil.h"
#include "MCCheater/BackTrackerCore.h"
#include "Geometry/GeometryGAr.h"

#include "TMath.h"



namespace gar{
    namespace cheat{

        //--------------------------------------------------------------------------
        BackTrackerCore::BackTrackerCore(fhicl::ParameterSet const& pset)
          : fClocks(nullptr), fGeo(nullptr) {

            fClocks = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();
            fGeo    = gar::providerFrom<geo::GeometryGAr>();

            fDisableRebuild           = pset.get<bool       >("DisableRebuild",       false);

            fG4ModuleLabel            = pset.get<std::string>("G4ModuleLabel",       "geant");

            fRawTPCDataLabel          = pset.get<std::string>("RawTPCDataLabel",     "daq");

            fRawECALDataLabel         = pset.get<std::string>("RawECALDataLabel",    "daqsipm");
            fRawECALDataInstance      = pset.get<std::string>("RawECALDataInstance", "ECAL"  );
            fRawMuIDDataLabel         = pset.get<std::string>("RawMuIDDataLabel",    "daqsipmmuid");
            fRawMuIDDataInstance      = pset.get<std::string>("RawMuIDDataInstance", "MuID"  );
            fECALtimeResolution       = pset.get<double     >("ECALtimeResolution",   1.0);
            fMinHitEnergyFraction     = pset.get<double     >("MinHitEnergyFraction", 0.1);
            fMinCaloHitEnergyFrac     = pset.get<double     >("MinCaloHitEnergyFrac", 0.1);

            fTrackLabel               = pset.get<std::string>("TrackLabel",          "track");
            fTPCClusterLabel          = pset.get<std::string>("TPCClusterLabel",     "tpccluster");
            fTrackFracMCP             = pset.get<double     >("TrackFracMCP",         0.8);

            fECALClusterLabel         = pset.get<std::string>("ECALClusterLabel",    "calocluster");
            fClusterECALInstance      = pset.get<std::string>("ClusterECALInstance", "ECAL"  );
            fMuIDClusterLabel         = pset.get<std::string>("MuIDClusterLabel",    "caloclustermuid");
            fClusterMuIDInstance      = pset.get<std::string>("ClusterMuIDInstance", "MuID"  );
            fClusterFracMCP           = pset.get<double     >("ClusterFracMCP",       0.8);

            fSplitEDeps               = pset.get<bool       >("SplitEDeps",           true);
            fDontChaseNeutrons        = pset.get<bool       >("DontChaseNeutrons",    true);

            fOdetTrackToIdetTrack      = new std::unordered_map<int, int>;
            return;
        }

        //--------------------------------------------------------------------------
        BackTrackerCore::~BackTrackerCore() {}





        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::TrackIDToParticle(int const& id) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::TrackIDToParticle")
                    << "Attempting to backtrack without MC truth information";
            }

            // What to do for negative trackID?  ParticleList::find(key) does
            // std::map<int,simb::MCParticle*>.find( abs(key) ) so a track ID
            // corresponding to a radiated photon for example gets the original
            // parent of that photon.  In most parts of GArSoft, trackIDs are
            // non-negative; the exception is in the EnergyDeposit class.
            auto part_it = fParticleList.find(id);

            if (part_it == fParticleList.end()) {
                MF_LOG_WARNING("BackTrackerCore::TrackIDToParticle")
                    << "can't find particle with track id " << id
                    << " in sim::ParticleList returning null pointer";
                    return nullptr;
            }

            return part_it->second;
        }



        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::FindEve(simb::MCParticle* const p) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindEve")
                    << "Attempting to backtrack without MC truth information";
            }

            int EveTrackId = fParticleList.EveId( p->TrackId() );
            return TrackIDToParticle(EveTrackId);
        }



        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::FindTPCEve(simb::MCParticle* const p) const {

            return FindTPCEve(p->TrackId());
        }



        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::FindTPCEve(int const trackID) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindTPCEve")
                    << "Attempting to backtrack without MC truth information";
            }

            int const neutronPDG = 2112;

            // If I had ever been here before I would probably know just what to do / Don't you?
            if ( fOdetTrackToIdetTrack->find(trackID) != fOdetTrackToIdetTrack->end() ) {
                return TrackIDToParticle( fOdetTrackToIdetTrack->at(trackID) );
            }

            // Negative track ID, as found in EnergyDeposit or CaloEnergyDeposit,
            // gets turned into positve track ID corresponding to parent of EM
            // that created that particular EnergyDeposit or CaloEnergyDeposit
            simb::MCParticle* walker  = TrackIDToParticle(trackID);

            if (walker==nullptr) return nullptr;
            while ( walker != nullptr ) {
                // Stopping for neutrons?
                if (walker->PdgCode()==neutronPDG && fDontChaseNeutrons) break;

                // Stop in TPC
                TVector3 verty(walker->Vx(),walker->Vy(),walker->Vz());
                if (fGeo->PointInGArTPC(verty)) break;

                // Walk up the ParticleList not the event store so someday somebody can
                // change the ParticleList and this will still work.
                 int mommaTID = fParticleList[walker->TrackId()]->Mother();
                 // Stop at top of tree!
                if (mommaTID == 0) break;

                walker = TrackIDToParticle(mommaTID);
            }
            (*fOdetTrackToIdetTrack)[trackID] = walker->TrackId();
            return walker;
         }



        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::FindECALEve(simb::MCParticle* const p) const {

            return FindECALEve(p->TrackId());
        }



        //----------------------------------------------------------------------
        simb::MCParticle*  BackTrackerCore::FindECALEve(int const trackID) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::FindECALEve")
                    << "Attempting to backtrack without MC truth information";
            }

            int const neutronPDG = 2112;

            // If I had ever been here before I would probably know just what to do / Don't you?
            if ( fOdetTrackToIdetTrack->find(trackID) != fOdetTrackToIdetTrack->end() ) {
                return TrackIDToParticle( fOdetTrackToIdetTrack->at(trackID) );
            }

            // Negative track ID, as found in EnergyDeposit or CaloEnergyDeposit,
            // gets turned into positve track ID corresponding to parent of EM
            // that created that particular EnergyDeposit or CaloEnergyDeposit
            simb::MCParticle* walker  = TrackIDToParticle(trackID);

            if (walker==nullptr) return nullptr;
            while ( walker != nullptr ) {

                // Stopping for neutrons?
                if (walker->PdgCode()==neutronPDG && fDontChaseNeutrons) break;

                // Stop in ECAL or TPC
                TVector3 verty(walker->Vx(),walker->Vy(),walker->Vz());
                if (fGeo->PointInGArTPC(verty)) break;
                if (fGeo->PointInECALBarrel(verty)) break;
                if (fGeo->PointInECALEndcap(verty)) break;

                // Walk up the ParticleList not the event store so someday somebody can
                // change the ParticleList and this will still work.
                int mommaTID = fParticleList[walker->TrackId()]->Mother();
                 // Stop at top of tree!
                if (mommaTID == 0) break;

                walker = TrackIDToParticle(mommaTID);
            }
            (*fOdetTrackToIdetTrack)[trackID] = walker->TrackId();
            return walker;
         }



        //----------------------------------------------------------------------
        bool BackTrackerCore::IsForebearOf(simb::MCParticle* const forebear,
                                           simb::MCParticle* const afterbear) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::IsForebearOf")
                    << "Attempting to backtrack without MC truth information";
            }

            if (forebear==nullptr || afterbear==nullptr) return false;

            int stopper = forebear->TrackId();
            int walker  = afterbear->TrackId();
            while ( walker!=stopper ) {
                // Walk up the ParticleList not the event store so someday somebody can
                // change the ParticleList and this will still work.
                int momma = fParticleList[walker]->Mother();
                if (momma == 0 || momma == -1) {        // -1 convention for edep-sim
                    return false;
                } else {
                    walker = momma;
                }
            }
            return true;
         }



        //----------------------------------------------------------------------
        art::Ptr<simb::MCTruth> const
        BackTrackerCore::ParticleToMCTruth(simb::MCParticle* const p) const {
            if (!fHasMC) {
                throw cet::exception("BackTrackerCore::ParticleToMCTruth")
                    << "Attempting to backtrack without MC truth information";
            }

            int id = p->TrackId();
            // find the entry in the MCTruth collection for this track id.
            // The std::map fTrackIDToMCTrutheIndex will not contain keys < 0
            size_t mct = fTrackIDToMCTruthIndex.find(std::abs(id))->second;

            if( mct > fMCTruthList.size() ) {
                throw cet::exception("BackTrackerCore::ParticleToMCTruth")
                    << "attempting to find MCTruth index for out-of-range value: "
                    << mct << "/" << fMCTruthList.size();
            }
            return fMCTruthList[mct];
        }

        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::HitToHitIDEs(art::Ptr<rec::Hit> const& hit) const {
            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitToHitIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            // By cutting on the StartTime and EndTime for this hit, we get only the
            // IDEs for this hit; ides.size() should basically be the number of
            // MCParticles that contributed to this hit.
            std::vector<HitIDE> ides;
            ides = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());
            return ides;
        }

        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::HitToHitIDEs(rec::Hit const& hit) const {
            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitToHitIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            // By cutting on the StartTime and EndTime for this hit, we get only the
            // IDEs for this hit; ides.size() should basically be the number of
            // MCParticles that contributed to this hit.
            std::vector<HitIDE> ides;
            ides = this->ChannelToHitIDEs(hit.Channel(),hit.StartTime(),hit.EndTime());
            return ides;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::Hit>>
        BackTrackerCore::ParticleToHits(simb::MCParticle* const p,
                                        std::vector<art::Ptr<rec::Hit>> const& allhits,
                                        bool checkNeutrals) const {
            // returns a subset of the hits in the allhits collection that match
            // the MC particle passed as an argument

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::ParticleToHits")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            std::vector<art::Ptr<rec::Hit>> hitList;
            int PDGpid = p->PdgCode();
            bool obviousNeutral = PDGpid==22 || PDGpid==2112 || PDGpid==111 ||
                                  abs(PDGpid)==14 || abs(PDGpid)==12;
            if (!checkNeutrals && obviousNeutral) return hitList;

            int tkID = p->TrackId();
            std::vector<HitIDE> hids;

            for (auto hit : allhits) {
                hids.clear();

                hids = this->ChannelToHitIDEs(hit->Channel(),hit->StartTime(),hit->EndTime());

                for (auto const& hid : hids) {
                    if ( hid.trackID==tkID && hid.energyFrac>fMinHitEnergyFraction ) {
                        hitList.push_back(hit);
                    }
                }
            }

            return hitList;
        }



        //----------------------------------------------------------------------
        std::vector<HitIDE>
        BackTrackerCore::ChannelToHitIDEs(raw::Channel_t const& channel,
                                          double const start, double const stop) const {

            // loop over the energy deposits in the channel and grab those in time window

            if(start<0 || stop <0) 
                MF_LOG_WARNING("BackTrackerCore::ChannelToHitIDEs") << "start and/or stop times are negative!";
            std::vector<HitIDE> hitIDEs;

            if(fChannelToEDepCol.size() < channel){
                MF_LOG_WARNING("BackTrackerCore::ChannelToHitIDEs")
                    << "Attempting to find energy deposits for channel " << channel
                    << " while there are only " << fChannelToEDepCol.size()
                    << " channels available, returning an empty vector.";
                return hitIDEs;
            }

            auto chanEDeps = fChannelToEDepCol[channel];

            if(chanEDeps.size() < 1){
                MF_LOG_WARNING("BackTrackerCore::ChannelToHitIDEs")
                    << "No sdp::EnergyDeposits for selected channel: " << channel
                    << ". There is no way to backtrack the given hit, returning"
                    << " an empty vector.";
                return hitIDEs;
            }

            std::map<int,size_t> tidmap;  // Maps track IDs to index into hitIDEs

            // first get the total energy represented by all track ids for
            // this channel and range of tdc values
            unsigned short tdc = 0;
            float chanpos[3]={0.,0.,0};
            fGeo->ChannelToPosition(channel, chanpos);

            // loop over all EnergyDeposits for this channel and
            // calculate a weighted average, in contributed energy, 
            // 4-position and 4-momentum
            double totalEallTracks = 0.0;
            size_t itraj=0;
            for (auto const& edepfrac : chanEDeps) {

                auto edep = edepfrac.first;
                float weight = edepfrac.second;
                if(weight<0 || weight>1) 
                    MF_LOG_WARNING("BackTrackerCore::ChannelToHitIDEs") << "bad weight (" << weight << ")";
                if (edep->TrackID() == sdp::NoParticleId) continue;

                tdc = fClocks->TPCG4Time2TDC
                        (edep->Time() +TMath::Abs(edep->X()-chanpos[0]) *fInverseVelocity);
                if ( (unsigned int)start <= tdc && tdc <= (unsigned int)stop) {
                    float energy = edep->Energy() * weight;
                    MF_LOG_DEBUG("BackTrackerCore::ChannelToHitIDEs") << "edep energy: " << edep->Energy() 
                                 << ", weight: " << weight << ", contributed energy: " << energy;
                    totalEallTracks += energy;

                    TLorentzVector simpos(edep->X(),edep->Y(),edep->Z(),edep->Time());
                    TLorentzVector simmom = EnergyDepositToMomentum(edep->TrackID(),simpos,itraj);
                    simpos *= energy;
                    simmom *= energy;


                    // EnergyDeposits can contain negative track IDs corresponding to tracks
                    // created in showers.  The absolute magnitude of the track is the track
                    // id of the parent of the EM showering process.  That is the particle
                    // to which we want to assign this energy for backtracking purposes.
                    int tid = abs( edep->TrackID() );
                    auto tidmapiter = tidmap.find(tid);
                    if (tidmapiter == tidmap.end()) {
                        hitIDEs.emplace_back(tid, 0.0, energy, simpos, simmom);
                        tidmap[tid] = hitIDEs.size()-1;
                    } else {
                        hitIDEs.at(tidmapiter->second).energyTot += energy;
                        hitIDEs.at(tidmapiter->second).position += simpos;
                        hitIDEs.at(tidmapiter->second).momentum += simmom;
                    }
                }
            }

            // loop over the hitIDEs to set the fractional energy for each TrackID
            for (size_t i = 0; i < hitIDEs.size(); ++i) {
                hitIDEs[i].energyFrac = hitIDEs[i].energyTot / totalEallTracks;
                hitIDEs[i].position *= 1.0/totalEallTracks;
                hitIDEs[i].momentum *= 1.0/totalEallTracks;
            }

            // remove duplicate HitIDEs (not that any should be there)
            size_t nbefore = hitIDEs.size();
            std::sort(hitIDEs.begin(),hitIDEs.end());
            hitIDEs.erase(std::unique(hitIDEs.begin(),hitIDEs.end()),hitIDEs.end());
            size_t nafter = hitIDEs.size();
            if(nbefore!=nafter) 
              MF_LOG_WARNING("BackTrackerCore::ChannelToHitIDEs") << "found and removed " 
                             << nbefore-nafter << " duplicate HitIDEs (not expected for single channel)";

            return hitIDEs;
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::HitPurity(simb::MCParticle* const p,
            std::vector<art::Ptr<rec::Hit> > const& hits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitCollectionPurity")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }


            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            float total   = 0.0;
            for (auto hit : hits) {
                std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                weight = (weightByCharge) ? hit->Signal() : 1.0;
                total += weight;

                for (auto iHitIDE : hitIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT trackIDs
                    if (iHitIDE.trackID    == trackIDin &&
                        iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            double purity      = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                purity      = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(purity,uncertainty);
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::HitEfficiency(simb::MCParticle* const p,
            std::vector<art::Ptr<rec::Hit> > const& hits,
            std::vector<art::Ptr<rec::Hit> > const& allhits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::HitCollectionEfficiency")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            for (auto hit : hits) {
                std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                weight = (weightByCharge) ? hit->Signal() : 1.0;

                for (auto iHitIDE : hitIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT track IDs
                    if (iHitIDE.trackID    == trackIDin &&
                        iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            // Next figure out how many hits in the whole collection are associated with this id
            float total   = 0.0;
            if (desired>0) {
                for (auto hit : allhits) {
                    std::vector<HitIDE> hitIDEs = this->HitToHitIDEs(hit);

                    weight = (weightByCharge) ? hit->Signal() : 1.0;

                    for (auto iHitIDE : hitIDEs) {
                        // Don't double count if this hit has > 1 of the desired GEANT track IDs
                        if (iHitIDE.trackID   == trackIDin &&
                           iHitIDE.energyFrac >= fMinHitEnergyFraction) {
                           total += weight;
                           break;
                        }
                    }
                } // end loop over all hits
            }

            double efficiency  = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                efficiency  = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(efficiency,uncertainty);
        }





        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CaloHitToCalIDEs(art::Ptr<rec::CaloHit> const& hit) const {
            if (!fHasMC || !fHasCalHits) {
                throw cet::exception("BackTrackerCore::CaloHitToCalIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            // Cut on Time for this hit +/- ECAL time resolution given in BackTracker.fcl
            std::vector<CalIDE> ides;
            ides = this->CellIDToCalIDEs(*hit);
            return ides;
        }

        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CaloHitToCalIDEs(rec::CaloHit const& hit) const {
            if (!fHasMC || !fHasCalHits) {
                throw cet::exception("BackTrackerCore::CaloHitToCalIDEs")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            // Cut on Time for this hit +/- ECALtimeResolution given in BackTracker.fcl
            std::vector<CalIDE> ides;
            ides = this->CellIDToCalIDEs(hit);
            return ides;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::CaloHit>>
        BackTrackerCore::ParticleToCaloHits(simb::MCParticle* const p,
                                            std::vector<art::Ptr<rec::CaloHit>> const& allhits) const {
            // returns a subset of the CaloHits in the allhits collection that match
            // the MC particle passed as an argument

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::ParticleToCaloHits")
                    << "Attempting to backtrack without MC truth or gar::rec::CaloHit information";
            }

            std::vector<art::Ptr<rec::CaloHit>> calHitList;
            int tkID = p->TrackId();
            std::vector<CalIDE> calhids;

            for (auto hit : allhits) {
                calhids.clear();

                calhids = this->CellIDToCalIDEs(*hit);

                for (auto const& hid : calhids) {
                    if ( hid.trackID==tkID && hid.energyFrac>fMinCaloHitEnergyFrac ) {
                        calHitList.push_back(hit);
                   }
                }
            }
            return calHitList;
        }



        //----------------------------------------------------------------------
        std::vector<CalIDE>
        BackTrackerCore::CellIDToCalIDEs(rec::CaloHit hit) const {

            // loop over the energy deposits in the channel and grab those in time window

            raw::CellID_t hitCell = hit.CellID();
            std::pair<float, float> hitTime = hit.Time();   // need time corrected for strip length
            float time;
            if (hitTime.second==0) {
                time = hitTime.first;
            } else {
                std::array<double, 3> point = {hit.Position()[0], hit.Position()[1], hit.Position()[2]};
                double stripLength = fGeo->getStripLength(point, hit.CellID()); // in cm
                float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
                TVector3 whereBe(hit.Position()[0], hit.Position()[1], hit.Position()[2]);    // C++ & ROOT at their finest
                if ( fGeo->PointInECALBarrel(whereBe) || fGeo->PointInECALEndcap(whereBe) ) {
                    c /= fGeo->getIofRECAL();
                } else {
                    c /= fGeo->getIofRMuID();
                }
                time = (hitTime.first + hitTime.second) / 2. - (stripLength / (2 * c));
            }

            std::vector<CalIDE> calIDEs;
            std::vector<const sdp::CaloDeposit*>  cellEDeps;
            try {
                cellEDeps = fCellIDToEDepCol.at(hitCell);
            }
            catch (std::out_of_range &) {
                MF_LOG_WARNING("BackTrackerCore::CellIDToCalIDEs")
                    << "No sdp::EnergyDeposits for selected ECAL cell: " << hitCell
                    << ". There is no way to backtrack the given calorimeter hit,"
                    << " returning an empty vector.";
                return calIDEs;
            }

            std::map<int,size_t> tidmap;  // Maps track IDs to index into calIDEs

            // Get the total energy produced by all track ids for this cell & time
            float start = time -fECALtimeResolution;
            float stop  = time +fECALtimeResolution;

            TVector3 verty(hit.Position());
            bool inYoke = ( fGeo->PointInMuIDBarrel(verty) ||
                            fGeo->PointInMuIDEndcap(verty) );

            double totalEallTracks = 0.0;
            for (auto const& edep : cellEDeps) {
                if (edep->TrackID() == sdp::NoParticleId) continue;

                if ( start<=edep->Time() && edep->Time()<=stop) {
                    float energy = edep->Energy();
                    totalEallTracks += energy;
                    int tid = edep->TrackID();
                    // Chase the tid up to the parent which started in the gas.
                    // tid < 0 is possible here, but FindTPCEve calls
                    // TrackIDToParticle, which fixes that little issue.
                    simb::MCParticle* DETEve;
                    DETEve = inYoke ? FindECALEve(tid) : FindTPCEve(tid);
                    int tidDET;
                    if (DETEve!=nullptr) {
                        tidDET = DETEve->TrackId();
                    } else {
                        tidDET = tid;
                    }					
					
                    auto tidmapiter = tidmap.find(tidDET);
                    if (tidmapiter == tidmap.end()) {
                        calIDEs.emplace_back(tidDET, 0.0, energy);
                        tidmap[tidDET] = calIDEs.size()-1;
                    } else {
                        calIDEs.at(tidmapiter->second).energyTot += energy;
                    }
                }
            }

            // loop over the hitIDEs to set the fractional energy for each TrackID
            for (size_t i = 0; i < calIDEs.size(); ++i) {
                calIDEs[i].energyFrac = calIDEs[i].energyTot / totalEallTracks;
            }
            return calIDEs;
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::CaloHitPurity(simb::MCParticle* const p,
            std::vector<art::Ptr<rec::CaloHit> > const& hits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::CaloHitPurity")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }


            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            float total   = 0.0;
            for (auto hit : hits) {
                std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                weight = (weightByCharge) ? hit->Energy() : 1.0;
                total += weight;

                for (auto iCalIDE : calIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT trackIDs
                    if (iCalIDE.trackID    == trackIDin &&
                        iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            double purity      = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                purity      = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(purity,uncertainty);
        }



        //----------------------------------------------------------------------
        std::pair<double,double>
        BackTrackerCore::CaloHitEfficiency(simb::MCParticle* const p,
            std::vector<art::Ptr<rec::CaloHit> > const& hits,
            std::vector<art::Ptr<rec::CaloHit> > const& allhits, bool weightByCharge) const {

            if (!fHasMC || !fHasHits) {
                throw cet::exception("BackTrackerCore::CaloHitEfficiency")
                    << "Attempting to backtrack without MC truth or gar::rec::Hit information";
            }

            int trackIDin = p->TrackId();

            float weight;
            float desired = 0.0;
            for (auto hit : hits) {
                std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                weight = (weightByCharge) ? hit->Energy() : 1.0;

                for (auto iCalIDE : calIDEs) {
                    // Don't double count if this hit has > 1 of the desired GEANT track IDs
                    if (iCalIDE.trackID    == trackIDin &&
                        iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                        desired += weight;
                        break;
                    }
                }
            } // end loop over hits

            // Next figure out how many hits in the whole collection are associated with this id
            float total   = 0.0;
            if (desired>0.0) {
                for (auto hit : allhits) {
                    std::vector<CalIDE> calIDEs = this->CaloHitToCalIDEs(hit);

                    weight = (weightByCharge) ? hit->Energy() : 1.0;

                    for (auto iCalIDE : calIDEs) {
                        // Don't double count if this hit has > 1 of the desired GEANT track IDs
                        if (iCalIDE.trackID    == trackIDin &&
                            iCalIDE.energyFrac >= fMinCaloHitEnergyFrac) {
                            total += weight;
                             break;
                        }
                    }
                 } // end loop over all hits
            }

            double efficiency  = 0.0;
            double uncertainty = 0.0;
            if (total > 0.0) {
                efficiency  = desired/total;
                uncertainty = std::sqrt( desired*(total -desired)/total )/total;
            }
            return std::make_pair(efficiency,uncertainty);
        }





        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::Hit>> const
        BackTrackerCore::TrackToHits(rec::Track* const t) {

            MF_LOG_DEBUG("BackTrackerCore::TrackToHits") << "BT:TrackToHits called, read map with " << 
                fTrackIDToHits.size() << " entries" << std::endl;

            // Evidently, use of std::unordered_map keeps this method from being const
            std::vector<art::Ptr<gar::rec::Hit>> retval;
            if ( !fTrackIDToHits.empty() ) {
                retval = fTrackIDToHits[ t->getIDNumber() ];
            }
            return retval;
        }
        //-----------------------------------------------------------------------
        std::vector<art::Ptr<rec::Hit>> const 
        BackTrackerCore::TPCClusterToHits(rec::TPCCluster* const clust) {

            // Evidently, use of std::unordered_map keeps this method from being const
            std::vector<art::Ptr<gar::rec::Hit>> retval;
            if ( !fTPCClusterIDToHits.empty() ) {
                retval = fTPCClusterIDToHits[ clust->getIDNumber() ];
            }
            return retval;
        }

        //-----------------------------------------------------------------------
        //std::vector<art::Ptr<rec::TPCCluster>> const 
        std::vector<const rec::TPCCluster*> const
        BackTrackerCore::TrackToTPCClusters(rec::Track* const t) {

            // Evidently, use of std::unordered_map keeps this method from being const
            //std::vector<art::Ptr<gar::rec::TPCCluster>> retval;
            std::vector<const gar::rec::TPCCluster*> retval;
            if ( !fTrackIDToTPCClusters.empty() ) {
                retval = fTrackIDToTPCClusters[ t->getIDNumber() ];
            }
            return retval;     
        }

        //-----------------------------------------------------------------------
        double BackTrackerCore::TrackToTotalEnergy(rec::Track* const t) {

            double etot = 0.; // total energy
            std::vector<art::Ptr<rec::Hit>> const hits = TrackToHits(t); // all of this track's hits
            std::vector<HitIDE> allIdes;

            // loop over hits, get IDEs, and sum up energy deposits
            for(auto const& hit : hits) {
                std::vector<HitIDE> ides = HitToHitIDEs(hit);
                allIdes.insert(allIdes.end(),ides.begin(),ides.end());

            }// for hits

            /*size_t pretotal = allIdes.size();

            // remove duplicate HitIdes
            // comment out for now since we could have meaningful duplicate HitIDEs 
            // e.g. a common sdp::EnergyDeposit split evenly among channels
            std::sort(allIdes.begin(),allIdes.end());
            allIdes.erase( std::unique(allIdes.begin(),allIdes.end()), allIdes.end());
            size_t posttotal = allIdes.size();

            if(pretotal!=posttotal) std::cout << "removed " << pretotal-posttotal << " duplicated IDEs" << std::endl;*/

            for(auto const& ide : allIdes) {
                etot += ide.energyTot;
            }// for ides

            return etot;
        }



        //----------------------------------------------------------------------
        std::vector<std::pair<simb::MCParticle*,float>>
        BackTrackerCore::TrackToMCParticles(rec::Track* const t) {
            // Can't be const because it uses TrackToHits

            std::vector<art::Ptr<gar::rec::Hit>> const hitCol = TrackToHits(t);
            std::vector<HitIDE> hitIDECol;
            for (auto hit : hitCol) {
                std::vector<HitIDE> IDEsThisHit = HitToHitIDEs(hit);
                hitIDECol.insert(hitIDECol.end(), IDEsThisHit.begin(),IDEsThisHit.end());
            }

            // A tree to accumulate the found energy for each GEANT/GENIE track
            std::map<int,float> lTrackIdToEnergy;
            for (auto hitIDE : hitIDECol) {
                lTrackIdToEnergy[hitIDE.trackID] += hitIDE.energyFrac * hitIDE.energyTot;
            }

            // There are still electrons in lTrackIdToEnergy that have parents in
            // lTrackIDToEnergy... these should be consolidated.  Also find total Energy
            std::map<int,float>::iterator iTrackIdTo = lTrackIdToEnergy.begin();
            for (; iTrackIdTo != lTrackIdToEnergy.end(); ++iTrackIdTo) {
                simb::MCParticle* iPart = TrackIDToParticle(iTrackIdTo->first);
                if ( iPart->PdgCode() == 11 && iPart->Mother()!=0 ) {
                    if ( lTrackIdToEnergy.count(iPart->Mother())==1 ) {
                        lTrackIdToEnergy[iPart->Mother()] += iTrackIdTo->second;
                        iTrackIdTo = lTrackIdToEnergy.erase(iTrackIdTo);
                        if (iTrackIdTo != lTrackIdToEnergy.begin()) --iTrackIdTo;
                    }
                }
            }

            // Next, sort the map by value.  If there are two entries with the exact same
            // ionization, only one is found in resulting set!  But the digitization scale
            // for edep is eVs (single molecular ionizations) i.e. 1e-5 or less so this should
            // be negligible for tracks of interesting ionizations.
            // Declaring the type of the sorting predicate
            typedef std::function<bool(std::pair<int,float>, std::pair<int,float>)> Comparator;
            // Declaring a set that will store the pairs using above comparison logic
            std::set<std::pair<int,float>, Comparator> setOfTracks(
                lTrackIdToEnergy.begin(), lTrackIdToEnergy.end(),
                [](std::pair<int,float> a ,std::pair<int,float> b) {
                    return a.second > b.second;
                }
            );

            // Make the returned vector
            float Etot = 0.0;
            auto iTrackSet = setOfTracks.begin();
            for (; iTrackSet!=setOfTracks.end(); ++iTrackSet) {
                Etot += iTrackSet->second;   // Can't do this in loop on lTrackIdToEnergy
            }

            std::pair<simb::MCParticle*,float> emptee(nullptr,0.0);
            std::vector<std::pair<simb::MCParticle*,float>> retval(setOfTracks.size(),emptee);
            std::vector<std::pair<simb::MCParticle*,float>>::iterator iRetval = retval.begin();
            iTrackSet = setOfTracks.begin();
            for (; iTrackSet != setOfTracks.end(); ++iTrackSet,++iRetval) {
                iRetval->first   = TrackIDToParticle(iTrackSet->first);
                iRetval->second  = iTrackSet->second / Etot;
            }

            return retval;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::Track>>
        BackTrackerCore::MCParticleToTracks(simb::MCParticle* const p,
            std::vector<art::Ptr<rec::Track>> const& tracks) {

            // Can't be const because it uses TrackToMCParticles
            std::vector<art::Ptr<rec::Track>> retval;
            std::vector<std::pair<simb::MCParticle*,float>> fromMatch;
            int desiredTrackId = p->TrackId();

            for (art::Ptr<rec::Track> iTrack : tracks) {
                // de-ref the art::Ptr and create a bare pointer for
                // TrackToMCParticles call, managing const-correctness as best we can
                rec::Track* argument = const_cast<rec::Track*>( &(*iTrack) );
                fromMatch = TrackToMCParticles( argument );
                for (std::pair<simb::MCParticle*,float> matches : fromMatch) {
                    if ( matches.first->TrackId() == desiredTrackId &&
                          matches.second           >  fTrackFracMCP) {
                        retval.push_back( iTrack );
                    }
                }
            }
            return retval;
        }





        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::CaloHit>> const
        BackTrackerCore::ClusterToCaloHits(rec::Cluster* const c) {

			// Someday, could convert the const std::vector<gar::rec::CaloHit*>&
			// produced by Cluster::CalorimeterHits() into a
			// std::vector<art::Ptr<gar::rec::CaloHit>> which is the required
			// signature here, delete the construction of fClusterIDToCaloHits,
			// and make this method const.  Someday.

            // Evidently, use of std::unordered_map keeps this method from being const
            std::vector<art::Ptr<gar::rec::CaloHit>> retval;
            if ( !fClusterIDToCaloHits.empty() ) {
                retval = fClusterIDToCaloHits[ c->getIDNumber() ];
            }
            return retval;
        }



        //----------------------------------------------------------------------
        std::vector<std::pair<simb::MCParticle*,float>>
        BackTrackerCore::ClusterToMCParticles(rec::Cluster* const c) {
            // Can't be const because it uses ClusterToCaloHits

            std::vector<art::Ptr<gar::rec::CaloHit>> const hitCol = ClusterToCaloHits(c);
            std::vector<CalIDE> hitIDECol;
            for (auto hit : hitCol) {
                std::vector<CalIDE> IDEsThisHit = CaloHitToCalIDEs(hit);
                hitIDECol.insert(hitIDECol.end(), IDEsThisHit.begin(),IDEsThisHit.end());
            }

            TVector3 verty(c->Position());
            bool inYoke = ( fGeo->PointInMuIDBarrel(verty) ||
                            fGeo->PointInMuIDEndcap(verty) );

            // A tree is a good way to accumulate the info we need
            std::map<int,float> lClusterIdToEnergy;
            for (auto calIDE : hitIDECol) {
                simb::MCParticle* thisTrack     = TrackIDToParticle(calIDE.trackID);
                simb::MCParticle* incomingTrack;
                incomingTrack = inYoke ? FindECALEve(thisTrack) : FindTPCEve(thisTrack);
                int iTrackID = incomingTrack->TrackId();
                lClusterIdToEnergy[iTrackID] += calIDE.energyFrac * calIDE.energyTot;
            }

            // Next, sort the map by value.
            // Declaring the type of the sorting predicate
            typedef std::function<bool(std::pair<int,float>, std::pair<int,float>)> Comparator;
            // Declaring a set that will store the pairs using above comparison logic
            std::set<std::pair<int,float>, Comparator> setOfTracks(
                lClusterIdToEnergy.begin(), lClusterIdToEnergy.end(),
                [](std::pair<int,float> a ,std::pair<int,float> b) {
                    return a.second > b.second;
                }
            );

            // Make the returned vector.  Could have no tracks at all?
            float Etot = 0.0;
            auto iTrackSet = setOfTracks.begin();
            for (; iTrackSet!=setOfTracks.end(); ++iTrackSet) {
                Etot += iTrackSet->second;   // Can't do this in loop on lClusterIdToEnergy
            }

            std::pair<simb::MCParticle*,float> emptee(nullptr,0.0);
            std::vector<std::pair<simb::MCParticle*,float>> retval(lClusterIdToEnergy.size(),emptee);

            std::vector<std::pair<simb::MCParticle*,float>>::iterator iRetval = retval.begin();
            iTrackSet = setOfTracks.begin();
            for (; iTrackSet != setOfTracks.end(); ++iTrackSet,++iRetval) {
                iRetval->first   = TrackIDToParticle(iTrackSet->first);
                iRetval->second  = iTrackSet->second / Etot;
            }

            return retval;
        }



        //----------------------------------------------------------------------
        std::vector<art::Ptr<rec::Cluster>>
        BackTrackerCore::MCParticleToClusters(simb::MCParticle* const p,
                             std::vector<art::Ptr<rec::Cluster>> const& clusters) {

            // Can't be const because it uses ClusterToMCParticles
            std::vector<art::Ptr<rec::Cluster>> retval;
            std::vector<std::pair<simb::MCParticle*,float>> fromMatch;
            int desiredTrackId = p->TrackId();

            for (art::Ptr<rec::Cluster> iCluster : clusters) {
               // de-ref the art::Ptr and create a bare pointer for
               // TrackToMCParticles call, managing const-correctness as best we can
               rec::Cluster* argument = const_cast<rec::Cluster*>( &(*iCluster) );
               fromMatch = ClusterToMCParticles( argument );
               for (std::pair<simb::MCParticle*,float> matches : fromMatch) {
                   if ( matches.first->TrackId() == desiredTrackId &&
                        matches.second           >  fClusterFracMCP) {
                       retval.push_back( iCluster );
                   }
               }
            }
            return retval;
        }



        //----------------------------------------------------------------------
        bool BackTrackerCore::ClusterCreatedMCParticle(simb::MCParticle* const p,
                                                       rec::Cluster* const c) {
            std::vector<simb::MCParticle*> partsIsParts = MCPartsInCluster(c);
            for (simb::MCParticle* clusPart : partsIsParts) {
                if (IsForebearOf(clusPart,p)) return true;
            }
            return false;
        }



        //----------------------------------------------------------------------
        bool BackTrackerCore::MCParticleCreatedCluster(simb::MCParticle* const p,
                                                       rec::Cluster* const c) {
            std::vector<simb::MCParticle*> partsIsParts = MCPartsInCluster(c);

            for (simb::MCParticle* clusPart : partsIsParts) {
                if (IsForebearOf(p,clusPart)) return true;
            }
            return false;
        }



        //----------------------------------------------------------------------
        std::vector<simb::MCParticle*>
        BackTrackerCore::MCPartsInCluster(rec::Cluster* const c) {

            std::vector<int> trackIdsInCluster;
            std::vector<art::Ptr<rec::CaloHit>> const caloHits = ClusterToCaloHits(c);
            for (art::Ptr iCaloHit : caloHits) {
                std::vector<CalIDE> IDEsThisHit = CaloHitToCalIDEs(iCaloHit);
                for ( CalIDE iIDE : IDEsThisHit ) trackIdsInCluster.push_back(iIDE.trackID);
            }

            std::sort(trackIdsInCluster.begin(),trackIdsInCluster.end());
            std::vector<int>::iterator itr =
                std::unique(trackIdsInCluster.begin(),trackIdsInCluster.end());
            trackIdsInCluster.resize(std::distance(trackIdsInCluster.begin(),itr));

            std::vector<simb::MCParticle*> retval;
            for ( int trackId : trackIdsInCluster) {
                simb::MCParticle* pahtay = TrackIDToParticle(trackId);
                retval.push_back(pahtay);
            }
            return retval;
        }





        //----------------------------------------------------------------------
        TLorentzVector BackTrackerCore::EnergyDepositToMomentum(const int& trackID, 
                              const TLorentzVector& position, size_t& startTrajIndex) const {

            auto const& mcp = TrackIDToParticle(trackID); // MCParticle
            float dr=FLT_MAX, dt=FLT_MAX; // spatial and temporal distance

            // stepping through the trajectory, the distance should decrease until
            // min dist found (could be starting point)
            for(;startTrajIndex<mcp->NumberTrajectoryPoints(); startTrajIndex++){

                if(   (position - mcp->Position(startTrajIndex)).Vect().Mag() < dr &&
                  abs((position - mcp->Position(startTrajIndex)).T()) < dt) {

                    dr = (position - mcp->Position(startTrajIndex)).Vect().Mag();
                    dt = abs((position - mcp->Position(startTrajIndex)).T());
                } else {
                    break;
                }

            } // for trajectory points

            return  mcp->Momentum(startTrajIndex);

        } // def EnergyDepositToMomentum()





    } // cheat
} // gar
