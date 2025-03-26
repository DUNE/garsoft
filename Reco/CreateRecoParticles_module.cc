////////////////////////////////////////////////////////////////////////
// Class:       CreateRecoParticles
// Plugin Type: producer (art v)
// File:        CreateRecoParticles_module.cc
//
// Takes Tracks and Clusters to construct RecoParticles, which contain
// information about the (calibrated) HPgTPC dE/dx, ECal muon scores
// and time-of-flight
////////////////////////////////////////////////////////////////////////

// C++ Includes
#include <iostream>
#include <memory>
#include <vector> // std::ostringstream

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "cetlib/search_path.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"
#include "Utilities/AssociationUtil.h"

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/RecoParticle.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/Vertex.h"

#include "RecoAlg/ECALMuonBDT.h"
#include "RecoAlg/ECALToFAlg.h"
#include "RecoAlg/TruncatedIonizationCalculator.h"

namespace gar {
  namespace rec {

    class CreateRecoParticles : public ::art::EDProducer {
    public:
      /// Standard constructor and destructor for an FMWK module.
      explicit CreateRecoParticles(fhicl::ParameterSet const& pset);
      virtual ~CreateRecoParticles();

      // Plugins should not be copied or assigned.
      CreateRecoParticles(CreateRecoParticles const&) = delete;
      CreateRecoParticles(CreateRecoParticles&&) = delete;
      CreateRecoParticles& operator=(CreateRecoParticles const&) = delete;
      CreateRecoParticles& operator=(CreateRecoParticles&&) = delete;

      void produce(art::Event& evt) override;

      void reconfigure(fhicl::ParameterSet const& pset);

    private:
      typedef int ClusterId;

      void PrepareClusterMap(art::Handle<std::vector<rec::Cluster>>& RecoClusterHandle,
                             std::unordered_map<ClusterId, Int_t>& IdToIndexMap);

      float FindTZero(art::Handle<std::vector<rec::Cluster>>& RecoECALHandle,
                      std::unordered_map<ClusterId, Int_t>& ECALIdToIndex,
                      art::FindManyP<gar::rec::CaloHit>* findManyECALRecoHit);

      void ParticlesFromTracks(
        const art::Event& evt,
        std::vector<rec::RecoParticle*>& particleVector,
        std::unordered_map<size_t, art::Ptr<rec::Track>>& MapRecoParticleTrackPtr,
        std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>>& MapRecoParticleECalPtr,
        std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>>& MapRecoParticleMuIDPtr);

      void GuessParticleDirection(std::vector<rec::RecoParticle*>& particleVector,
                                  art::Handle<std::vector<rec::Track>>& RecoTrackHandle);

      std::string fTrackLabel;         ///< module label for TPC Tracks rec:Track
      std::string fClusterLabel;       ///< module label for calo clusters rec::Cluster
      std::string fECALAssnLabel;      ///< module label for track-clusters associations
      std::string fClusterLabel_MuID;  ///< module label for calo clusters rec::Cluster
      std::string fECALAssnLabel_MuID; ///< module label for track-clusters associations
      std::string fVertexLabel;        ///< module label for vertexes rec:Vertex

      std::string fInstanceLabelCalo;      ///< Instance name for ECAL
      std::string fInstanceLabelCalo_MuID; ///< Instance name for MuID

      const gar::geo::GeometryCore* fGeo; ///< geometry information

      std::unique_ptr<gar::rec::alg::TruncatedIonizationCalculator>
        fTPCIonizationAlg; ///< algorithm to compute TPC mean ionization
      std::unique_ptr<gar::rec::alg::ECALMuonBDT>
        fECALMuonBDT; ///< algorithm to compute ECAL/MuID muon scores
      std::unique_ptr<gar::rec::alg::ECALToFAlg>
        fECALToFAlg; ///< algorithm to compute ECAL arrival times
    };

  } // namespace rec

  namespace rec {

    //----------------------------------------------------------------------
    // Constructor
    CreateRecoParticles::CreateRecoParticles(fhicl::ParameterSet const& pset)
      : art::EDProducer{pset}
    {
      fGeo = gar::providerFrom<geo::GeometryGAr>();

      this->reconfigure(pset);

      consumes<std::vector<rec::Track>>(fTrackLabel);
      consumes<art::Assns<rec::Track, rec::Vertex>>(fVertexLabel);

      produces<std::vector<gar::rec::RecoParticle>>();

      produces<art::Assns<rec::RecoParticle, rec::Track>>();
      produces<art::Assns<rec::RecoParticle, rec::Cluster>>(fInstanceLabelCalo);
      produces<art::Assns<rec::RecoParticle, rec::Cluster>>(fInstanceLabelCalo_MuID);

      return;
    }

    //----------------------------------------------------------------------
    // Destructor
    CreateRecoParticles::~CreateRecoParticles() {}

    //----------------------------------------------------------------------
    void CreateRecoParticles::reconfigure(fhicl::ParameterSet const& pset)
    {
      MF_LOG_DEBUG("CreateRecoParticles") << "Debug: CreateRecoParticles()";

      fTrackLabel = pset.get<std::string>("TrackLabel", "track");
      fClusterLabel = pset.get<std::string>("ClusterLabel", "calocluster");
      fECALAssnLabel = pset.get<std::string>("ECALAssnLabel", "trkecalassn");
      fInstanceLabelCalo = pset.get<std::string>("InstanceLabelCalo", "ECAL");
      fClusterLabel_MuID = pset.get<std::string>("ClusterLabel_MuID", "caloclustermuid");
      fECALAssnLabel_MuID = pset.get<std::string>("ECALAssnLabel_MuID", "trkecalassnmuid");
      fInstanceLabelCalo_MuID = pset.get<std::string>("InstanceLabelCalo_MuID", "MuID");
      fVertexLabel = pset.get<std::string>("VertexLabel", "vertex");

      // Configuring TPCIonizationAlg
      auto TPCIonizationAlgPars = pset.get<fhicl::ParameterSet>("TPCIonizationAlgPars");
      fTPCIonizationAlg =
        std::make_unique<gar::rec::alg::TruncatedIonizationCalculator>(TPCIonizationAlgPars);

      // Configuring ECALMuonBDT
      auto ECALMuonBDTPars = pset.get<fhicl::ParameterSet>("ECALMuonBDTPars");
      fECALMuonBDT = std::make_unique<gar::rec::alg::ECALMuonBDT>(ECALMuonBDTPars, fGeo);

      // Configuring ECALToFAlg
      auto ECALToFAlgPars = pset.get<fhicl::ParameterSet>("ECALToFAlgPars");
      fECALToFAlg = std::make_unique<gar::rec::alg::ECALToFAlg>(ECALToFAlgPars, fGeo);

      return;
    }

    //--------------------------------------------------------------------------
    void CreateRecoParticles::produce(::art::Event& evt)
    {
      MF_LOG_DEBUG("CreateRecoParticles") << "produce()";

      // Output collections
      std::unique_ptr<std::vector<rec::RecoParticle>> RecoParticleCol(
        new std::vector<rec::RecoParticle>);

      std::unique_ptr<art::Assns<rec::RecoParticle, rec::Track>> RecoParticleTrkAssns(
        new ::art::Assns<rec::RecoParticle, rec::Track>);
      std::unique_ptr<art::Assns<rec::RecoParticle, rec::Cluster>> RecoParticleECalAssns(
        new ::art::Assns<rec::RecoParticle, rec::Cluster>);
      std::unique_ptr<art::Assns<rec::RecoParticle, rec::Cluster>> RecoParticleMuIDAssns(
        new ::art::Assns<rec::RecoParticle, rec::Cluster>);

      auto const particlePtrMaker = art::PtrMaker<rec::RecoParticle>(evt);

      // Create unordered maps to collect associations
      std::unordered_map<size_t, art::Ptr<rec::Track>> MapRecoParticleTrackPtr;
      std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>> MapRecoParticleECalPtr;
      std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>> MapRecoParticleMuIDPtr;

      std::vector<rec::RecoParticle*> particleVector;
      this->ParticlesFromTracks(evt,
                                particleVector,
                                MapRecoParticleTrackPtr,
                                MapRecoParticleECalPtr,
                                MapRecoParticleMuIDPtr);

      size_t iRecoParticle = 0;
      for (auto const& it : particleVector) {
        RecoParticleCol->push_back(*it);

        if (MapRecoParticleTrackPtr.find(iRecoParticle) != MapRecoParticleTrackPtr.end()) {
          auto const particlePtr = particlePtrMaker(iRecoParticle);
          auto const trackPtr = MapRecoParticleTrackPtr[iRecoParticle];
          RecoParticleTrkAssns->addSingle(particlePtr, trackPtr);
        }

        if (MapRecoParticleECalPtr.find(iRecoParticle) != MapRecoParticleECalPtr.end()) {
          auto const particlePtr = particlePtrMaker(iRecoParticle);
          size_t n_ecal_assns = MapRecoParticleECalPtr[iRecoParticle].size();
          for (size_t iAssnECal = 0; iAssnECal < n_ecal_assns; ++iAssnECal) {
            auto const ecalPtr = MapRecoParticleECalPtr[iRecoParticle].at(iAssnECal);
            RecoParticleECalAssns->addSingle(particlePtr, ecalPtr);
          }
        }

        if (MapRecoParticleMuIDPtr.find(iRecoParticle) != MapRecoParticleMuIDPtr.end()) {
          auto const particlePtr = particlePtrMaker(iRecoParticle);
          size_t n_muid_assns = MapRecoParticleMuIDPtr[iRecoParticle].size();
          for (size_t iAssnMuID = 0; iAssnMuID < n_muid_assns; ++iAssnMuID) {
            auto const muidPtr = MapRecoParticleMuIDPtr[iRecoParticle].at(iAssnMuID);
            RecoParticleMuIDAssns->addSingle(particlePtr, muidPtr);
          }
        }

        iRecoParticle++;
      }

      evt.put(std::move(RecoParticleCol));

      evt.put(std::move(RecoParticleTrkAssns));
      evt.put(std::move(RecoParticleECalAssns), fInstanceLabelCalo);
      evt.put(std::move(RecoParticleMuIDAssns), fInstanceLabelCalo_MuID);

      return;
    } // CreateRecoParticles::produce()

    //--------------------------------------------------------------------------
    void CreateRecoParticles::PrepareClusterMap(
      art::Handle<std::vector<rec::Cluster>>& RecoClusterHandle,
      std::unordered_map<ClusterId, Int_t>& IdToIndexMap)
    {

      Int_t iCluster = 0;
      for (auto const& cluster : (*RecoClusterHandle)) {
        // Fill the map between IDs and ECAL indeces
        int Id = cluster.getIDNumber();
        IdToIndexMap[Id] = iCluster++;
      } // end loop over Clusters
    }

    //--------------------------------------------------------------------------
    float CreateRecoParticles::FindTZero(art::Handle<std::vector<rec::Cluster>>& RecoECALHandle,
                                         std::unordered_map<ClusterId, Int_t>& ECALIdToIndex,
                                         art::FindManyP<gar::rec::CaloHit>* findManyECALRecoHit)
    {

      // Get index of cluster with earliest hit time
      float t0 = 99999.0;
      int t0_cluster_index = -1;
      for (auto const& ecal_cluster : (*RecoECALHandle)) {

        int ECALIndex = ECALIdToIndex[ecal_cluster.getIDNumber()];

        if (findManyECALRecoHit->isValid()) {

          int nECALClusterHit = findManyECALRecoHit->at(ECALIndex).size();

          for (int iECALClusterHit = 0; iECALClusterHit < nECALClusterHit; ++iECALClusterHit) {

            const rec::CaloHit ecal_hit = *(findManyECALRecoHit->at(ECALIndex).at(iECALClusterHit));
            float time = ecal_hit.Time().first;

            if (time <= t0) {
              t0 = time;
              t0_cluster_index = ECALIndex;
            }
          }
        }
      } // end loop over Clusters

      if (t0_cluster_index != -1) { return t0; }
      else {
        return 0.0; // in case there is no ECAL hit with time less than initial guess
      }
    }

    //--------------------------------------------------------------------------
    void CreateRecoParticles::ParticlesFromTracks(
      const art::Event& evt,
      std::vector<rec::RecoParticle*>& particleVector,
      std::unordered_map<size_t, art::Ptr<rec::Track>>& MapRecoParticleTrackPtr,
      std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>>& MapRecoParticleECalPtr,
      std::unordered_map<size_t, std::vector<art::Ptr<rec::Cluster>>>& MapRecoParticleMuIDPtr)
    {

      // Get handle for tracks
      auto TrackHandle = evt.getHandle<std::vector<rec::Track>>(fTrackLabel);
      if (!TrackHandle) {
        throw cet::exception("CreateRecoParticles")
          << " No rec::Track branch."
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Get the tags for the ECAL clusters and associations
      art::InputTag ecalclustertag(fClusterLabel, fInstanceLabelCalo);
      art::InputTag ecalassntag(fECALAssnLabel, fInstanceLabelCalo);

      // Get the tags for the MuID clusters and associations
      art::InputTag ecalclustertagmuid(fClusterLabel_MuID, fInstanceLabelCalo_MuID);
      art::InputTag ecalassntagmuid(fECALAssnLabel_MuID, fInstanceLabelCalo_MuID);

      // Get associations between tracks and track ionization objects
      art::FindOneP<rec::TrackIoniz>* findIonization =
        new art::FindOneP<rec::TrackIoniz>(TrackHandle, evt, fTrackLabel);

      // Get handle for ECAL clusters
      art::Handle<std::vector<rec::Cluster>> RecoECALHandle;

      RecoECALHandle = evt.getHandle<std::vector<rec::Cluster>>(ecalclustertag);
      if (!RecoECALHandle) {
        throw cet::exception("ECALAna")
          << " No rec::Cluster branch."
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Get handle for MuID clusters
      art::Handle<std::vector<rec::Cluster>> RecoMuIDHandle;

      RecoMuIDHandle = evt.getHandle<std::vector<rec::Cluster>>(ecalclustertagmuid);
      if (!RecoMuIDHandle) {
        throw cet::exception("ECALAna")
          << " No rec::Cluster MuID branch."
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Get associations between tracks and ECAL clusters and between ECAL clusters and ECAL hits
      art::FindManyP<rec::Cluster, rec::TrackEnd>* findManyTrackECAL =
        new art::FindManyP<rec::Cluster, rec::TrackEnd>(TrackHandle, evt, ecalassntag);
      art::FindManyP<gar::rec::CaloHit>* findManyECALRecoHit =
        new art::FindManyP<gar::rec::CaloHit>(RecoECALHandle, evt, ecalclustertag);

      // Get associations between tracks and MuID clusters and between MuID clusters and MuID hits
      art::FindManyP<rec::Cluster, rec::TrackEnd>* findManyTrackMuID =
        new art::FindManyP<rec::Cluster, rec::TrackEnd>(TrackHandle, evt, ecalassntagmuid);
      art::FindManyP<gar::rec::CaloHit>* findManyMuIDRecoHit =
        new art::FindManyP<gar::rec::CaloHit>(RecoMuIDHandle, evt, ecalclustertagmuid);

      // Get associations between Tracks and Vertices
      art::FindManyP<rec::Vertex, rec::TrackEnd>* findManyTracksVertices = NULL;
      findManyTracksVertices =
        new art::FindManyP<rec::Vertex, rec::TrackEnd>(TrackHandle, evt, fVertexLabel);

      // Fill maps between Id and index for ECAL and MuID clusters
      std::unordered_map<ClusterId, Int_t> ECALIdToIndex;
      PrepareClusterMap(RecoECALHandle, ECALIdToIndex);

      std::unordered_map<ClusterId, Int_t> MuIDIdToIndex;
      PrepareClusterMap(RecoMuIDHandle, MuIDIdToIndex);

      // Compute t0 estimate
      float t0 = FindTZero(RecoECALHandle, ECALIdToIndex, findManyECALRecoHit);

      // Start track loop
      for (size_t iTrack = 0; iTrack < TrackHandle->size(); ++iTrack) {

        const art::Ptr<rec::Track> track(TrackHandle, iTrack);
        const rec::Track* track_ptr = track.get();

        MapRecoParticleTrackPtr[iTrack] = track;

        rec::RecoParticle* particle = new rec::RecoParticle();

        float momentum =
          0.5 * (track->Momentum_beg() + track->Momentum_end()); // get mean of both fits

        particle->setMomentum(momentum);

        if (findIonization->isValid()) {
          const rec::TrackIoniz* ionization_ptr = (findIonization->at(iTrack)).get();

          fTPCIonizationAlg->PrepareAlgo(track_ptr, ionization_ptr);
          fTPCIonizationAlg->ComputeMeanIonization();
          std::pair<float, float> IonizationInfo = fTPCIonizationAlg->GetIonization();

          particle->setTotalCaloEnergy(IonizationInfo.first);
          particle->setMeanCaloEnergy(IonizationInfo.second);

          particle->setProtondEdxScore(fTPCIonizationAlg->GetdEdxProtonScore());
        }

        fECALMuonBDT->PrepareAlgo(track_ptr);

        int nECALedTracks = 0;
        rec::TrackEnd iEnd = gar::rec::TrackEndBeg;
        if (findManyTrackECAL->isValid()) {
          nECALedTracks = findManyTrackECAL->at(iTrack).size();
          if (nECALedTracks > 0) {
            iEnd = *(findManyTrackECAL->data(iTrack).at(
              0)); // in trackecalassns2 the used track end is always the same for one track...
            particle->setTrackEndECALed(
              iEnd); // just fill this if there are associations, if not it's just a -1
          }
        }

        // A ToF measurement is only possible if we can propagate the track to the ECAL
        // This is checked when preparing the ToF algorithm
        bool tof_possible = fECALToFAlg->PrepareAlgo(track_ptr, iEnd, t0);

        // Start loop over Track - ECAL cluster associations
        for (int iECALedTrack = 0; iECALedTrack < nECALedTracks; ++iECALedTrack) {

          const art::Ptr<rec::Cluster> ecal_cluster =
            findManyTrackECAL->at(iTrack).at(iECALedTrack);
          const rec::Cluster* ecal_cluster_ptr = ecal_cluster.get();

          // Add ECal cluster pointer to map to add associations later
          MapRecoParticleECalPtr[iTrack].push_back(ecal_cluster);

          //Get the associated reco hits
          size_t ECALIndex = ECALIdToIndex[ecal_cluster_ptr->getIDNumber()];
          std::vector<const rec::CaloHit*> ecal_hit_ptr_vec;

          if (findManyECALRecoHit->isValid()) {
            int nECALClusterHit = findManyECALRecoHit->at(ECALIndex).size();

            for (int iECALClusterHit = 0; iECALClusterHit < nECALClusterHit; ++iECALClusterHit) {

              const rec::CaloHit* ecal_hit_ptr =
                findManyECALRecoHit->at(ECALIndex).at(iECALClusterHit).get();
              ecal_hit_ptr_vec.push_back(ecal_hit_ptr);
            }
          }

          fECALMuonBDT->AddECALHits(ecal_cluster_ptr, ecal_hit_ptr_vec);
          if (tof_possible) fECALToFAlg->AddHits(ecal_hit_ptr_vec);

        } // end loop over Track - ECAL cluster associations

        // For the ToF measurement we need hits in the inner layers of the ECAL
        if (tof_possible) {
          fECALToFAlg->ComputeArrivalTime();

          particle->setECALToFTime(fECALToFAlg->GetTime());
          particle->setECALToFBeta(fECALToFAlg->GetBeta());
          particle->setECALToFMass(fECALToFAlg->GetMass());

          particle->setProtonToFScore(fECALToFAlg->GetToFProtonScore());
        }

        int nMuIDedTracks = 0;
        if (findManyTrackMuID->isValid()) { nMuIDedTracks = findManyTrackMuID->at(iTrack).size(); }

        // Start loop over Track - MuID cluster associations
        for (int iMuIDedTrack = 0; iMuIDedTrack < nMuIDedTracks; ++iMuIDedTrack) {

          const art::Ptr<rec::Cluster> muid_cluster =
            findManyTrackMuID->at(iTrack).at(iMuIDedTrack);
          const rec::Cluster* muid_cluster_ptr = muid_cluster.get();

          // Add MuID cluster pointer to map to add associations later
          MapRecoParticleMuIDPtr[iTrack].push_back(muid_cluster);

          //Get the associated reco hits
          size_t MuIDIndex = MuIDIdToIndex[muid_cluster_ptr->getIDNumber()];
          std::vector<const rec::CaloHit*> muid_hit_ptr_vec;

          if (findManyMuIDRecoHit->isValid()) {
            int nMuIDClusterHit = findManyMuIDRecoHit->at(MuIDIndex).size();

            for (int iMuIDClusterHit = 0; iMuIDClusterHit < nMuIDClusterHit; ++iMuIDClusterHit) {

              const rec::CaloHit* muid_hit_ptr =
                findManyMuIDRecoHit->at(MuIDIndex).at(iMuIDClusterHit).get();
              muid_hit_ptr_vec.push_back(muid_hit_ptr);
            }
          }

          fECALMuonBDT->AddMuIDHits(muid_cluster_ptr, muid_hit_ptr_vec);
        } // end loop over Track - MuID cluster associations

        fECALMuonBDT->ComputeFeatures();
        fECALMuonBDT->ApplyClassifier();

        std::pair<float, int> ecal_energy = fECALMuonBDT->GetECALEnergy();
        particle->setTotalECALEnergy(ecal_energy.first);
        particle->setNHitsECAL(ecal_energy.second);

        std::pair<float, int> muid_energy = fECALMuonBDT->GetMuIDEnergy();
        particle->setTotalMuIDEnergy(muid_energy.first);
        particle->setNHitsMuID(muid_energy.second);

        particle->setMuonScore(fECALMuonBDT->GetScore());
        particleVector.push_back(particle);

        if (findManyTracksVertices->isValid()) {

          int nVertex = findManyTracksVertices->at(iTrack).size();

          // This is not the best choice probably, a track can be associated to more than one vertex (methinks)
          // Or at least it should be possible that both its track ends get vertexed
          if (nVertex > 0)
            particle->setTrackEndVertexed(*(findManyTracksVertices->data(iTrack).at(
              0))); // just fill this if there are associations, if not it's just a -1

        } // end if Tracks - Vertices Associations isValid
      }

      GuessParticleDirection(particleVector, TrackHandle);
    }

    //--------------------------------------------------------------------------
    void CreateRecoParticles::GuessParticleDirection(
      std::vector<rec::RecoParticle*>& particleVector,
      art::Handle<std::vector<rec::Track>>& RecoTrackHandle)
    {
      // We need to find the TrackEnd position of the highest momentum ECaled particle
      float highest_momentum_ecaled = -1.0;
      float position_trackend_ecaled[3] = {0.0};
      // Similar but for the highest momentum vertexed particle
      float highest_momentum_vertexed = -1.0;
      float position_trackend_vertexed[3] = {0.0};

      size_t iTrack = 0;

      // Start loop over RecoParticles
      for (auto const& particle : particleVector) {

        // Check first for ECaled particles...
        if ((particle->NHitsECAL() != 0) && (particle->Momentum() >= highest_momentum_ecaled)) {

          const art::Ptr<rec::Track> track(RecoTrackHandle, iTrack);
          const rec::Track* track_ptr = track.get();

          // If the track end ECALed is the End (0), then the true begin is the Begin (0)
          if (particle->TrackEndECALed() == gar::rec::TrackEndBeg) {

            for (int i = 0; i < 3; ++i)
              position_trackend_ecaled[i] = track_ptr->End()[i];

            // Else, if the end ECALed is the Begin (1), the true begin is the End (0)
          }
          else if (particle->TrackEndECALed() == gar::rec::TrackEndEnd) {
            for (int i = 0; i < 3; ++i)
              position_trackend_ecaled[i] = track_ptr->Vertex()[i];
          }

          highest_momentum_ecaled = particle->Momentum();
        }

        // ...and then for vertexed particles
        if ((particle->TrackEndVertexed() != -1) &&
            (particle->Momentum() >= highest_momentum_vertexed)) {

          const art::Ptr<rec::Track> track(RecoTrackHandle, iTrack);
          const rec::Track* track_ptr = track.get();

          // If the track end Vertexed is the Begin (1), then the true begin is the Begin (1)
          if (particle->TrackEndVertexed() == gar::rec::TrackEndBeg) {
            for (int i = 0; i < 3; ++i)
              position_trackend_vertexed[i] = track_ptr->Vertex()[i];

            // Else, if the end Vertexed is the End (0), the true begin is the End (0)
          }
          else if (particle->TrackEndVertexed() == gar::rec::TrackEndEnd) {
            for (int i = 0; i < 3; ++i)
              position_trackend_vertexed[i] = track_ptr->End()[i];
          }

          highest_momentum_vertexed = particle->Momentum();
        }

        iTrack++;

      } // end loop over RecoParticles

      float candidate_trackend_vertexed[3];
      // Use as reference position the true begin of the highest momentum ECaled particle
      if (highest_momentum_ecaled > 0.0) {
        for (int i = 0; i < 3; ++i)
          candidate_trackend_vertexed[i] = position_trackend_ecaled[i];
        // If no particle was ECaled use the position of the highest momentum vertexed particle
      }
      else if (highest_momentum_vertexed > 0.0) {
        for (int i = 0; i < 3; ++i)
          candidate_trackend_vertexed[i] = position_trackend_vertexed[i];
        // If no particle in the event was ECaled or vertexed we cannot guess the direction
        // of the particles, therefore we cannot know their charge!
      }
      else {
        return;
      }

      // Use reference position to compute 3D distance to the tracks start and end points
      // Then assign charge based on what point is closer to reference position
      iTrack = 0;

      for (auto const& particle : particleVector) {

        const art::Ptr<rec::Track> track(RecoTrackHandle, iTrack);
        const rec::Track* track_ptr = track.get();

        float distance_track_begin =
          std::hypot(candidate_trackend_vertexed[0] - track_ptr->Vertex()[0],
                     candidate_trackend_vertexed[1] - track_ptr->Vertex()[1],
                     candidate_trackend_vertexed[2] - track_ptr->Vertex()[2]);

        float distance_track_end = std::hypot(candidate_trackend_vertexed[0] - track_ptr->End()[0],
                                              candidate_trackend_vertexed[1] - track_ptr->End()[1],
                                              candidate_trackend_vertexed[2] - track_ptr->End()[2]);

        if (distance_track_begin <= distance_track_end) {
          particle->setCharge(track_ptr->ChargeBeg());
          particle->setStart(track_ptr->Vertex());
          particle->setEnd(track_ptr->End());
          particle->setDirection(track_ptr->VtxDir());
        }
        else {
          particle->setCharge(track_ptr->ChargeEnd());
          particle->setStart(track_ptr->End());
          particle->setEnd(track_ptr->Vertex());
          particle->setDirection(track_ptr->EndDir());
        }

        iTrack++;

      } // end loop over RecoParticles
    }

  } // namespace rec

  namespace rec {

    DEFINE_ART_MODULE(CreateRecoParticles)

  } // namespace rec
} // gar
