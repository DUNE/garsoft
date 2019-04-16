////////////////////////////////////////////////////////////////////////
/// \file  CosmicsGen_plugin.cc
/// \brief Generator for cosmic-rays
///
/// Module to produce cosmic ray MC using CRY
///
/// \version $Id: CosmicsGen.cxx,v 1.5 2010/04/23 18:37:36 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVGEN_COSMICSGEN_H
#define EVGEN_COSMICSGEN_H

// ROOT includes
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "nutools/EventGeneratorBase/CRY/CRYHelper.h"
#include "nutools/RandomUtils/NuRandomService.h"

// garsoft includes
#include "Geometry/Geometry.h"
#include "SummaryDataProducts/RunData.h"
#include "CoreUtils/ServiceUtil.h"

namespace gar{
  namespace evgen {
    
    /// A module to check the results from the Monte Carlo generator
    class CosmicsGen : public ::art::EDProducer {
    public:
      explicit CosmicsGen(fhicl::ParameterSet const& pset);
      virtual ~CosmicsGen();
      
      
      void produce(::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      void reconfigure(fhicl::ParameterSet const& p);
      
    private:
      
      evgb::CRYHelper* fCRYHelp; ///< CRY generator object
      
      std::vector<double> fbuffbox;
      
      TH2F* fPhotonAngles;       ///< Photon rate vs angle
      TH2F* fPhotonAnglesLo;     ///< Photon rate vs angle, low momenta
      TH2F* fPhotonAnglesMi;     ///< Photon rate vs angle, middle momenta
      TH2F* fPhotonAnglesHi;     ///< Photon rate vs angle, high momenta
      TH1F* fPhotonCosQ;         ///< Photon rate vs cos(Q)
      TH1F* fPhotonEnergy;       ///< Photon energy (GeV)
      TH1F* fPhotonsPerSample;   ///< number of photons in the sampled time window
      TH1F* fPhotonsInTPC;       ///< number of photons in the tpc during
                                 ///< the sampled time window
      
      TH2F* fElectronAngles;     ///< Electron rate vs angle
      TH2F* fElectronAnglesLo;   ///< Electron rate vs angle, low momenta
      TH2F* fElectronAnglesMi;   ///< Electron rate vs angle, middle momenta
      TH2F* fElectronAnglesHi;   ///< Electron rate vs angle, high momenta
      TH1F* fElectronCosQ;       ///< Electron rate vs cos(Q)
      TH1F* fElectronEnergy;     ///< Electron energy (GeV)
      TH1F* fElectronsPerSample; ///< number of electrons in the sampled time window
      TH1F* fElectronsInTPC;     ///< number of electrons in the tpc during
                                 ///< the sampled time window
      
      TH2F* fMuonAngles;         ///< Muon rate vs angle
      TH2F* fMuonAnglesLo;       ///< Muon rate vs angle, low momenta
      TH2F* fMuonAnglesMi;       ///< Muon rate vs angle, middle momenta
      TH2F* fMuonAnglesHi;       ///< Muon rate vs angle, high momenta
      TH1F* fMuonCosQ;           ///< Muon rate vs cos(Q)
      TH1F* fMuonEnergy;         ///< Muon energy (GeV)
      TH1F* fMuonsPerSample;     ///< number of muons in the sampled time window
      TH1F* fMuonsInTPC;         ///< number of muons in the tpc during
                                 ///< the sampled time window

      rndm::NuRandomService::seed_t fSeed;  ///< override seed with a fcl parameter not equal to zero
      
    };
  }
  
  namespace evgen{
    
    //____________________________________________________________________________
    CosmicsGen::CosmicsGen(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset}, fCRYHelp(0)
    {
      
      //the buffer box bounds specified here will extend on the cryostat boundaries
      fbuffbox = pset.get< std::vector<double> >("BufferBox",{0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
      fSeed          = pset.get< rndm::NuRandomService::seed_t >("Seed",0);
      this->reconfigure(pset);
      
      produces< std::vector<simb::MCTruth> >();
      produces< gar::sumdata::RunData, ::art::InRun >();
    }
    
    //____________________________________________________________________________
    CosmicsGen::~CosmicsGen()
    {
      if(fCRYHelp) delete fCRYHelp;
    }
    
    //____________________________________________________________________________
    void CosmicsGen::reconfigure(fhicl::ParameterSet const& p)
    {
      if(fCRYHelp){
        delete fCRYHelp;
        fCRYHelp = 0;
      }

    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this);
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    auto& engine = rng->getEngine(art::ScheduleID::first(),
                                  p.get<std::string>("module_label"));
    if (fSeed != 0) {
      engine.setSeed(fSeed, 0 /* dummy? */);
    }

      

      auto geo = gar::providerFrom<gar::geo::Geometry>();
      
      fCRYHelp = new evgb::CRYHelper(p, engine, geo->GetWorldVolumeName());
      
      return;
    }
    
    //____________________________________________________________________________
    void CosmicsGen::beginJob()
    {
      ::art::ServiceHandle<::art::TFileService> tfs;
      
      fPhotonAngles       = tfs->make<TH2F>("fPhotonAngles",        ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fPhotonAnglesLo     = tfs->make<TH2F>("fPhotonAnglesLo",      ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fPhotonAnglesMi     = tfs->make<TH2F>("fPhotonAnglesMi",      ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fPhotonAnglesHi     = tfs->make<TH2F>("fPhotonAnglesHi",      ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fPhotonCosQ         = tfs->make<TH1F>("fPhotonCosQ",          ";cos#theta;tracks",         50,-1.0,1.0);
      fPhotonEnergy       = tfs->make<TH1F>("fPhotonEnergy",        ";E (GeV)",                  5000,0.0,1000.0);
      fPhotonsPerSample   = tfs->make<TH1F>("fPhotonsPerSample",    ";Number Photons;Samples",   100, 0, 1000);
      fPhotonsInTPC       = tfs->make<TH1F>("fPhotonsInTPC",        ";Number Photons;Samples",   100, 0, 1000);
      
      fElectronAngles     = tfs->make<TH2F>("fElectronAngles",      ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fElectronAnglesLo   = tfs->make<TH2F>("fElectronAnglesLo",    ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fElectronAnglesMi   = tfs->make<TH2F>("fElectronAnglesMi",    ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fElectronAnglesHi   = tfs->make<TH2F>("fElectronAnglesHi",    ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fElectronCosQ       = tfs->make<TH1F>("fElectronCosQ",        ";cos#theta;tracks",         50,-1.0,1.0);
      fElectronEnergy     = tfs->make<TH1F>("fElectronEnergy",      ";E (GeV)",                  5000,0.0,1000.0);
      fElectronsPerSample = tfs->make<TH1F>("fElectronsPerSample",  ";Number Electrons;Samples", 100, 0, 1000);
      fElectronsInTPC     = tfs->make<TH1F>("fElectronsInTPC",      ";Number Electrons;Samples", 100, 0, 1000);
      
      fMuonAngles         = tfs->make<TH2F>("fMuonAngles",          ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fMuonAnglesLo       = tfs->make<TH2F>("fMuonAnglesLo",        ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fMuonAnglesMi       = tfs->make<TH2F>("fMuonAnglesMi",        ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fMuonAnglesHi       = tfs->make<TH2F>("fMuonAnglesHi",        ";#phi;cos#theta",           36,-180.0,180.0,50,-1.0,1.0);
      fMuonCosQ           = tfs->make<TH1F>("fMuonCosQ",            ";cos#theta;tracks",         50,-1.0,1.0);
      fMuonEnergy         = tfs->make<TH1F>("fMuonEnergy",          ";E (GeV)",                  5000,0.0,1000.0);
      fMuonsPerSample     = tfs->make<TH1F>("fMuonsPerSample",      ";Number Muons;Samples",     100, 0, 1000);
      fMuonsInTPC         = tfs->make<TH1F>("fMuonsInTPC",          ";Number Muons;Samples",     100, 0, 1000);
      
    }
    
    //____________________________________________________________________________
    void CosmicsGen::beginRun(::art::Run& run)
    {
      // grab the geometry object to see what geometry we are using
      auto geo = gar::providerFrom<geo::Geometry>();
      
      std::unique_ptr<gar::sumdata::RunData> runcol(new gar::sumdata::RunData(geo->DetectorName()));
      
      run.put(std::move(runcol));
      
      return;
    }
    
    //____________________________________________________________________________
    void CosmicsGen::produce(::art::Event& evt)
    {
      std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
      
      // fill some histograms about this event
      auto geom = gar::providerFrom<geo::Geometry>();
      
      // int nCrossCryostat = 0;   // no cryostat here.
      
      simb::MCTruth truth;
      
      // don't need to have a particle on every event

      //while(nCrossCryostat < 1){
        
        simb::MCTruth pretruth;
        truth.SetOrigin(simb::kCosmicRay);
	//std::cout << " calling CRY Helper: " << geom->SurfaceY() << " " << geom->DetLength() << std::endl;
        fCRYHelp->Sample(pretruth,
                         geom->SurfaceY(),
                         geom->DetLength(),
                         0);
        
        int numPhotons   = 0;
        int numElectrons = 0;
        int numMuons     = 0;
        int allPhotons   = 0;
        int allElectrons = 0;
        int allMuons     = 0;
        
          // loop over particles in the truth object
        for(int i = 0; i < pretruth.NParticles(); ++i){
          simb::MCParticle particle = pretruth.GetParticle(i);
          const TLorentzVector& p4 = particle.Momentum();
          
          if      (std::abs(particle.PdgCode())==13) ++allMuons;
          else if (std::abs(particle.PdgCode())==22) ++allPhotons;
          else if (std::abs(particle.PdgCode())==11) ++allElectrons;
          
          TH1F* hCosQ     = 0;
          TH2F* hAngles   = 0;
          TH2F* hAnglesLo = 0;
          TH2F* hAnglesMi = 0;
          TH2F* hAnglesHi = 0;
          TH1F* hEnergy   = 0;
          if (std::abs(particle.PdgCode())==13) {
            hCosQ     = fMuonCosQ;
            hAngles   = fMuonAngles;
            hAnglesLo = fMuonAnglesLo;
            hAnglesMi = fMuonAnglesMi;
            hAnglesHi = fMuonAnglesHi;
            hEnergy   = fMuonEnergy;
          }
          else if (std::abs(particle.PdgCode())==22) {
            hCosQ     = fPhotonCosQ;
            hAngles   = fPhotonAngles;
            hAnglesLo = fPhotonAnglesLo;
            hAnglesMi = fPhotonAnglesMi;
            hAnglesHi = fPhotonAnglesHi;
            hEnergy   = fPhotonEnergy;
          }
          else if (std::abs(particle.PdgCode())==11) {
            hCosQ     = fElectronCosQ;
            hAngles   = fElectronAngles;
            hAnglesLo = fElectronAnglesLo;
            hAnglesMi = fElectronAnglesMi;
            hAnglesHi = fElectronAnglesHi;
            hEnergy   = fElectronEnergy;
          }
          
          // TODO: check if the particle goes through the detector, if it does
          // add it to the list
          
          truth.Add(particle);
          
          if      (std::abs(particle.PdgCode())==13) ++numMuons;
          else if (std::abs(particle.PdgCode())==22) ++numPhotons;
          else if (std::abs(particle.PdgCode())==11) ++numElectrons;
          
          if (hCosQ!=0) {
            double cosq = -p4.Py()/p4.P();
            double phi  = std::atan2(p4.Pz(),p4.Px());
            phi *= 180/M_PI;
            hCosQ->Fill(cosq);
            hAngles->Fill(phi,cosq);
            if      (p4.E()<1.0)  hAnglesLo->Fill(phi,cosq);
            else if (p4.E()<10.0) hAnglesMi->Fill(phi,cosq);
            else                  hAnglesHi->Fill(phi,cosq);
            hEnergy->Fill(p4.E());
          }//end if there is a cos(theta) histogram
          
        }// loop on particles
        
	//}
      
      truthcol->push_back(truth);
      evt.put(std::move(truthcol));
      
      return;
    }// end produce
    
  }// end namespace
  
  
  namespace evgen{
    
    DEFINE_ART_MODULE(CosmicsGen)
    
  }
} // namespace gar

#endif // EVGEN_COSMICSGEN_H
////////////////////////////////////////////////////////////////////////
