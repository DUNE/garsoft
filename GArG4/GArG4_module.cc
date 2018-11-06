////////////////////////////////////////////////////////////////////////
/// \file  GArG4.cxx
/// \brief Use Geant4 to run the GArSoft detector simulation
///
/// \version $Id: GArG4.cxx,v 1.22 2010/07/20 06:08:30 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////

/// This a module.  It has the following functions:
///
/// - Initialize Geant4 physics, detector geometry, and other
///   processing.
///
/// - Accept sim::MCTruth objects from the MC branch of the FMWK
///   Event structure.
///
/// - Pass the primary particles to the Geant4 simulation to calculate
///   "truth" information for the detector response.
///
/// - Pass the truth information to the DetSim branch of the FMWK event.

#ifndef GARG4GARG4H
#define GARG4GARG4H 1

#include "nutools/G4Base/G4Helper.h"
#include "nutools/G4Base/ConvertMCTruthToG4.h"

// C++ Includes
#include <sstream> // std::ostringstream
#include <vector> // std::ostringstream
#include <map> // std::ostringstream
#include <set> // std::ostringstream
#include <iostream>
// #include <cstring>
#include <sys/stat.h>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

// nutools extensions
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/G4Base/DetectorConstruction.h"
#include "nutools/G4Base/UserActionManager.h"
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "GArG4/AuxDetAction.h"
#include "GArG4/EnergyDepositAction.h"
#include "GArG4/PhysicsList.h"
#include "GArG4/ParticleListAction.h"
#include "GArG4/MaterialPropertyLoader.h"
#include "GArG4/ParticleFilters.h" // garg4::PositionInVolumeFilter
#include "GArG4/G4SimulationParameters.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "SimulationDataProducts/LArDeposit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// G4 Includes
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4VUserDetectorConstruction.hh"
#include "Geant4/G4VUserPrimaryGeneratorAction.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4UserRunAction.hh"
#include "Geant4/G4UserEventAction.hh"
#include "Geant4/G4UserTrackingAction.hh"
#include "Geant4/G4UserSteppingAction.hh"
#include "Geant4/G4UserStackingAction.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/G4VUserPhysicsList.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4LogicalVolumeStore.hh"
#include "Geant4/G4RegionStore.hh"
#include "Geant4/Randomize.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/globals.hh"
#include "Geant4/G4ProductionCuts.hh"

// ROOT Includes
#include "TGeoManager.h"

// Forward declarations
class G4RunManager;
class G4UImanager;
class G4VisExecutive;

///Geant4 interface
namespace gar {
  namespace garg4 {

    // Forward declarations within namespace.
    class ParticleListAction;

    /**
     * @brief Runs Geant4 simulation and propagation of electrons and photons to readout
     *
     *
     * Randomness
     * -----------
     *
     * The random number generators used by this process are:
     * - 'GEANT' instance: used by Geant4
     * - 'propagation' instance: used in electron propagation
     * - 'radio' instance: used for radiological decay
     *
     *
     * Configuration parameters
     * -------------------------
     *
     * - <b>G4PhysListName</b> (string, default "garg4::PhysicsList"):
     *     whether to use the G4 overlap checker, which catches different issues than ROOT
     * - <b>CheckOverlaps</b> (bool, default false):
     *     whether to use the G4 overlap checker
     * - <b>DumpParticleList</b> (bool, default false):
     *     whether to print all MCParticles tracked
     * - <b>SmartStacking</b> (int, default 0):
     *     whether to use class to dictate how tracks are put on stack (nonzero is on)
     * - <b>KeepParticlesInVolumes</b> (vector<string>, default empty):
     *     list of volumes in which to keep MCParticles (empty keeps all)
     * - <b>GeantCommandFile</b> (string, required):
     *     G4 macro file to pass to G4Helper for setting G4 command
     * - <b>Seed</b> (pset key, not defined by default): if defined, override the seed for
     *     random number generator used in Geant4 simulation (which is obtained from
     *     NuRandomService by default)
     * - <b>PropagationSeed</b> (pset key, not defined by default): if defined,
     *     override the seed for the random generator used for electrons propagation
     *     to the wire planes (obtained from the NuRandomService by default)
     * - <b>RadioSeed</b> (pset key, not defined by default): if defined,
     *     override the seed for the random generator used for radiological decay
     *     (obtained from the NuRandomService by default)
     * - <b>InputLabels</b> (vector<string>, defualt unnecessary):
     *     optional list of generator labels which produce MCTruth;
     *     otherwise look for anything that has made MCTruth
     */
    class GArG4 : public ::art::EDProducer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit GArG4(fhicl::ParameterSet const& pset);
      virtual ~GArG4();

      /// The main routine of this module: Fetch the primary particles
      /// from the event, simulate their evolution in the detctor, and
      /// produce the detector response.
      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);

    private:
      g4b::G4Helper*              fG4Help;             ///< G4 interface object
      garg4::EnergyDepositAction* fEDepAction;         ///< Geant4 user action to handle GAr energy depositions
      garg4::AuxDetAction*        fAuxDetAction;       ///< Geant4 user action to handle GAr energy depositions
      garg4::ParticleListAction*  fParticleListAction; ///< Geant4 user action to particle information.
      fhicl::ParameterSet         fEDepActionPSet;     ///< configuration for GArAction
      fhicl::ParameterSet         fAuxDetActionPSet;   ///< configuration for AuxAction
      std::vector< std::string >  fSensDetVolumeName;   ///< Name of the sensitive volumes
      std::string                 fG4PhysListName;     ///< predefined physics list to use if not making a custom one
      std::string                 fG4MacroPath;        ///< directory path for Geant4 macro file to be
                                                       ///< executed before main MC processing.
      bool                        fCheckOverlaps;      ///< Whether to use the G4 overlap checker
      bool                        fdumpParticleList;   ///< Whether each event's sim::ParticleList will be displayed.
      int                         fSmartStacking;      ///< Whether to instantiate and use class to
                                                       ///< dictate how tracks are put on stack.
      float                       fMaxStepSize;        ///< maximum argon step size in cm
      float                       fProductionCut;      ///< G4 will check if a produced particle should travel this far and drop it if not
      std::vector<std::string>    fInputLabels;
      std::vector<std::string>    fKeepParticlesInVolumes; ///<Only write particles that have trajectories through these volumes

      /// Configures and returns a particle filter
      std::unique_ptr<PositionInVolumeFilter> CreateParticleVolumeFilter
      (std::set<std::string> const& vol_names) const;

      ///Set the Step limits of the volumes
      void SetVolumeStepLimit(G4LogicalVolumeStore *store, std::vector<std::string> SensDetVolumeName);

      ///Set the production cuts per region
      void SetProductionCuts();

    };

  } // namespace garg4

  namespace garg4 {

    //----------------------------------------------------------------------
    // Constructor
    GArG4::GArG4(fhicl::ParameterSet const& pset)
    : fG4Help                (nullptr)
    , fEDepAction            (nullptr)
    , fAuxDetAction          (nullptr)
    , fParticleListAction    (nullptr)
    , fEDepActionPSet        (pset.get<fhicl::ParameterSet>("EDepActionPSet")                         )
    , fAuxDetActionPSet      (pset.get<fhicl::ParameterSet>("AuxDetActionPSet")                       )
    , fSensDetVolumeName      (pset.get< std::vector< std::string > >("AuxDetVolumeName",{})           )
    , fG4PhysListName        (pset.get< std::string       >("G4PhysListName",   "garg4::PhysicsList") )
    , fCheckOverlaps         (pset.get< bool              >("CheckOverlaps",    false)                )
    , fdumpParticleList      (pset.get< bool              >("DumpParticleList", false)                )
    , fSmartStacking         (pset.get< int               >("SmartStacking",    0)                    )
    , fMaxStepSize           (pset.get< float             >("MaxStepSize",      0.2)                  )//in mm
    , fProductionCut         (pset.get< float             >("ProductionCut",    0.05)                 )//in mm
    , fKeepParticlesInVolumes(pset.get< std::vector< std::string > >("KeepParticlesInVolumes",{})     )
    {
      // initialize the GArSimulationParameters singleton
      G4SimulationParameters::CreateInstance(pset.get<fhicl::ParameterSet>("GArSimParsPSet"));

      LOG_DEBUG("GArG4") << "Debug: GArG4()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;

      if (pset.has_key("Seed")) {
        throw ::art::Exception(::art::errors::Configuration)
        << "The configuration of GArG4 module has the discontinued 'Seed' parameter.\n"
        "Seeds are now controlled by three parameters: 'GEANTSeed', 'PropagationSeed' and 'RadioSeed'.";
      }

      // setup the random number service for Geant4, the "G4Engine" label is a
      // special tag setting up a global engine for use by Geant4/CLHEP;
      // obtain the random seed from NuRandomService,
      // unless overridden in configuration with key "Seed" or "GEANTSeed"
      // same thing for the propagation engine:
      // and again for radio decay
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "G4Engine",       "GEANT",       pset, "GEANTSeed");
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "propagation", pset, "PropagationSeed");
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "radio",       pset, "RadioSeed");

      //get a list of generators to use, otherwise, we'll end up looking for anything that's
      //made an MCTruth object
      bool useInputLabels = pset.get_if_present< std::vector<std::string> >("InputLabels", fInputLabels);
      if(!useInputLabels) fInputLabels.resize(0);

      produces< std::vector<simb::MCParticle>                 >();
      produces< std::vector<sdp::EnergyDeposit>               >();
      produces< std::vector<sdp::CaloDeposit>                 >();
      produces< std::vector<sdp::LArDeposit>                 >();
      produces< ::art::Assns<simb::MCTruth, simb::MCParticle> >();

      // constructor decides if initialized value is a path or an environment variable
      cet::search_path sp("FW_SEARCH_PATH");

      sp.find_file(pset.get< std::string >("GeantCommandFile"), fG4MacroPath);
      struct stat sb;
      if (fG4MacroPath.empty() || stat(fG4MacroPath.c_str(), &sb)!=0)
        // failed to resolve the file name
        throw cet::exception("NoG4Macro")
        << "G4 macro file "
        << fG4MacroPath
        << " not found!\n";

    }

    //----------------------------------------------------------------------
    // Destructor
    GArG4::~GArG4()
    {
      if(fG4Help) delete fG4Help;
    }

    //----------------------------------------------------------------------
    void GArG4::SetVolumeStepLimit(G4LogicalVolumeStore *store, std::vector<std::string> SensDetVolumeName)
    {
        // Set the step size limits for the gaseous argon volume
        // get the logical volume for the desired volume name
        for(auto &name : SensDetVolumeName)
        {
            G4LogicalVolume *logVol_vec = store->GetVolume(name);

            if(logVol_vec != nullptr){

                LOG_INFO("GArG4")
                << "setting step limit size to be "
                << fMaxStepSize
                << " mm in volume "
                << name
                << " volume has address "
                << logVol_vec;

                fG4Help->SetVolumeStepLimit(name, fMaxStepSize);
            } else {
                LOG_DEBUG("GArG4")
                << "Volume "
                << name
                << " is null pointer.. "
                << "Verify the name of the sensitive volume";
            }
        }
    }

    //----------------------------------------------------------------------
    void GArG4::SetProductionCuts()
    {
        // Regions {TPC, Inner ECAL, Outer ECAL, Endcap ECAL}
        std::string regName[] = {"volGArTPC","InnerBarrelECal_vol", "OuterBarrelECal_vol", "InnerEndcapECal_vol"};

        for(G4int i = 0; i < 4; i++)
        {
            G4Region* aRegion = new G4Region(regName[i]);
            G4LogicalVolume* calorLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(regName[i]);
            calorLogical->SetRegion(aRegion);
            aRegion->AddRootLogicalVolume(calorLogical);

            G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName[i]);

            LOG_INFO("GArG4")
            << "Setting production cuts for region "
            << reg->GetName()
            << " production cut set to "
            << fProductionCut
            << " mm";

            G4ProductionCuts* cuts = new G4ProductionCuts();
            cuts->SetProductionCut(fProductionCut);
            reg->SetProductionCuts(cuts);
        }
    }

    //----------------------------------------------------------------------
    void GArG4::beginJob()
    {
      auto geo = gar::providerFrom<geo::Geometry>();
      auto* rng = &*(::art::ServiceHandle<::art::RandomNumberGenerator>());

      fG4Help = new g4b::G4Helper(fG4MacroPath, fG4PhysListName);
      if(fCheckOverlaps) fG4Help->SetOverlapCheck(true);
      fG4Help->ConstructDetector(geo->GDMLFile());

      G4LogicalVolumeStore *store = G4LogicalVolumeStore::GetInstance();
      if(store == nullptr) {
          LOG_INFO("GArG4")
          << "G4LogicalVolumeStore::GetInstance() not available";
      }

      // Get the logical volume store and assign material properties
      garg4::MaterialPropertyLoader* MPL = new garg4::MaterialPropertyLoader();
      MPL->GetPropertiesFromServices();
      MPL->UpdateGeometry(store);

      this->SetVolumeStepLimit(store, fSensDetVolumeName);

      this->SetProductionCuts();

      // Intialize G4 physics and primary generator action
      fG4Help->InitPhysics();

      // Use the UserActionManager to handle all the Geant4 user hooks.
      g4b::UserActionManager* uaManager = g4b::UserActionManager::Instance();

      // User-action class for accumulating particles and trajectories
      // produced in the detector.
      auto g4SimPars = gar::garg4::G4SimulationParameters::Instance();

      fParticleListAction = new garg4::ParticleListAction(g4SimPars->KineticEnergyCut(),
                                                          g4SimPars->StoreTrajectories(),
                                                          g4SimPars->KeepEMShowerDaughters());

      uaManager->AddAndAdoptAction(fParticleListAction);


      // add UserAction for handling steps in gaseous argon
      fEDepAction = new garg4::EnergyDepositAction(&rng->getEngine("propagation"),
                                                   fEDepActionPSet);
      uaManager->AddAndAdoptAction(fEDepAction);

      // add UserAction for handling steps in auxiliary detectors
      fAuxDetAction = new garg4::AuxDetAction(&rng->getEngine("propagation"),
                                              fAuxDetActionPSet);
      uaManager->AddAndAdoptAction(fAuxDetAction);


      // UserActionManager is now configured so continue G4 initialization
      fG4Help->SetUserAction();

      return;
    }

    //--------------------------------------------------------------------------
    void GArG4::beginRun(::art::Run& run)
    {
      // prepare the filter object (null if no filtering)

      std::set<std::string> volnameset(fKeepParticlesInVolumes.begin(), fKeepParticlesInVolumes.end());
      fParticleListAction->ParticleFilter(CreateParticleVolumeFilter(volnameset));

      return;
    }

    //--------------------------------------------------------------------------
    std::unique_ptr<PositionInVolumeFilter> GArG4::CreateParticleVolumeFilter(std::set<std::string> const& vol_names) const
    {

      // if we don't have favourite volumes, don't even bother creating a filter
      if (vol_names.empty()) return {};

      auto geom = providerFrom<geo::Geometry>();

      std::vector<std::vector<TGeoNode const*>> node_paths = geom->FindAllVolumePaths(vol_names);

      // collection of interesting volumes
      PositionInVolumeFilter::AllVolumeInfo_t GeoVolumePairs;
      GeoVolumePairs.reserve(node_paths.size()); // because we are obsessed

      //for each interesting volume, follow the node path and collect
      //total rotations and translations
      for(auto const& path : node_paths){

        TGeoTranslation* pTransl = new TGeoTranslation(0.,0.,0.);
        TGeoRotation* pRot = new TGeoRotation();
        for (TGeoNode const* node: path) {
          TGeoTranslation thistranslate(*node->GetMatrix());
          TGeoRotation thisrotate(*node->GetMatrix());
          pTransl->Add(&thistranslate);
          *pRot=*pRot * thisrotate;
        }

        //for some reason, pRot and pTransl don't have tr and rot bits set correctly
        //make new translations and rotations so bits are set correctly
        TGeoTranslation* pTransl2 = new TGeoTranslation(pTransl->GetTranslation()[0],
                                                        pTransl->GetTranslation()[1],
                                                        pTransl->GetTranslation()[2]);
        double phi=0.,theta=0.,psi=0.;
        pRot->GetAngles(phi,theta,psi);
        TGeoRotation* pRot2 = new TGeoRotation();
        pRot2->SetAngles(phi,theta,psi);

        TGeoCombiTrans* pTransf = new TGeoCombiTrans(*pTransl2,*pRot2);

        GeoVolumePairs.emplace_back(path.back()->GetVolume(), pTransf);

      }

      return std::make_unique<PositionInVolumeFilter>(std::move(GeoVolumePairs));

    } // CreateParticleVolumeFilter()


    //--------------------------------------------------------------------------
    void GArG4::produce(::art::Event& evt)
    {
      LOG_DEBUG("GArG4") << "produce()";

      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<simb::MCParticle> >                 partCol(new std::vector<simb::MCParticle>                );
      std::unique_ptr< std::vector<sdp::EnergyDeposit>  >              TPCCol  (new std::vector<sdp::EnergyDeposit>              );
      std::unique_ptr< ::art::Assns<simb::MCTruth, simb::MCParticle> > tpassn (new ::art::Assns<simb::MCTruth, simb::MCParticle>);
      std::unique_ptr< std::vector< sdp::CaloDeposit > >          ECALCol  (new std::vector<sdp::CaloDeposit>           );
      std::unique_ptr< std::vector< sdp::LArDeposit > >           LArCol  (new std::vector<sdp::LArDeposit>           );

      // reset the track ID offset as we have a new collection of interactions
      fParticleListAction->ResetTrackIDOffset();

      // look to see if there is any MCTruth information for this
      // event
      std::vector< ::art::Handle< std::vector<simb::MCTruth> > > mclists;
      if(fInputLabels.size() < 1)
        evt.getManyByType(mclists);
      else{
        mclists.resize(fInputLabels.size());
        for(size_t i = 0; i < fInputLabels.size(); ++i)
          evt.getByLabel(fInputLabels[i], mclists[i]);
      }

      // create a vector of art::Ptrs to MCTruths for use with the
      // association creation and a vector of const* for use with
      // G4Helper
      std::vector< ::art::Ptr<simb::MCTruth> > mctPtrs;
      std::vector< const simb::MCTruth*      > mcts;

      for(size_t mc = 0; mc < mclists.size(); ++mc){
        for(size_t i = 0; i < mclists[mc]->size(); ++i){
          art::Ptr<simb::MCTruth> ptr(mclists[mc], i);
          mctPtrs.push_back(ptr);
          mcts   .push_back(ptr.get());
        }
      }

      // Process Geant4 simulation for the mctruths
      fG4Help->G4Run(mcts);

      // receive the particle list
      auto       particleList = fParticleListAction->YieldList();
      auto const trackIDToMCT = fParticleListAction->TrackIDToMCTruthIndexMap();

      int    trackID             = std::numeric_limits<int>::max();
      size_t mctidx              = 0;
      size_t nGeneratedParticles = 0;

      // Has the user request a detailed dump of the output objects?
      if (fdumpParticleList){
        LOG_INFO("GArG4")
        << "Dump sim::ParticleList; size() = "
        << particleList.size()
        << "\n"
        << particleList;
      }

      auto iPartPair = particleList.begin();
      while (iPartPair != particleList.end()) {
        simb::MCParticle& p = *(iPartPair->second);

        trackID = p.TrackId();

        partCol->push_back(std::move(p));
        if( trackIDToMCT.count(trackID) > 0){

          mctidx = trackIDToMCT.find(trackID)->second;

          util::CreateAssn(*this,
                           evt,
                           *partCol,
                           mctPtrs[mctidx],
                           *tpassn,
                           nGeneratedParticles);
        }
        else
          throw cet::exception("GArG4")
          << "Cannot find MCTruth for Track Id: "
          << trackID
          << " to create association between Particle and MCTruth";

        // FIXME workaround until https://cdcvs.fnal.gov/redmine/issues/12067
        // is solved and adopted in GArSoft, after which moving will suffice
        // to avoid dramatic memory usage spikes;
        // for now, we immediately disposed of used particles
        iPartPair = particleList.erase(iPartPair);

        ++nGeneratedParticles;
      } // while(particleList)


      // Now for the sdp::EnergyDepositions
      for(auto const& tpc : fEDepAction->EnergyDeposits()){
        LOG_DEBUG("GArG4")
        << "adding TPC deposits for track id: "
        << tpc.TrackID();

        TPCCol->emplace_back(tpc);
      }

      // And finally the sdp::CaloDepositions
      for(auto const& hit : fAuxDetAction->CaloDeposits())
      {
        LOG_DEBUG("GArG4")
        << "adding calo deposits for track id: "
        << hit.TrackID();

        ECALCol->emplace_back(hit);
      }

      // And finally the sdp::CaloDepositions
      for(auto const& hit : fAuxDetAction->LArDeposits())
      {
        LOG_DEBUG("GArG4")
        << "adding LAr deposits for track id: "
        << hit.TrackID();

        LArCol->emplace_back(hit);
      }

      evt.put(std::move(TPCCol));
      evt.put(std::move(ECALCol));
      evt.put(std::move(LArCol));
      evt.put(std::move(partCol));
      evt.put(std::move(tpassn));

      return;
    } // GArG4::produce()

  } // namespace garg4

  namespace garg4 {

    DEFINE_ART_MODULE(GArG4)

  } // namespace garg4
} // gar
#endif // GARG4GARG4H
