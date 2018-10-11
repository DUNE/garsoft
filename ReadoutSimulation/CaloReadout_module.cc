////////////////////////////////////////////////////////////////////////
/// \file  CaloReadout_module.cc
/// \brief Create the readout from the calorimeter signals
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_CALOREADOUT
#define GAR_READOUTSIMULATION_CALOREADOUT

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

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
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "DetectorInfo/DetectorClocksService.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReadoutSimulation/IonizationAndScintillation.h"
#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "Geometry/Geometry.h"
#include "CoreUtils/ServiceUtil.h"

// Forward declarations

///Geant4 interface
namespace gar {
  namespace rosim {

    /**
    * @brief Runs simulation including smearing of energy and SiPM saturation for the ECAL readout
    *
    *
    * Randomness
    * -----------
    *
    * The random number generators used by this process are:
    * - 'GEANT' instance: used by Geant4
    * - 'propagation' instance: used in electron propagation
    *
    */
    class CaloReadout : public ::art::EDProducer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit CaloReadout(fhicl::ParameterSet const& pset);
      virtual ~CaloReadout();

      void produce (::art::Event& evt);
      void beginJob();
      void beginRun(::art::Run& run);
      void reconfigure(fhicl::ParameterSet const& pset);

    private:

      void CreateSignalDigit(std::vector<raw::CaloRawDigit> & digCol, ::art::ValidHandle< std::vector<sdp::CaloDeposit> > & eDepCol, ::art::Event & evt);

      std::string                         fG4Label;    ///< label of G4 module
      const gar::geo::GeometryCore*       fGeo;        ///< geometry information
      std::unique_ptr<ECALReadoutSimAlg>  fROSimAlg;   ///< algorithm to simulate the electronics
    };

  } // namespace rosim

  namespace rosim {

    //----------------------------------------------------------------------
    // Constructor
    CaloReadout::CaloReadout(fhicl::ParameterSet const& pset)
    {
      fGeo = gar::providerFrom<geo::Geometry>();

      // setup the random number service
      // obtain the random seed from NuRandomService
      ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "sipm", pset, "SiPMSeed");

      this->reconfigure(pset);

      produces< std::vector<raw::CaloRawDigit>                      >();
      produces< ::art::Assns<sdp::CaloDeposit, raw::CaloRawDigit>   >();

      return;
    }

    //----------------------------------------------------------------------
    // Destructor
    CaloReadout::~CaloReadout()
    {
    }

    //----------------------------------------------------------------------
    void CaloReadout::reconfigure(fhicl::ParameterSet const& pset)
    {
      LOG_DEBUG("CaloReadout") << "Debug: CaloReadout()";
      ::art::ServiceHandle<::art::RandomNumberGenerator> rng;

      fG4Label    = pset.get<std::string        >("G4ModuleLabel",      "geant");

      auto ECALROAlgPars = pset.get<fhicl::ParameterSet>("ECALReadoutSimAlgPars");
      auto ECALROAlgName = ECALROAlgPars.get<std::string>("ECALReadoutSimType");

      if(ECALROAlgName.compare("Standard") == 0)
      fROSimAlg = std::make_unique<gar::rosim::ECALReadoutSimStandardAlg>(rng->getEngine("sipm"),
      ECALROAlgPars);
      else
      throw cet::exception("CaloReadout")
      << "Unable to determine which ECAL readout simulation algorithm to use, bail";

      return;
    }

    //----------------------------------------------------------------------
    void CaloReadout::beginJob()
    {
      return;
    }

    //--------------------------------------------------------------------------
    void CaloReadout::beginRun(::art::Run& run)
    {
      return;
    }

    //--------------------------------------------------------------------------
    void CaloReadout::produce(::art::Event& evt)
    {
      LOG_DEBUG("CaloReadout") << "produce()";

      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<raw::CaloRawDigit>                    > rdCol (new std::vector<raw::CaloRawDigit>                   );
      // std::unique_ptr< ::art::Assns<sdp::CaloDeposit, raw::CaloRawDigit> > erassn(new ::art::Assns<sdp::CaloDeposit, raw::CaloRawDigit>);

      auto eDepCol = evt.getValidHandle< std::vector<sdp::CaloDeposit> >(fG4Label);

      if(eDepCol->size() > 0){
        //Create the digitized signal for the Calo hits
        this->CreateSignalDigit(*rdCol, eDepCol, evt);

        // for(size_t ihit = 0; eDepCol->size(); ++ihit)
        // {
        //   art::Ptr<sdp::CaloDeposit> CaloPtr(eDepCol, ihit);
        //   util::CreateAssn(*this, evt, *rdCol, CaloPtr, *erassn);
        // }
      }

      evt.put(std::move(rdCol));
      // evt.put(std::move(erassn));

      return;
    } // CaloReadout::produce()

    //--------------------------------------------------------------------------
    void CaloReadout::CreateSignalDigit(std::vector<raw::CaloRawDigit> & digCol, ::art::ValidHandle< std::vector<sdp::CaloDeposit> > & eDepCol, ::art::Event & evt)
    {
      std::vector<sdp::CaloDeposit> const& CaloVec(*eDepCol);

      fROSimAlg->CreateCaloRawDigits(CaloVec, digCol);
      fROSimAlg->CreateNoiseDigits(digCol);

      return;
    }

  } // namespace rosim

  namespace rosim {

    DEFINE_ART_MODULE(CaloReadout)

  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_CALOREADOUT
