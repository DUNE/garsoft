////////////////////////////////////////////////////////////////////////
/// \file  SiPMReadout_module.cc
/// \brief Create the readout from the sipm signals
///
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#ifndef GAR_READOUTSIMULATION_SIPMREADOUT
#define GAR_READOUTSIMULATION_SIPMREADOUT

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
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "art/Persistency/Common/PtrMaker.h"

// nutools extensions
#include "nurandom/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/CaloDeposit.h"
#include "ReadoutSimulation/ECALReadoutSimStandardAlg.h"
#include "ReadoutSimulation/MuIDReadoutSimStandardAlg.h"
#include "RawDataProducts/CaloRawDigit.h"
#include "Geometry/GeometryGAr.h"
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
    class SiPMReadout : public ::art::EDProducer{
    public:

      /// Standard constructor and destructor for an FMWK module.
      explicit SiPMReadout(fhicl::ParameterSet const& pset);
      virtual ~SiPMReadout();

      // Plugins should not be copied or assigned.
      SiPMReadout(SiPMReadout const &) = delete;
      SiPMReadout(SiPMReadout &&) = delete;
      SiPMReadout & operator = (SiPMReadout const &) = delete;
      SiPMReadout & operator = (SiPMReadout &&) = delete;

      void produce (::art::Event& evt) override;

      void reconfigure(fhicl::ParameterSet const& pset);

    private:

      void CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

      std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

      std::string                         fG4Label;    ///< label of G4 module
      std::string                         fInstanceLabelName; ///< product instance name

      const gar::geo::GeometryCore*       fGeo;        ///< geometry information
      std::unique_ptr<SiPMReadoutSimAlg>  fROSimAlg;   ///< algorithm to simulate the electronics

      CLHEP::HepRandomEngine              &fEngine;  ///< random engine
    };

  } // namespace rosim

  namespace rosim {

    //----------------------------------------------------------------------
    // Constructor
    SiPMReadout::SiPMReadout(fhicl::ParameterSet const& pset)
      : art::EDProducer{pset},
      fEngine(art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this,pset,"Seed"))
      {
        fGeo = gar::providerFrom<geo::GeometryGAr>();

        this->reconfigure(pset);

        art::InputTag tag(fG4Label, fInstanceLabelName);
        consumes< std::vector<sdp::CaloDeposit> >(tag);
        produces< std::vector<raw::CaloRawDigit> >(fInstanceLabelName);
        produces< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>  >(fInstanceLabelName);

        return;
      }

    //----------------------------------------------------------------------
    // Destructor
    SiPMReadout::~SiPMReadout()
    {
    }

    //----------------------------------------------------------------------
    void SiPMReadout::reconfigure(fhicl::ParameterSet const& pset)
    {
      MF_LOG_DEBUG("SiPMReadout") << "Debug: SiPMReadout()";

      fG4Label = pset.get<std::string >("G4ModuleLabel", "geant");
      fInstanceLabelName =  pset.get<std::string >("InstanceLabelName", "");

      auto ROAlgPars = pset.get<fhicl::ParameterSet>("ReadoutSimAlgPars");
      auto ROAlgName = ROAlgPars.get<std::string>("ReadoutSimType");

      if(ROAlgName.compare("Standard_ECAL") == 0) {
        fROSimAlg = std::make_unique<gar::rosim::ECALReadoutSimStandardAlg>(fEngine, ROAlgPars);
      } 
      else if(ROAlgName.compare("Standard_MuID") == 0) {
        fROSimAlg = std::make_unique<gar::rosim::MuIDReadoutSimStandardAlg>(fEngine, ROAlgPars);
      }
      else {
        throw cet::exception("SiPMReadout")
          << "Unable to determine which readout simulation algorithm to use, bail";
      }

      return;
    }

    //--------------------------------------------------------------------------
    void SiPMReadout::produce(::art::Event& evt)
    {
      MF_LOG_DEBUG("SiPMReadout") << "produce()";

      //Collect the hits to be passed to the algo
      std::vector< art::Ptr<sdp::CaloDeposit> > artHits;
      this->CollectHits(evt, fG4Label, fInstanceLabelName, artHits);

      //Get contributions per cellID for association to digi hit
      std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > m_cIDMapArtPtrVec = this->MakeCellIDMapArtPtr(artHits);

      //Pass the sim hits to the algo
      fROSimAlg->PrepareAlgo(artHits);

      //Perform the digitization
      fROSimAlg->DoDigitization();

      //Get the digitized hits
      std::vector< raw::CaloRawDigit* > digiVec = fROSimAlg->GetDigitizedHits();

      // loop over the lists and put the particles and voxels into the event as collections
      std::unique_ptr< std::vector<raw::CaloRawDigit> > digitCol (new std::vector<raw::CaloRawDigit>);
      std::unique_ptr< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit> > DigiSimHitsAssns(new art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>);

      art::PtrMaker<raw::CaloRawDigit> makeDigiPtr(evt, fInstanceLabelName);

      for(auto const &it : digiVec)
        {
          digitCol->emplace_back(*it);
          art::Ptr<raw::CaloRawDigit> digiPtr = makeDigiPtr(digitCol->size() - 1);
          //get the associated vector of art ptr based on cellID
          std::vector< art::Ptr<sdp::CaloDeposit> > simPtrVec = m_cIDMapArtPtrVec[it->CellID()];
          for(auto hitpointer : simPtrVec)
            DigiSimHitsAssns->addSingle(digiPtr, hitpointer);
        }

      evt.put(std::move(digitCol), fInstanceLabelName);
      evt.put(std::move(DigiSimHitsAssns), fInstanceLabelName);

      return;
    } // CaloReadout::produce()

    //--------------------------------------------------------------------------
    void SiPMReadout::CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
    {
      art::InputTag itag(label, instance);
      auto theHits = evt.getHandle< std::vector<sdp::CaloDeposit> >(itag);
      if (!theHits)
        return;

      for (unsigned int i = 0; i < theHits->size(); ++i)
        {
          const art::Ptr<sdp::CaloDeposit> hit(theHits, i);
          hitVector.push_back(hit);
        }
    }

    //--------------------------------------------------------------------------
    std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > SiPMReadout::MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
    {
      std::map<raw::CellID_t, std::vector< art::Ptr<sdp::CaloDeposit> > > cIDMapArtPtrVec;

      for (std::vector< art::Ptr<sdp::CaloDeposit> >::const_iterator iter = hitVector.begin(), iterEnd = hitVector.end(); iter != iterEnd; ++iter)
        {
          art::Ptr<sdp::CaloDeposit> hitPtr = *iter;
          const sdp::CaloDeposit *hit = hitPtr.get();

          if( cIDMapArtPtrVec.count(hit->CellID()) == 0 )
            {
              std::vector< art::Ptr<sdp::CaloDeposit> > vecArtPtr;
              vecArtPtr.push_back(hitPtr);
              cIDMapArtPtrVec.insert( std::make_pair(hit->CellID(), vecArtPtr) );
            }
          else
            {
              cIDMapArtPtrVec[hit->CellID()].push_back(hitPtr);
            }
        }

      return cIDMapArtPtrVec;
    }


  } // namespace rosim

  namespace rosim {

    DEFINE_ART_MODULE(SiPMReadout)

  } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_SIPMREADOUT
