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
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"
#include "art/Persistency/Common/PtrMaker.h"

// nutools extensions
#include "nutools/RandomUtils/NuRandomService.h"

// GArSoft Includes
#include "DetectorInfo/DetectorClocksService.h"
#include "Utilities/AssociationUtil.h"
#include "SimulationDataProducts/CaloDeposit.h"
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

            // Plugins should not be copied or assigned.
            CaloReadout(CaloReadout const &) = delete;
            CaloReadout(CaloReadout &&) = delete;
            CaloReadout & operator = (CaloReadout const &) = delete;
            CaloReadout & operator = (CaloReadout &&) = delete;

            void produce (::art::Event& evt) override;

            void reconfigure(fhicl::ParameterSet const& pset);

        private:

            void CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

            std::map<long long int, std::vector< art::Ptr<sdp::CaloDeposit> > > MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector);

            std::string                         fG4Label;    ///< label of G4 module
            const gar::geo::GeometryCore*       fGeo;        ///< geometry information
            std::unique_ptr<ECALReadoutSimAlg>  fROSimAlg;   ///< algorithm to simulate the electronics
        };

    } // namespace rosim

    namespace rosim {

        //----------------------------------------------------------------------
        // Constructor
        CaloReadout::CaloReadout(fhicl::ParameterSet const& pset) : art::EDProducer{pset}
        {
            fGeo = gar::providerFrom<geo::Geometry>();

            // setup the random number service
            // obtain the random seed from NuRandomService
            ::art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "sipm", pset, "SiPMSeed");

            this->reconfigure(pset);

            produces< std::vector<raw::CaloRawDigit> >();
            produces< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>  >();

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

            fG4Label = pset.get<std::string >("G4ModuleLabel", "geant");

            auto ECALROAlgPars = pset.get<fhicl::ParameterSet>("ECALReadoutSimAlgPars");
            auto ECALROAlgName = ECALROAlgPars.get<std::string>("ECALReadoutSimType");

            if(ECALROAlgName.compare("Standard") == 0)
            fROSimAlg = std::make_unique<gar::rosim::ECALReadoutSimStandardAlg>(rng->getEngine(art::ScheduleID::first(), pset.get<std::string>("module_label"), "sipm"), ECALROAlgPars);
            else
            throw cet::exception("CaloReadout")
            << "Unable to determine which ECAL readout simulation algorithm to use, bail";

            return;
        }

        //--------------------------------------------------------------------------
        void CaloReadout::produce(::art::Event& evt)
        {
            LOG_DEBUG("CaloReadout") << "produce()";

            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<sdp::CaloDeposit> > artHits;
            this->CollectHits(evt, fG4Label, artHits);

            //Get contributions per cellID for association to digi hit
            std::map<long long int, std::vector< art::Ptr<sdp::CaloDeposit> > > m_cIDMapArtPtrVec = this->MakeCellIDMapArtPtr(artHits);

            //Pass the sim hits to the algo
            fROSimAlg->PrepareAlgo(artHits);

            //Perform the digitization
            fROSimAlg->DoDigitization();

            //Get the digitized hits
            std::vector<raw::CaloRawDigit*> digiVec = fROSimAlg->GetDigitizedHits();

            // loop over the lists and put the particles and voxels into the event as collections
            std::unique_ptr< std::vector<raw::CaloRawDigit> > digitCol (new std::vector<raw::CaloRawDigit>);
            std::unique_ptr< art::Assns<raw::CaloRawDigit, sdp::CaloDeposit> > DigiSimHitsAssns(new art::Assns<raw::CaloRawDigit, sdp::CaloDeposit>);

            art::PtrMaker<raw::CaloRawDigit> makeDigiPtr(evt);

            for(auto it : digiVec)
            {
                raw::CaloRawDigit digihit = raw::CaloRawDigit(it->ADC(), it->Time(), it->X(), it->Y(), it->Z(), it->CellID());
                digitCol->push_back(digihit);

                art::Ptr<raw::CaloRawDigit> digiPtr = makeDigiPtr(digitCol->size() - 1);
                //get the associated vector of art ptr based on cellID
                std::vector< art::Ptr<sdp::CaloDeposit> > simPtrVec = m_cIDMapArtPtrVec[it->CellID()];
                for(auto hitpointer : simPtrVec)
                DigiSimHitsAssns->addSingle(digiPtr, hitpointer);
            }

            evt.put(std::move(digitCol));
            evt.put(std::move(DigiSimHitsAssns));

            return;
        } // CaloReadout::produce()

        //--------------------------------------------------------------------------
        void CaloReadout::CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
        {
            art::Handle< std::vector<sdp::CaloDeposit> > theHits;
            evt.getByLabel(label, theHits);

            if (!theHits.isValid())
            return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<sdp::CaloDeposit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        std::map<long long int, std::vector< art::Ptr<sdp::CaloDeposit> > > CaloReadout::MakeCellIDMapArtPtr(std::vector< art::Ptr<sdp::CaloDeposit> > &hitVector)
        {
            std::map<long long int, std::vector< art::Ptr<sdp::CaloDeposit> > > cIDMapArtPtrVec;

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

        DEFINE_ART_MODULE(CaloReadout)

    } // namespace rosim
} // gar
#endif // GAR_READOUTSIMULATION_CALOREADOUT
