////////////////////////////////////////////////////////////////////////
// Class:       CaloHitFinder
// Plugin Type: producer (art v2_11_02)
// File:        CompressedHitFinder_module.cc
//
// Generated at Wed Jun 13 15:27:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Report the compressed raw digit blocks as hits -- just a threshold
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "Utilities/ECALUtils.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include <memory>

namespace gar {
    namespace rec {

        class CaloHitFinder : public art::EDProducer {
        public:
            explicit CaloHitFinder(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            CaloHitFinder(CaloHitFinder const &) = delete;
            CaloHitFinder(CaloHitFinder &&) = delete;
            CaloHitFinder & operator = (CaloHitFinder const &) = delete;
            CaloHitFinder & operator = (CaloHitFinder &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

            float CalibrateToMIP(unsigned int ADC);

            float CalibratetoMeV(float &x, float &y, float &z, double MIP, const long long int& cID);

        private:

            // Declare member data here.

            float fMIPThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
            float fMIPtoMeV; ///< in MeV / MIP for 5 mm
            bool fDesaturation; ///< flag to perform the SiPM desaturation

            std::string fRawDigitLabel;  ///< label to find the right raw digits

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
            std::unique_ptr<util::ECALUtils>          fECALUtils;    ///< pointer to the Util fcn for the ECAL containing the desaturation function
        };


        CaloHitFinder::CaloHitFinder(fhicl::ParameterSet const & p)
        // :
        {
            fMIPThreshold = p.get<float>("MIPThreshold", 0.25);
            fRawDigitLabel = p.get<std::string>("RawDigitLabel", "daqecal");
            fMIPtoMeV     = p.get<float>("MIPtoMeV", 0.814);
            fDesaturation = p.get<bool>("Desaturation", false);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fECALUtils = std::make_unique<util::ECALUtils>(fDetProp->EffectivePixel(), 0.95);

            produces< std::vector<rec::CaloHit> >();
            produces< art::Assns<rec::CaloHit, raw::CaloRawDigit>  >();
        }

        void CaloHitFinder::produce(art::Event & e)
        {
            // create an emtpy output hit collection -- to add to.
            std::unique_ptr< std::vector<rec::CaloHit> > hitCol ( new std::vector< rec::CaloHit > );
            std::unique_ptr< art::Assns<rec::CaloHit, raw::CaloRawDigit> > RecoDigiHitsAssns( new art::Assns<rec::CaloHit, raw::CaloRawDigit> );

            art::PtrMaker<rec::CaloHit> makeCaloHitPtr(e);

            // the input raw digits
            auto digiCol = e.getValidHandle< std::vector<raw::CaloRawDigit> >(fRawDigitLabel);

            for (size_t idigit = 0; idigit < digiCol->size(); ++idigit)
            {
                const raw::CaloRawDigit& digitHit = (*digiCol)[idigit];

                unsigned int hitADC = digitHit.ADC();
                float hitTime = digitHit.Time();
                float x = digitHit.X();
                float y = digitHit.Y();
                float z = digitHit.Z();
                long long int cellID = digitHit.CellID();

                //Do Calibration of the hit in MIPs
                float hitMIP = this->CalibrateToMIP(hitADC);

                if(hitMIP < fMIPThreshold)
                {
                    LOG_DEBUG("CaloHitFinder") << "Signal under the " << fMIPThreshold << " MIP threshold" << std::endl;
                    continue;
                }

                float energy = 0.;

                if(fDesaturation)
                {
                    double sat_px = hitMIP * fDetProp->LightYield();
                    //DeSaturate
                    double unsat_px = fECALUtils->DeSaturate(sat_px);

                    //Calibrate to the MeV scale
                    double unsat_energy = unsat_px / fDetProp->LightYield();
                    energy = this->CalibratetoMeV(x, y, z, unsat_energy, cellID);
                }
                else{
                    //Calibrate to the MeV scale
                    energy = this->CalibratetoMeV(x, y, z, hitMIP, cellID);
                }

                float pos[3] = {x, y, z};

                //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                rec::CaloHit hit(energy * CLHEP::MeV / CLHEP::GeV, hitTime, pos, cellID);
                hitCol->emplace_back(hit);

                //Make association between digi hits and reco hits
                art::Ptr<rec::CaloHit> hitPtr = makeCaloHitPtr(hitCol->size() - 1);
                art::Ptr<raw::CaloRawDigit> digiArtPtr = art::Ptr<raw::CaloRawDigit>(digiCol, idigit);
                RecoDigiHitsAssns->addSingle(hitPtr, digiArtPtr);
            }

            //move the reco hit collection
            e.put(std::move(hitCol));
            e.put(std::move(RecoDigiHitsAssns));
            return;
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CalibrateToMIP(unsigned int ADC)
        {
            if(ADC <= 0) return ADC;

            return ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CalibratetoMeV(float &x, float &y, float &z, double MIP, const long long int& cID)
        {
            if(MIP <= 0) return MIP;

            TVector3 point(x, y, z);
            float factor = fGeo->GetSensVolumeThickness(point) / 0.5;

            return MIP * fMIPtoMeV * factor; // MeV
        }

        DEFINE_ART_MODULE(CaloHitFinder)

    } // namespace rec
} // namespace gar
