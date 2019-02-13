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

#include "RawDataProducts/CaloRawDigit.h"
#include "RawDataProducts/raw.h"
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
            float CalibratetoMeV(double MIP, const long long int& cID);

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
        }

        void CaloHitFinder::produce(art::Event & e)
        {
            // create an emtpy output hit collection -- to add to.
            std::unique_ptr<std::vector<CaloHit> > hitCol (new std::vector<CaloHit> );

            // the input raw digits
            auto rdCol = e.getValidHandle< std::vector<raw::CaloRawDigit> >(fRawDigitLabel);

            for (size_t ird = 0; ird < rdCol->size(); ++ ird)
            {
                auto const& rd = (*rdCol)[ird];
                unsigned int hitADC = rd.ADC();
                float hitTime = rd.Time();
                float x = rd.X();
                float y = rd.Y();
                float z = rd.Z();
                long long int cellID = rd.CellID();

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
                    energy = this->CalibratetoMeV(unsat_energy, cellID);
                }
                else{
                    //Calibrate to the MeV scale
                    energy = this->CalibratetoMeV(hitMIP, cellID);
                }

                float pos[3] = {x, y, z};

                //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                hitCol->emplace_back(energy * CLHEP::MeV / CLHEP::GeV, hitTime, pos, cellID);
            }

            //move the reco hit collection
            e.put(std::move(hitCol));
            return;
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CalibrateToMIP(unsigned int ADC)
        {
            if(ADC <= 0) return ADC;

            return ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CalibratetoMeV(double MIP, const long long int& cID)
        {
            if(MIP <= 0) return MIP;

            if(fGeo->isTile(cID)){
                return MIP * fMIPtoMeV * 2; // MeV
            }
            else{
                return MIP * fMIPtoMeV; // MeV
            }
        }

        DEFINE_ART_MODULE(CaloHitFinder)

    } // namespace rec
} // namespace gar
