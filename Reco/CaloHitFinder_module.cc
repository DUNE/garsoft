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
#include "Geometry/LocalTransformation.h"

#include "Utilities/ECALUtils.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TGeoNode.h"
#include "TGeoManager.h"

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

        protected:
            float CalibrateToMIP(unsigned int ADC);

            float CalibratetoMeV(float x, float y, float z, double MIP);

            std::array<double, 3U> CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, long long int cID);

            float CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, long long int cID);

        private:

            // Declare member data here.
            float fMIPThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
            bool fDesaturation; ///< flag to perform the SiPM desaturation
            std::string fRawDigitLabel;  ///< label to find the right raw digits

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
            std::unique_ptr<util::ECALUtils>          fECALUtils;    ///< pointer to the Util fcn for the ECAL containing the desaturation function
            TGeoManager* fGeoManager;
        };


        CaloHitFinder::CaloHitFinder(fhicl::ParameterSet const & p)
        // :
        {
            fMIPThreshold = p.get<float>("MIPThreshold", 0.25);
            fRawDigitLabel = p.get<std::string>("RawDigitLabel", "daqecal");
            fDesaturation = p.get<bool>("Desaturation", false);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fGeoManager = fGeo->ROOTGeoManager();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fECALUtils = std::make_unique<util::ECALUtils>(fDetProp->EffectivePixel());

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
                std::pair<float, float> hitTime = digitHit.Time();
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

                //Desaturation
                float energy = 0.;
                if(fDesaturation)
                {
                    double sat_px = hitMIP * fDetProp->LightYield();
                    //DeSaturate
                    double unsat_px = fECALUtils->DeSaturate(sat_px);

                    //Calibrate to the MeV scale
                    double unsat_energy = unsat_px / fDetProp->LightYield();
                    energy = this->CalibratetoMeV(x, y, z, unsat_energy);
                }
                else{
                    //Calibrate to the MeV scale
                    energy = this->CalibratetoMeV(x, y, z, hitMIP);
                }

                //Position reconstruction based on time for the strips
                float pos[3] = {0., 0., 0.};
                float time = 0.;

                if(fGeo->isTile(cellID))
                {
                    pos[0] = x;
                    pos[1] = y;
                    pos[2] = z;
                    time = hitTime.first;
                }
                else{
                    std::array<double, 3U> strip_pos = this->CalculateStripHitPosition(x, y, z, hitTime, cellID);
                    pos[0] = strip_pos[0];
                    pos[1] = strip_pos[1];
                    pos[2] = strip_pos[2];

                    //Correct the time based on the strip length
                    time = this->CorrectStripHitTime(x, y, z, hitTime, cellID);
                }

                //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                rec::CaloHit hit(energy * CLHEP::MeV / CLHEP::GeV, time, pos, cellID);
                hitCol->emplace_back(hit);

                LOG_DEBUG("CaloHitFinder") << "recohit " << &hit
                << " with cellID " << cellID
                << " has energy " << energy * CLHEP::MeV / CLHEP::GeV
                << " time " << time << " ns"
                << " pos (" << pos[0] << ", " <<  pos[1] << ", " << pos[2] << ")";

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
            float calibrated_ADC = 0.;
            if(ADC <= 0) {
                return calibrated_ADC;
            }

            calibrated_ADC = ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
            return calibrated_ADC;
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CalibratetoMeV(float x, float y, float z, double MIP)
        {
            float calibrated_MIP = 0.;
            if(MIP <= 0) {
                return calibrated_MIP;
            }

            TVector3 point(x, y, z);
            float factor = fGeo->GetSensVolumeThickness(point) / 0.5;

            calibrated_MIP = MIP * fDetProp->MeVtoMIP() * factor;
            return calibrated_MIP; // MeV
        }

        //----------------------------------------------------------------------------
        std::array<double, 3U> CaloHitFinder::CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, long long int cID)
        {
            // TGeoNode *node = fGeoManager->FindNode(x, y, z);//Node in cm...
            // double stripLength = fGeo->getStripLength(node, cID); // in mm

            //Find the volume path
            TVector3 point(x, y, z);
            std::string name = fGeo->VolumeName(point);
            auto const& path = fGeo->FindVolumePath(name);
            if (path.empty())
            {
                throw cet::exception("CaloHitFinder") << "CalculateStripHitPosition(): can't find volume '" << name << "'\n";
            }

            //Change to local frame
            gar::geo::LocalTransformation<TGeoHMatrix> trans(path, path.size() - 1);
            std::array<double, 3U> world{ {x, y, z} }, local;
            trans.WorldToLocal(world.data(), local.data());

            //Calculate the position of the hit based on the time of both SiPM along the strip
            // pos along strip is
            // x = c * (t1 - t2) / 2 (0 is the center of the strip)
            float c = CLHEP::c_light * CLHEP::mm / CLHEP::ns;
            float xlocal = c * ( hitTime.first - hitTime.second ) / 2.;

            if( hitTime.first - hitTime.second < 0)
            xlocal = - xlocal;

            std::array<double, 3U> local_back = fGeo->ReconstructStripHitPosition(local, xlocal / CLHEP::cm, cID);
            std::array<double, 3U> world_back;

            trans.LocalToWorld(local_back.data(), world_back.data());

            return world_back;
        }

        //----------------------------------------------------------------------------
        float CaloHitFinder::CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, long long int cID)
        {
            TGeoNode *node = fGeoManager->FindNode(x, y, z);//Node in cm...
            double stripLength = fGeo->getStripLength(node, cID); // in mm

            float c = CLHEP::c_light * CLHEP::mm / CLHEP::ns;
            float time = (hitTime.first + hitTime.second) / 2. - (stripLength / (2 * c));

            return time;
        }

        DEFINE_ART_MODULE(CaloHitFinder)

    } // namespace rec
} // namespace gar
