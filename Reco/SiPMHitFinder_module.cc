////////////////////////////////////////////////////////////////////////
// Class:       SiPMHitFinder
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
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "RawDataProducts/CaloRawDigit.h"
#include "ReconstructionDataProducts/CaloHit.h"
#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Geometry/LocalTransformation.h"
#include "Geometry/BitFieldCoder.h"

#include "RecoAlg/SiPMUtils.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "TGeoNode.h"
#include "TGeoManager.h"

#include <memory>

namespace gar {
    namespace rec {

        class SiPMHitFinder : public art::EDProducer {
        public:
            explicit SiPMHitFinder(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            SiPMHitFinder(SiPMHitFinder const &) = delete;
            SiPMHitFinder(SiPMHitFinder &&) = delete;
            SiPMHitFinder & operator = (SiPMHitFinder const &) = delete;
            SiPMHitFinder & operator = (SiPMHitFinder &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        protected:
            void CollectDigiHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::raw::CaloRawDigit> > &hitVector);

            float CalibrateToMIP(unsigned int ADC);

            float CalibratetoMeV(float x, float y, float z, double MIP);

            std::array<double, 3U> CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID);

            float CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID);

        private:

            // Declare member data here.
            std::string fRawDigitLabel;  ///< label to find the right raw digits
            std::string fInstanceLabelName; ///< product instance name

            float fMIPThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
            bool  fDesaturation; ///< flag to perform the SiPM desaturation
            float fSamplingCorrection;
            bool  fUseTimePositionReco;

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
            std::unique_ptr<util::SiPMUtils>    fSiPMUtils;    ///< pointer to the Util fcn containing the SiPM desaturation function
            gar::geo::BitFieldCoder const* fFieldDecoder;
        };


        SiPMHitFinder::SiPMHitFinder(fhicl::ParameterSet const & p) : EDProducer{p}
        {
            fRawDigitLabel       = p.get<std::string>("RawDigitLabel", "daqsipm");
            fInstanceLabelName   = p.get<std::string>("InstanceLabelName", "");

            fMIPThreshold        = p.get<float>("MIPThreshold", 0.25);
            fDesaturation        = p.get<bool>("Desaturation", false);
            fSamplingCorrection  = p.get<float>("ECALSamplingFactorGeV", 1.0);
            fUseTimePositionReco = p.get<bool>("UseTimePositionReco", true);

            fGeo       = gar::providerFrom<geo::GeometryGAr>();
            fDetProp   = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fSiPMUtils = std::make_unique<util::SiPMUtils>(fDetProp->EffectivePixel());

            std::string fEncoding = fGeo->GetECALCellIDEncoding();
            art::InputTag tag(fRawDigitLabel, fInstanceLabelName);
            consumes< std::vector<raw::CaloRawDigit> >(tag);
            produces< std::vector<rec::CaloHit> >(fInstanceLabelName);
            produces< art::Assns<rec::CaloHit, raw::CaloRawDigit>  >(fInstanceLabelName);

            if(fInstanceLabelName.compare("MuID") == 0) {
                fEncoding = fGeo->GetMuIDCellIDEncoding();
                fSamplingCorrection = p.get<float>("MuIDSamplingFactorGeV", 1.0);
            }

            fFieldDecoder = new gar::geo::BitFieldCoder( fEncoding );
        }

        void SiPMHitFinder::produce(art::Event & e)
        {
            // create an emtpy output hit collection -- to add to.
            std::unique_ptr< std::vector<rec::CaloHit> > hitCol ( new std::vector< rec::CaloHit > );
            std::unique_ptr< art::Assns<rec::CaloHit, raw::CaloRawDigit> > RecoDigiHitsAssns( new art::Assns<rec::CaloHit, raw::CaloRawDigit> );

            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<gar::raw::CaloRawDigit> > artDigiHits;
            this->CollectDigiHits(e, fRawDigitLabel, fInstanceLabelName, artDigiHits);

            art::PtrMaker<gar::rec::CaloHit> makeHitPtr(e, fInstanceLabelName);

            for (std::vector< art::Ptr<gar::raw::CaloRawDigit> >::const_iterator iter = artDigiHits.begin(), iterEnd = artDigiHits.end(); iter != iterEnd; ++iter)
            {
                art::Ptr<gar::raw::CaloRawDigit> digiArtPtr = *iter;
                const gar::raw::CaloRawDigit &digitHit = *(digiArtPtr.get());

                unsigned int hitADC = digitHit.ADC().first;
                std::pair<float, float> hitTime = digitHit.Time();
                float x = digitHit.X();
                float y = digitHit.Y();
                float z = digitHit.Z();
                raw::CellID_t cellID = digitHit.CellID();

                //Do Calibration of the hit in MIPs
                float hitMIP = this->CalibrateToMIP(hitADC);

                if (hitMIP < fMIPThreshold)
                {
                    MF_LOG_DEBUG("SiPMHitFinder") << "Signal under the " << fMIPThreshold << " MIP threshold" << std::endl;
                    continue;
                }

                //Desaturation
                float energy = 0.;
                if (fDesaturation)
                {
                    double sat_px = hitMIP * fDetProp->LightYield();
                    //DeSaturate
                    double unsat_px = fSiPMUtils->DeSaturate(sat_px);

                    //Calibrate to the MeV scale
                    double unsat_energy = unsat_px / fDetProp->LightYield();
                    energy = this->CalibratetoMeV(x, y, z, unsat_energy);
                }
                else {
                    //Calibrate to the MeV scale
                    energy = this->CalibratetoMeV(x, y, z, hitMIP);
                }

                //Position reconstruction based on time for the strips
                float pos[3] = {0., 0., 0.};
                std::pair<float, float> time;
                const std::array<double, 3> point = { x, y, z };

                if (fGeo->isTile(point, cellID))
                {
                    pos[0] = x;
                    pos[1] = y;
                    pos[2] = z;
                    time = hitTime;
                }
                else {
                    if( fUseTimePositionReco) {
                        std::array<double, 3> strip_pos = this->CalculateStripHitPosition(x, y, z, hitTime, cellID);
                        pos[0] = strip_pos[0];
                        pos[1] = strip_pos[1];
                        pos[2] = strip_pos[2];
                        time = std::make_pair( this->CorrectStripHitTime(pos[0], pos[1], pos[2], hitTime, cellID), 0. );
                    } else {
                        pos[0] = x;
                        pos[1] = y;
                        pos[2] = z;
                        time = hitTime;
                    }
                }

                //Store the hit (energy in GeV, time in ns, pos in cm and cellID)
                unsigned int layer = fFieldDecoder->get(cellID, "layer");
                rec::CaloHit hit(fSamplingCorrection * energy * CLHEP::MeV / CLHEP::GeV, time, pos, cellID, layer);
                hitCol->emplace_back(hit);

                MF_LOG_DEBUG("SiPMHitFinder") << "recohit " << &hit
                << " with cellID " << cellID
                << " has energy " << fSamplingCorrection * energy * CLHEP::MeV / CLHEP::GeV
                << " pos (" << pos[0] << ", " <<  pos[1] << ", " << pos[2] << ")";

                //Make association between digi hits and reco hits
                art::Ptr<rec::CaloHit> hitPtr = makeHitPtr(hitCol->size() - 1);
                RecoDigiHitsAssns->addSingle(hitPtr, digiArtPtr);
            }

            //move the reco hit collection
            e.put(std::move(hitCol), fInstanceLabelName);
            e.put(std::move(RecoDigiHitsAssns), fInstanceLabelName);

            return;
        }

        //----------------------------------------------------------------------------
        void SiPMHitFinder::CollectDigiHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::raw::CaloRawDigit> > &hitVector)
        {
	    art::InputTag itag(label,instance);
	    auto theHits = evt.getHandle< std::vector<gar::raw::CaloRawDigit> >(itag);
            if (!theHits) return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::raw::CaloRawDigit> hit(theHits, i);
                hitVector.push_back(hit);
            }
            return;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CalibrateToMIP(unsigned int ADC)
        {
            float calibrated_ADC = 0.;
            if (ADC <= 0) {
                return calibrated_ADC;
            }

            calibrated_ADC = ADC / (fDetProp->LightYield() * fDetProp->SiPMGain());
            return calibrated_ADC;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CalibratetoMeV(float x, float y, float z, double MIP)
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
        std::array<double, 3> SiPMHitFinder::CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            std::array<double, 3> pointLocal;
            gar::geo::LocalTransformation<TGeoHMatrix> trans;
            fGeo->WorldToLocal(point, pointLocal, trans);

            //Calculate the position of the hit based on the time of both SiPM along the strip
            // pos along strip is
            // x = c * (t1 - t2) / 2 (0 is the center of the strip)
            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float xlocal = c * ( hitTime.first - hitTime.second ) / 2.;

            std::array<double, 3> local_back = fGeo->ReconstructStripHitPosition(point, pointLocal, xlocal, cID);
            std::array<double, 3> world_back;
            fGeo->LocalToWorld(local_back, world_back, trans);

            return world_back;
        }

        //----------------------------------------------------------------------------
        float SiPMHitFinder::CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            double stripLength = fGeo->getStripLength(point, cID); // in cm

            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float time = (hitTime.first + hitTime.second) / 2. - (stripLength / (2 * c));

            return time;
        }

        DEFINE_ART_MODULE(SiPMHitFinder)

    } // namespace rec
} // namespace gar
