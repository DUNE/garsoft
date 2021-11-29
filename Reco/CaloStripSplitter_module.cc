////////////////////////////////////////////////////////////////////////
// Class:       CaloStripSplitter
// Plugin Type: producer (art v2_11_02)
// File:        CaloStripSplitter_module.cc
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/CaloHit.h"

#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include "RecoAlg/StripSplitterAlg.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include <memory>

namespace gar {
    namespace rec {

        namespace alg{
            class StripSplitterAlg;
        }

        class CaloStripSplitter : public art::EDProducer {
        public:
            explicit CaloStripSplitter(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            CaloStripSplitter(CaloStripSplitter const &) = delete;
            CaloStripSplitter(CaloStripSplitter &&) = delete;
            CaloStripSplitter & operator = (CaloStripSplitter const &) = delete;
            CaloStripSplitter & operator = (CaloStripSplitter &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        private:

            // Declare member data here.
            void CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);
            std::array<double, 3> CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> time, raw::CellID_t cID);
            float CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID);

            std::string fCaloHitLabel;  ///< label to find the right reco calo hits
            std::string fInstanceName; ///< product instance name

            bool fSaveStripEndsOnly;
            bool fSaveUnsplitHits;
            bool fSaveStripEnds;

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry

            std::unique_ptr<rec::alg::StripSplitterAlg> fSSAAlgo; //Cluster algorithm
        };

        //----------------------------------------------------------------------------
        CaloStripSplitter::CaloStripSplitter(fhicl::ParameterSet const & p) : EDProducer{p}
        {
            fCaloHitLabel = p.get<std::string>("CaloHitLabel", "calohit");
            fInstanceName =  p.get<std::string >("InstanceLabelName", "");

            fGeo     = gar::providerFrom<geo::GeometryGAr>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
            fSaveStripEndsOnly = p.get<bool>("SaveStripEndsOnly", false);
            fSaveUnsplitHits   = p.get<bool>("SaveUnsplitHits", true);

            //configure the cluster algorithm
            auto fSSAAlgoPars = p.get<fhicl::ParameterSet>("SSAAlgPars");
            fSSAAlgo = std::make_unique<rec::alg::StripSplitterAlg>(fSSAAlgoPars);
            fSaveStripEnds = fSSAAlgo->GetSaveStripEndsFlag();

            if(fSaveStripEnds) {
                MF_LOG_INFO("CaloStripSplitter_module")
                << " Saving Strip ends flag turned on! ";
            }

            art::InputTag tag(fCaloHitLabel, fInstanceName);
            consumes<std::vector<gar::rec::CaloHit>>(tag);
            produces<std::vector<gar::rec::CaloHit>>(fInstanceName);
        }

        //----------------------------------------------------------------------------
        void CaloStripSplitter::produce(art::Event & e)
        {
            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
            this->CollectHits(e, fCaloHitLabel, fInstanceName, artHits);

            //Prepare the hits for SSA
            fSSAAlgo->PrepareAlgo(artHits);
            //Perform the SS
            fSSAAlgo->DoStripSplitting();
            //Get back the vector of split hits
            std::vector<const gar::rec::CaloHit*> splitHits = fSSAAlgo->getSplitHits();
            //Get back the vector of unsplit hits
            std::vector<const gar::rec::CaloHit*> unsplitHits = fSSAAlgo->getUnSplitHits();

            // make an art::PtrVector of the hits
            std::unique_ptr< std::vector<gar::rec::CaloHit> > HitCol(new std::vector<gar::rec::CaloHit>);

            if(not fSaveStripEndsOnly) {
                //Copy the split hits to the collection
                for(auto const &it : splitHits) {

                    //Need to correct the time for the strips
                    float energy = it->Energy();
                    std::pair<float, float> time = it->Time();
                    raw::CellID_t cellID = it->CellID();
                    const float *pos = it->Position();
                    float newpos[3] = { pos[0], pos[1], pos[2] };
                    float newtime = time.first;
                    unsigned int layer = it->Layer();
                    const std::array<double, 3> point = { pos[0], pos[1], pos[2] };

                    if(not fGeo->isTile(point, cellID)) {
                        newtime = this->CorrectStripHitTime(newpos[0], newpos[1], newpos[2], time, cellID);
                    }

                    rec::CaloHit hit(energy, newtime, newpos, cellID, layer);
                    HitCol->emplace_back(hit);
                }

                if( fSaveUnsplitHits ) {
                    //Copy the unsplit hits to the collection
                    for(auto const &it : unsplitHits) {

                        //Need to correct the position of the un-used strips
                        float energy = it->Energy();
                        std::pair<float, float> time = it->Time();
                        raw::CellID_t cellID = it->CellID();
                        const float *pos = it->Position();
                        float newpos[3] = { pos[0], pos[1], pos[2] };
                        float newtime = time.first;
                        unsigned int layer = it->Layer();
                        const std::array<double, 3> point = { pos[0], pos[1], pos[2] };

                        if(not fGeo->isTile(point, cellID)) {
                            // Need to correct for the position of these based on time information
                            std::array<double, 3> strip_pos = this->CalculateStripHitPosition(newpos[0], newpos[1], newpos[2], time, cellID);
                            newpos[0] = strip_pos[0];
                            newpos[1] = strip_pos[1];
                            newpos[2] = strip_pos[2];
                            newtime = this->CorrectStripHitTime(newpos[0], newpos[1], newpos[2], time, cellID);
                        }

                        rec::CaloHit hit(energy, newtime, newpos, cellID, layer);

                        HitCol->emplace_back(hit);
                    }
                }
            }

            //Copy the strip end hits to the collection
            if(fSaveStripEnds) {
                std::vector<const gar::rec::CaloHit*> stripEndsHits = fSSAAlgo->getStripEndsHits();
                for(auto const &it : stripEndsHits)
                HitCol->emplace_back(*it);

                std::vector<const gar::rec::CaloHit*> stripInterHits = fSSAAlgo->getIntersectionHits();
                for(auto const &it : stripInterHits)
                HitCol->emplace_back(*it);
            }

            e.put(std::move(HitCol), fInstanceName);

            return;
        }

        //----------------------------------------------------------------------------
        void CaloStripSplitter::CollectHits(const art::Event &evt, const std::string &label, const std::string &instance, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
        {
	    art::InputTag itag(label,instance);
	    auto theHits = evt.getHandle< std::vector<gar::rec::CaloHit> >(itag);
            if (!theHits) return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        //----------------------------------------------------------------------------
        std::array<double, 3> CaloStripSplitter::CalculateStripHitPosition(float x, float y, float z, std::pair<float, float> time, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            std::array<double, 3> pointLocal;
            gar::geo::LocalTransformation<TGeoHMatrix> trans;
            fGeo->WorldToLocal(point, pointLocal, trans);

            //Calculate the position of the hit based on the corrected time along the strip
            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float xlocal = c * ( time.first - time.second ) / 2.;

            std::array<double, 3> local_back = fGeo->ReconstructStripHitPosition(point, pointLocal, xlocal, cID);
            std::array<double, 3> world_back;
            fGeo->LocalToWorld(local_back, world_back, trans);

            return world_back;
        }

        //----------------------------------------------------------------------------
        float CaloStripSplitter::CorrectStripHitTime(float x, float y, float z, std::pair<float, float> hitTime, raw::CellID_t cID)
        {
            std::array<double, 3> point = {x, y, z};
            double stripLength = fGeo->getStripLength(point, cID); // in cm

            float c = (CLHEP::c_light * CLHEP::mm / CLHEP::ns) / CLHEP::cm; // in cm/ns
            float time = (hitTime.first + hitTime.second) / 2. - (stripLength / (2 * c));

            return time;
        }

        DEFINE_ART_MODULE(CaloStripSplitter)

    } // namespace rec
} // namespace gar
