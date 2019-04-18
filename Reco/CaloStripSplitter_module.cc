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

#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
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
            void CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector);

            std::string fCaloHitLabel;  ///< label to find the right reco calo hits
            int fVerbosity;

            const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
            const geo::GeometryCore*            fGeo;          ///< pointer to the geometry

            std::unique_ptr<rec::alg::StripSplitterAlg> fSSAAlgo; //Cluster algorithm
        };


        CaloStripSplitter::CaloStripSplitter(fhicl::ParameterSet const & p)
        {
            fCaloHitLabel = p.get<std::string>("CaloHitLabel", "calohit");
            fVerbosity = p.get<int>("Verbosity", 0);

            fGeo     = gar::providerFrom<geo::Geometry>();
            fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            //configure the cluster algorithm
            auto fSSAAlgoPars = p.get<fhicl::ParameterSet>("SSAAlgPars");
            fSSAAlgo = std::make_unique<rec::alg::StripSplitterAlg>(fSSAAlgoPars);

            consumes< std::vector<gar::rec::CaloHit> >(fCaloHitLabel);
            produces< std::vector<gar::rec::CaloHit> >();
        }

        void CaloStripSplitter::produce(art::Event & e)
        {
            //Collect the hits to be passed to the algo
            std::vector< art::Ptr<gar::rec::CaloHit> > artHits;
            this->CollectHits(e, fCaloHitLabel, artHits);

            //Prepare the hits for SSA
            fSSAAlgo->PrepareAlgo(artHits);
            //Perform the SS
            fSSAAlgo->DoStripSplitting();
            //Get back the vector of split hits
            std::vector<gar::rec::CaloHit*> splitHits = fSSAAlgo->getSplitHits();
            //Get back the vector of unsplit hits
            std::vector<gar::rec::CaloHit*> unsplitHits = fSSAAlgo->getUnSplitHits();

            // make an art::PtrVector of the hits
            std::unique_ptr< std::vector<gar::rec::CaloHit> > HitCol(new std::vector<gar::rec::CaloHit>);

            //Copy the split hits to the collection
            for(auto it : splitHits)
            {
                gar::rec::CaloHit hit(*it);
                HitCol->push_back(hit);
            }

            //Copy the unsplit hits to the collection
            for(auto it : unsplitHits)
            {
                gar::rec::CaloHit hit(*it);
                HitCol->push_back(hit);
            }

            std::cout << "Before " << artHits.size() << std::endl;
            std::cout << "After " << HitCol->size() << std::endl;

            e.put(std::move(HitCol));

            return;
        }

        void CaloStripSplitter::CollectHits(const art::Event &evt, const std::string &label, std::vector< art::Ptr<gar::rec::CaloHit> > &hitVector)
        {
            art::Handle< std::vector<gar::rec::CaloHit> > theHits;
            evt.getByLabel(label, theHits);

            if (!theHits.isValid())
            return;

            for (unsigned int i = 0; i < theHits->size(); ++i)
            {
                const art::Ptr<gar::rec::CaloHit> hit(theHits, i);
                hitVector.push_back(hit);
            }
        }

        DEFINE_ART_MODULE(CaloStripSplitter)

    } // namespace rec
} // namespace gar
