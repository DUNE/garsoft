////////////////////////////////////////////////////////////////////////
//
// GENIEEventFilter_module class:
// Algoritm to produce a filtered event file having events with user-defined 
// GENIE stauscode==1 (not GEANT!) particles in MCTruth
//
// eldwan.brianne@desy.de
//
////////////////////////////////////////////////////////////////////////

//Framework Includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"

//nusim Includes
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace gar{
    namespace filt {

        class FSPEventFilter : public art::EDFilter  {

        public:

            explicit FSPEventFilter(fhicl::ParameterSet const& );

            // Filters should not be copied or assigned.
            FSPEventFilter(FSPEventFilter const &) = delete;
            FSPEventFilter(FSPEventFilter &&) = delete;
            FSPEventFilter & operator = (FSPEventFilter const &) = delete;
            FSPEventFilter & operator = (FSPEventFilter &&) = delete;

            bool filter(art::Event& evt);
            void beginJob();

        private:

            std::string fGeneratorLabel;
            std::vector<int> fPDG;

            bool isMatched(std::vector<int> const& a, std::vector<int> const& b) const;

        }; // class GENIEEventFilter

        //-------------------------------------------------
        FSPEventFilter::FSPEventFilter(fhicl::ParameterSet const & pset)
        : EDFilter{pset}
        {
            fGeneratorLabel = pset.get< std::string      >("GeneratorModuleLabel", "generator");
            fPDG            = pset.get< std::vector<int> >("PDG");
        }

        //-------------------------------------------------
        void FSPEventFilter::beginJob()
        {

        }

        //-------------------------------------------------
        bool FSPEventFilter::filter(art::Event &evt)
        {
	    auto mcthandlelist = evt.getHandle< std::vector<simb::MCTruth> >(fGeneratorLabel);
            if (!mcthandlelist) {
                throw cet::exception("FSPEventFilter") << " No simb::MCTruth branch."
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }

            art::Ptr<simb::MCTruth> mcp(mcthandlelist, 0);
            std::vector<int> FSP;

            for(int i = 0; i < mcp->NParticles(); ++i){
                simb::MCParticle part(mcp->GetParticle(i));

                if (part.StatusCode()== 1)
                    FSP.push_back(part.PdgCode());
            }

            return isMatched(fPDG, FSP); // returns true if the user-defined fPDG exist(s) in the final state particles
        }

        //------------------------------------------------
        bool FSPEventFilter::isMatched(std::vector<int> const& a, std::vector<int> const& b) const
        {
            for (auto const a_int : a) {
                for (auto const b_int : b) {
                    if (a_int == b_int) {
                        return true;
                    }
                }
            }
            return false;
        }

    } //namespace filt
} //namespace gar

//--------------------------------------------------
namespace gar {
    namespace filt {
        DEFINE_ART_MODULE(FSPEventFilter)
    } //namespace gar
} //namespace filt
