////////////////////////////////////////////////////////////////////////
// Class:       EventInit
// Plugin Type: producer (art v3_00_0)
// File:        EventInit_module.cc
//
// A little module for anything that has to be done at the start of
// processing every event.
// Leo Bellantoni, 8 Aug 2019
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

//#include "ReconstructionDataProducts/Track.h"
//#include "ReconstructionDataProducts/Vertex.h"
//#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/IDNumberGen.h"



namespace gar {
    namespace rec {



        class EventInit : public art::EDProducer {
        public:
            explicit EventInit(fhicl::ParameterSet const & p);
            // The compiler-generated destructor is fine for non-base
            // classes without bare pointers or other resource use.

            // Plugins should not be copied or assigned.
            EventInit(EventInit const &) = delete;
            EventInit(EventInit &&) = delete;
            EventInit & operator = (EventInit const &) = delete;
            EventInit & operator = (EventInit &&) = delete;

            // Required functions.
            void produce(art::Event & e) override;

        private:
            int fVerbosity;
        };



        EventInit::EventInit(fhicl::ParameterSet const & p) {
            fVerbosity = p.get<int>("Verbosity", 0);
        }



        void EventInit::produce(art::Event & e) {
            IDNumberGen::create()->newEventReset();
            return;
        }



        DEFINE_ART_MODULE(EventInit)

    } // namespace rec
} // namespace gar
