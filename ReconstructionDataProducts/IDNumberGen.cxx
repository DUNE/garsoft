//  IDNumberGen.cxx
//  garsoft-mrb
//
//  Created by Leo Bellantoni on 2 Aug, 2019
//  With lots of help from Paul Russo :)



#include "ReconstructionDataProducts/IDNumberGen.h"



namespace gar {
    namespace rec {



        IDNumberGen* IDNumberGen::create(IDNumber iniValue) { 

            static IDNumberGen* tmp = new IDNumberGen();
            // If creating an instance with an initial value
			if (iniValue != std::numeric_limits<IDNumber>::max()) {
                // and if the iniValue isn't one of the previously used ones
                std::set<IDNumber>::iterator itr = previousInitializers.find(iniValue);
                if (itr==previousInitializers.end()) {
                    previousInitializers.insert(iniValue);
                    nextOneToMake.store(iniValue);
                }
            }
            // Even if nextOneToMake.store() not called, still fine that tmp is
            // uninitialized as getNewOne() only pulls info from a static field.
            return tmp;
        }



        IDNumberGen::IDNumber IDNumberGen::getNewOne() {
            IDNumberGen::IDNumber retval = nextOneToMake;
            nextOneToMake.store(nextOneToMake+1);
            return retval;
        }



        void IDNumberGen::newEventReset() {
            previousInitializers.clear();
            nextOneToMake.store(std::numeric_limits<IDNumber>::max());
            return;
        }



        std::atomic<IDNumberGen::IDNumber> IDNumberGen::nextOneToMake = std::numeric_limits<IDNumber>::max();
        std::set<IDNumberGen::IDNumber>    IDNumberGen::previousInitializers = {};

        IDNumberGen::IDNumberGen() {}

    } // rec
} // gar
