////////////////////////////////////////////////////////////////////////
//
// GENIEEventFilter_module class:
// Algorithm to produce a filtered event file having events with user-defined 
// GENIE stauscode==1 (not GEANT!) particles in MCTruth.  MCTruth 
// particles are from the primary vertex.
// For some reason, this module does not produce a usable output when
// put into a GENIEGen job.
//
// eldwan.brianne@desy.de
//
// Expanded for location-of-interaction and multiple GENIE & GEANT
// generators
// Leo Bellantoni 1 Dec 2023
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
//#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "CoreUtils/ServiceUtil.h"
#include "Geometry/GeometryGAr.h"



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
            std::vector<std::string> fGeneratorLabels;
            std::vector<std::string> fGENIEGeneratorLabels;

            std::vector<int>         fPDG;
            std::vector<std::string> fVolumes;
            double                   fPminCut;

            const geo::GeometryCore* fGeo; ///< pointer to the geometry

            bool isMatched(std::vector<int> const& a, std::vector<int> const& b) const;

        }; // class GENIEEventFilter



        //-------------------------------------------------
        FSPEventFilter::FSPEventFilter(fhicl::ParameterSet const & pset)
        : EDFilter{pset} {
            fGeo     = gar::providerFrom<geo::GeometryGAr>();

            bool usegenlabels =
                pset.get_if_present<std::vector<std::string> >("GeneratorLabels",fGeneratorLabels);
            if (!usegenlabels) fGeneratorLabels.clear();

            bool usegeniegenlabels  =
                pset.get_if_present<std::vector<std::string> >("GENIEGeneratorLabels",fGENIEGeneratorLabels);
            if (!usegeniegenlabels) fGENIEGeneratorLabels.clear();

            fPDG     = pset.get< std::vector<int> >("PDG");
            fVolumes = pset.get< std::vector<std::string> >("volumes");
            fPminCut = pset.get<double>("PminCut");
            return;
        }



        //-------------------------------------------------
        void FSPEventFilter::beginJob() {
        }



        //-------------------------------------------------
        bool FSPEventFilter::filter(art::Event &evt) {

            std::vector< art::Handle< std::vector<simb::MCTruth> > > mcthandlelist;

            if (fGeneratorLabels.size()<1) {
                mcthandlelist = evt.getMany<std::vector<simb::MCTruth> >(); // get them all (even if there are none)
            } else {
                mcthandlelist.resize(fGeneratorLabels.size());
                for (size_t i=0; i< fGeneratorLabels.size(); ++i) {
                    // complain if we wanted a specific one but didn't find it
                    mcthandlelist.at(i) = evt.getHandle<std::vector<simb::MCTruth> >(fGeneratorLabels.at(i));
                    if (!mcthandlelist.at(i)) {
                        throw cet::exception("anatree") << " No simb::MCTruth branch."
                            << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
                    }
                }
            }


            // examine MCTruth info
            for (size_t imchl = 0; imchl < mcthandlelist.size(); ++imchl) {
                for ( simb::MCTruth const& mct : (*mcthandlelist.at(imchl)) ) {
                    bool inVolume = fVolumes.size()==0;
                    if (!inVolume) {
                        // Look for neutrino interactions in the defined place
						if (mct.NeutrinoSet()) {
                            simb::MCNeutrino nu = mct.GetNeutrino();
                            TVector3 place(nu.Nu().EndX(),nu.Nu().EndY(),nu.Nu().EndZ());
                            for (auto vol : fVolumes) {
                               if (vol=="GArTPC" &&  fGeo->PointInGArTPC(place))       inVolume = true;
                               if (vol=="ECAL"   && (fGeo->PointInECALBarrel(place) ||
                                                     fGeo->PointInECALEndcap(place)) ) inVolume = true;
                               if (vol=="MuID"   && (fGeo->PointInMuIDBarrel(place) ||
                                                     fGeo->PointInMuIDEndcap(place)) ) inVolume = true;
                            }
						}
                    }
                    if (!inVolume) continue;

                    // The interaction is in the right place... does it have the right FSP, with right momentum?
                    std::vector<int> FSP;

                    for (int i = 0; i < mct.NParticles(); ++i){
                        simb::MCParticle part(mct.GetParticle(i));
                        if (part.StatusCode()== 1 &&
                            part.P()>fPminCut) FSP.push_back(part.PdgCode());
                    }

                    if (isMatched(fPDG,FSP)) {
                        return true;
                    }
                } // end of loop over MCTruths
            } // end of loop over generators

            bool testval = false;
            return testval;
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
