////////////////////////////////////////////////////////////////////////
/// \file  TrackCheater_module.cc
/// \brief Create rec::Track products using MC truth information
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GAR_MCCHEATER_TRACKCHEATER
#define GAR_MCCHEATER_TRACKCHEATER

// C++ Includes
#include <memory>
#include <vector> // std::ostringstream
#include <iostream>

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "cetlib_except/exception.h"

#include "nusimdata/SimulationBase/MCParticle.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"

// ROOT includes

#include "TMath.h"
#include "TDatabasePDG.h"

// Forward declarations

namespace gar {
  namespace cheat {
    
    class TrackCheater : public ::art::EDProducer{
    public:
      
      // Standard constructor and destructor for an FMWK module.
      explicit TrackCheater(fhicl::ParameterSet const& pset);
      virtual ~TrackCheater();
      
      void produce (::art::Event& evt);
      void reconfigure(fhicl::ParameterSet const& pset);
      
    private:
      
      std::string fHitLabel; ///< label of module creating hits
      std::string fG4Label;  ///< label of module creating mc particles
      
    };
    
  } // cheat
  
  namespace cheat {
    
    //--------------------------------------------------------------------------
    TrackCheater::TrackCheater(fhicl::ParameterSet const& pset)
    {
      this->reconfigure(pset);
      
      produces< std::vector<rec::Track>                    >();
      produces< ::art::Assns<rec::Hit,         rec::Track> >();
      produces< ::art::Assns<simb::MCParticle, rec::Track> >();
      
      return;
    }
    
    //--------------------------------------------------------------------------
    TrackCheater::~TrackCheater()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    void TrackCheater::reconfigure(fhicl::ParameterSet const& pset)
    {
      fHitLabel = pset.get<std::string >("HitLabel",      "hit"  );
      fG4Label  = pset.get<std::string >("G4ModuleLabel", "geant");
      return;
    }
    
    //--------------------------------------------------------------------------
    void TrackCheater::produce(::art::Event & evt)
    {
      if( evt.isRealData() ){
        throw cet::exception("TrackCheater")
	  << "Attempting to cheat the track reconstruction for real data, "
	  << "that will never work";
        return;
      }
      
      std::unique_ptr< std::vector<rec::Track>                    > trkCol     (new std::vector<rec::Track>                   );
      std::unique_ptr< ::art::Assns<rec::Hit,         rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,         rec::Track>);
      std::unique_ptr< ::art::Assns<simb::MCParticle, rec::Track> > pTrkAssns  (new ::art::Assns<simb::MCParticle, rec::Track>);
      
      // get the MCParticles from the event and make a vector of Ptrs of them
      auto partCol = evt.getValidHandle< std::vector<simb::MCParticle> >(fG4Label);
      std::vector<art::Ptr<simb::MCParticle> > partVec;
      art::fill_ptr_vector(partVec, partCol);

      // get the FindMany for the digit to EnergyDeposit association
      ::art::FindManyP<rec::Hit> fmphit(partVec, evt, fHitLabel);
      
      // test that the FindMany is valid - don't worry about the digCol, the
      // getValidHandle throws if it can't find a valid handle
      if(!fmphit.isValid() ){
        throw cet::exception("TrackCheater")
	  << "Unable to find valid FindMany<Hit> "
	  << fmphit.isValid()
	  << " this is a problem for cheating";
      }
      
      std::vector< art::Ptr<rec::Hit> > hits;
      simb::MCParticle const* part      = nullptr;
      float                   length    = 0.;
      float                   momentum  = 0.;
      float                   end_momentum  = 0.;
      float                   vtx[3]    = {0.};
      float                   end[3]    = {0.};
      float                   vtxDir[3] = {0.};
      float                   endDir[3] = {0.};
      
      for(size_t p = 0; p < partVec.size(); ++p){
        fmphit.get(p, hits);
        if(hits.size() < 1) continue;
        
        part = partVec[p].get();
        
	// skip this -- we want to track electrons in the TPC and pizeros don't make tracks
        // ignore if we have an EM shower
        //if(std::abs(part->PdgCode()) == 11 ||
        //   std::abs(part->PdgCode()) == 111) continue;
        
        // find the length of the track by getting the distance between each hit
        length = 0.;
        for(size_t h = 1; h < hits.size(); ++h)
          length += std::sqrt(TMath::Sq(hits[h]->Position()[0] - hits[h-1]->Position()[0]) +
                              TMath::Sq(hits[h]->Position()[1] - hits[h-1]->Position()[1]) +
                              TMath::Sq(hits[h]->Position()[2] - hits[h-1]->Position()[2]));
        
        momentum  = part->Momentum().Mag();
        end_momentum  = part->EndMomentum().Mag();
        vtx[0]    = part->Vx();
        vtx[1]    = part->Vx();
        vtx[2]    = part->Vx();
        vtxDir[0] = part->Px() / momentum;
        vtxDir[1] = part->Py() / momentum;
        vtxDir[2] = part->Pz() / momentum;
        end[0]    = part->EndX();
        end[1]    = part->EndY();
        end[2]    = part->EndZ();
        endDir[0] = part->EndPx() / end_momentum;
        endDir[1] = part->EndPy() / end_momentum;
        endDir[2] = part->EndPz() / end_momentum;
	int pdgid = part->PdgCode();
	TDatabasePDG *PDGdb = TDatabasePDG::Instance(); 
	int charge = PDGdb->GetParticle(pdgid)->Charge();   // nb we set the charge to +- 1 in the tracker, but get the full charge here
	trkCol->emplace_back(length,
                             momentum,
			     end_momentum,
                             vtx,
                             end,
                             vtxDir,
                             endDir,
			     hits.size(),
			     charge,
			     0);  // zero timestamp for now

        // make the hit to MCParticle association
        util::CreateAssn(*this, evt, *trkCol, partVec[p], *pTrkAssns);

        // make the hit to track associations
        for(auto hit : hits) util::CreateAssn(*this, evt, *trkCol, hit, *hitTrkAssns);

      } // end loop over particles
      
      evt.put(std::move(trkCol     ));
      evt.put(std::move(hitTrkAssns));
      evt.put(std::move(pTrkAssns  ));
      
      return;
    }
    
    
  } // cheat
  
  namespace cheat {
    
    DEFINE_ART_MODULE(TrackCheater)
    
  } // cheat
} // gar
#endif // GAR_MCCHEATER_TRACKCHEATER
