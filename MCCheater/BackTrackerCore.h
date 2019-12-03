//
//  BackTrackerCore.h
//
//  Created by Brian Rebel on 3/13/17.
//

#ifndef GAR_CHEAT_BackTrackerCore_h
#define GAR_CHEAT_BackTrackerCore_h

#include "canvas/Persistency/Common/Ptr.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/ParticleNavigation/ParticleList.h"
#include "nutools/ParticleNavigation/EveIdCalculator.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "DetectorInfo/DetectorClocksService.h"
#include "SimulationDataProducts/EnergyDeposit.h"
#include "ReconstructionDataProducts/Hit.h"
#include "RawDataProducts/RawDigit.h"
#include "Geometry/Geometry.h"

///code to link reconstructed objects back to the MC truth information

namespace sim {
  class ParticleList;
  class EveIdCalculator;
}

namespace gar{
  namespace cheat{

      // make convenience structs for doing the backtracking
    struct HitIDE {
      
      int   trackID;      ///< Geant4 supplied trackID
      float energyFrac;   ///< fraction of hit energy from the particle with this trackID
      float numElectrons; ///< number of ionization electrons for this trackID
      HitIDE()
      {}
      
      HitIDE(int   id,
             float ef,
             float e)
      : trackID(id)
      , energyFrac(ef)
      , numElectrons(e)
      {}
      
    };
    
    class BackTrackerCore {
    
      public:
      
      explicit BackTrackerCore(fhicl::ParameterSet const& pset);
      ~BackTrackerCore();
      
      void reconfigure(fhicl::ParameterSet const& pset);
    
      // Get a reference to the ParticleList
      sim::ParticleList const& ParticleList()  const { return fParticleList; }
      
      // Set the EveIdCalculator for the owned ParticleList
      void AdoptEveIdCalculator(sim::EveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }
      
      // Return a pointer to the simb::MCParticle object corresponding to
      // the given TrackID
      const simb::MCParticle* TrackIDToParticle(int const& id)       const;
      const simb::MCParticle* TrackIDToMotherParticle(int const& id) const;
      
      // Get ::art::Ptr<> to simb::MCTruth and related information
      ::art::Ptr<simb::MCTruth>                const& TrackIDToMCTruth(int const& id)                          const;
      ::art::Ptr<simb::MCTruth>                const& ParticleToMCTruth(const simb::MCParticle* p)             const;
      std::vector<const simb::MCParticle*>            MCTruthToParticles(::art::Ptr<simb::MCTruth> const& mct) const;
      std::vector< ::art::Ptr<simb::MCTruth> > const& MCTruthVector() const { return fMCTruthList; };
      
      // this method determines the appropriate EnergyDeposits for a given TDC
      // on a channel
      std::map<size_t, const sdp::EnergyDeposit*> ChannelTDCToEnergyDeposit(raw::Channel_t channel) const;
      
      // this method will return the Geant4 track IDs of
      // the particles contributing ionization electrons to the identified hit
      std::vector<HitIDE> HitToTrackID(::art::Ptr<gar::rec::Hit> const& hit) const;
      
      // method to return a subset of allhits that are matched to a list of TrackIDs
      std::vector<std::vector<::art::Ptr<gar::rec::Hit>>> const TrackIDsToHits(std::vector<::art::Ptr<gar::rec::Hit>> const& allhits,
                                                                               std::vector<int>                       const& tkIDs) const;
      
      // method to return the EveIDs of particles contributing ionization
      // electrons to the identified hit
      std::vector<HitIDE> HitToEveID(::art::Ptr<gar::rec::Hit> const& hit) const;
      
      // method to return the XYZ position of the weighted average energy deposition for a given hit
      std::vector<float>  HitToXYZ(::art::Ptr<gar::rec::Hit> const& hit) const;
      
      // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
      double HitCollectionPurity(std::set<int>                                   trackIDs,
                                 std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                 bool                                            weightByCharge=false) const;
      
      // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
      // represented in a collection of hits
      double HitCollectionEfficiency(std::set<int>                                   trackIDs,
                                     std::vector< ::art::Ptr<gar::rec::Hit> > const& hits,
                                     std::vector< ::art::Ptr<gar::rec::Hit> > const& allhits,
                                     bool                                            weightByCharge=false) const;
      
      // method to return all EveIDs corresponding to the current sim::ParticleList
      std::set<int> GetSetOfEveIDs() const;
      
      // method to return all TrackIDs corresponding to the current sim::ParticleList
      std::set<int> GetSetOfTrackIDs() const;
      
      // method to return all EveIDs corresponding to the given list of hits
      std::set<int> GetSetOfEveIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits) const;
      
      // method to return all TrackIDs corresponding to the given list of hits
      std::set<int> GetSetOfTrackIDs(std::vector< ::art::Ptr<gar::rec::Hit> > const& hits) const;
      
    private:
      
      void ChannelToTrackID(std::vector<HitIDE>      & hitIDEs,
                            gar::raw::Channel_t const& channel,
                            double              const  start,
                            double              const  end) const;
      
    protected:
      
      const detinfo::DetectorClocks*                        fClocks;                ///< Detector clock information
      sim::ParticleList                                     fParticleList;          ///< ParticleList to map track ID to sim::Particle
      std::vector< std::vector<const sdp::EnergyDeposit*> > fChannelToEDepCol;      ///< convenience collections of EnergyDeposits for each channel
      std::vector< ::art::Ptr<simb::MCTruth> >              fMCTruthList;           ///< all the MCTruths for the event
      std::map<int, int>                                    fTrackIDToMCTruthIndex; ///< map of track ids to MCTruthList entry
      std::string                                           fG4ModuleLabel;         ///< label for geant4 module
      std::string                                           fRawDataLabel;          ///< label for geant4 module
      double                                                fMinHitEnergyFraction;  ///< minimum fraction of energy a track id has to
                                                                                    ///< contribute to a hit to be counted in
                                                                                    ///< purity and efficiency calculations
                                                                                    ///< based on hit collections
      double                                                fInverseVelocity;       ///< inverse drift velocity
      double                                                fLongDiffConst;         ///< longitudinal diffusion constant

      const geo::GeometryCore*                              fGeo;                   ///< pointer to the geometry
      bool                                                  fDisableRebuild;        ///< for switching off backtracker's rebuild of the MCParticle tables
    };
  } // namespace
}


#endif /* GAR_CHEAT_BackTrackerCore_h */
