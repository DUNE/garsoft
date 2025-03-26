////////////////////////////////////////////////////////////////////////
/// \file  GArG4Ana.h
/// \brief Check of Geant4 to run the GArSoft detector simulation
///
/// \version $Id: GArG4.h,v 1.11 2010/06/04 21:47:27 bjpjones Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GARG4MNSP_H
#define GARG4MNSP_H

#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ParticleWithCuts.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4ProcessVector.hh"
#include "Geant4/G4WrapperProcess.hh"
#include "Geant4/globals.hh"

#include "Geant4/G4Material.hh"
#include "Geant4/G4MaterialTable.hh"
#include "Geant4/G4MuonNuclearProcess.hh"
#include "Geant4/G4ios.hh"

//#include "Geant4/G4DataQuestionaire.hh" //depreciated

namespace gar {
  namespace garg4 {

    class MuNuclearSplittingProcess : public G4WrapperProcess {
      // Override PostStepDoIt method
    public:
      MuNuclearSplittingProcess(){};
      ~MuNuclearSplittingProcess(){};

      void SetNSplit(G4int nTrx) { fNSplit = nTrx; };
      void SetIsActive(G4bool doIt) { fActive = doIt; };

    private:
      // Data members
      G4int fNSplit;
      G4bool fActive;
      G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step& step);
    };

  } // end namespace
} // gar

#endif // MNSP
