////////////////////////////////////////////////////////////////////////
/// \file  PhysicsList.h
/// \brief Create the physics lists to be used by Geant4.
///
/// \version $Id: PhysicsList.cxx,v 1.3 2010/07/26 15:16:35 bjpjones Exp $
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
//
// Don't be too confused by the names.  PhysicsList.h and
// PhysicsList.cxx define what the name "garg4::PhysicList" means.
// However, that definition is mainly in terms of
// garg4::ModularPhysicsList, a class that inherits from
// G4VModularPhysicsList.

#include "GArG4/PhysicsList.h"

#include "Geant4/G4ChargeExchange.hh"
#include "Geant4/G4ChargeExchangeProcess.hh"
#include "Geant4/G4ParallelWorldScoringProcess.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4VModularPhysicsList.hh"
#include "Geant4/G4VUserPhysicsList.hh"
//#include "nug4/G4Base/G4PhysListFactorySingleton.hh"
///PHYSLISTREG3(garg4,PhysicsList,garg4::PhysicsList)
#include "Geant4/G4PhysListStamper.hh"
G4_DECLARE_PHYSLIST_FACTORY_NS(gar::garg4::PhysicsList, garg4, PhysicsList);

namespace gar {

#define G4MT_physicsVector ((G4VMPLsubInstanceManager.offset[g4vmplInstanceID]).physicsVector)

  namespace garg4 {

    // Constructor: call the G4 constructor.
    ModularPhysicsList::ModularPhysicsList() : G4VModularPhysicsList() {}

    // Destructor; C++ will automatically call the G4VModulePhysicsList
    // destructor, so we need do nothing here.
    ModularPhysicsList::~ModularPhysicsList() {}

    // This is the method we have to modify to use the Geant4 parallel geometries.
    void ModularPhysicsList::ConstructProcess()
    {
      // We don't need to modify G4VModularPhysicsList's
      // AddTransportation method.  Just invoke it directly.
      G4VModularPhysicsList::AddTransportation();

      // This code is also unchanged from
      // G4VModularPhysicsList::ConstructProcess(); it means "activate
      // the physics processes and particle combinations we've specified
      // in physicsVector."  The physicsVector is built up by the
      // pre-supplied Geant4 physics list; see PhysicsList.h for the
      // name of that list.  "physicsVector" is defined in
      // G4VModularPhysicsList.hh.
      G4PhysConstVector::iterator itr;
      for (itr = G4MT_physicsVector->begin(); itr != G4MT_physicsVector->end(); ++itr) {
        (*itr)->ConstructProcess();
      }
    } // end ConstructProcess

  } // namespace garg4
} // namespace gar
