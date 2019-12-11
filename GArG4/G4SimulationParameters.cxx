//
//  G4SimulationParameters.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/3/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#include "GArG4/G4SimulationParameters.h"


namespace gar {
  namespace garg4 {
    
    static G4SimulationParameters* gInstance = nullptr;
    
    //--------------------------------------------------------------------------
    G4SimulationParameters* G4SimulationParameters::CreateInstance(fhicl::ParameterSet const& pset)
    {
      if(!gInstance) gInstance = new G4SimulationParameters(pset);
      return gInstance;
    }

    //--------------------------------------------------------------------------
    G4SimulationParameters* G4SimulationParameters::Instance()
    {
      // the instance must have been created already by CreateInstance()
      if(!gInstance)
        throw cet::exception("G4SimulationParameters")
        << "instance pointer is null, that is bad";
      
      return gInstance;
    }

    //--------------------------------------------------------------------------
    G4SimulationParameters::G4SimulationParameters(fhicl::ParameterSet const& pset)
    {
      fEnabledPhysics             = pset.get<std::vector<std::string> >("EnabledPhysics"              );
      fKeepEMShowerDaughters      = pset.get<bool                     >("KeepEMShowerDaughters", false);
      fEMShowerDaughterMatRegex   = pset.get<std::string              >("EMShowerDaughterMatRegex", ".*");
      fStoreTrajectories          = pset.get<bool                     >("StoreTrajectories",     true );
      fKineticEnergyCut           = pset.get<float                    >("KineticEnergyCut",      1.e-4);
      
      return;
    }

    //--------------------------------------------------------------------------
    G4SimulationParameters::~G4SimulationParameters()
    {
      return;
    }

    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    
    
  } // garg4
} // gar
