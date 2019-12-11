//
//  G4SimulationParameters.h
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/3/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef G4SimulationParameters_hpp
#define G4SimulationParameters_hpp

#include <vector>
#include <iostream>

#include "fhiclcpp/ParameterSet.h"


namespace gar {
  namespace garg4 {
    
    class G4SimulationParameters
    {
    public:
      
      static G4SimulationParameters* CreateInstance(fhicl::ParameterSet const& pset);
      static G4SimulationParameters* Instance();
      
      std::vector<std::string> const& EnabledPhysics()        const { return fEnabledPhysics;        }
      double                          KineticEnergyCut()      const { return fKineticEnergyCut;      }
      bool                            StoreTrajectories()     const { return fStoreTrajectories;     }
      bool                            KeepEMShowerDaughters() const { return fKeepEMShowerDaughters; }
      std::string                     EMShowerDaughterMatRegex() const { return fEMShowerDaughterMatRegex; }
      
    private:
      
      explicit G4SimulationParameters(fhicl::ParameterSet const& pset);
      ~G4SimulationParameters();
      
      std::vector<std::string> fEnabledPhysics;        ///< list of enabled physics processes
      float                    fKineticEnergyCut;      ///< don't track particles below this energy
      bool                     fStoreTrajectories;     ///< keep particle trajectory information or not
      bool                     fKeepEMShowerDaughters; ///< keep EM shower daughters or not
      std::string              fEMShowerDaughterMatRegex; ///< keep EM shower daughters only in these materials
    };
    
  } // garg4
} // gar


#endif /* G4SimulationParameters_hpp */
