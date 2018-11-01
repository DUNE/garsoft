//
//  ElectronDriftAlg.h
//
//  Algorithm base class for drifting electrons in an electric field
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef GAR_READOUTSIMULATION_ElectronDriftAlg_h
#define GAR_READOUTSIMULATION_ElectronDriftAlg_h

#include <vector>

#include "CLHEP/Random/RandGauss.h"


namespace fhicl {
  class ParameterSet;
}

namespace gar {
  
  namespace sdp{
    class EnergyDeposit;
  }
  
  namespace rosim{
    
    class ElectronDriftInfo{
      
    public:
      
      ElectronDriftInfo();
      
      void Reset(std::vector<double> & xPos,
                 std::vector<double> & yPos,
                 std::vector<double> & zPos,
                 std::vector<double> & time,
                 std::vector<int   > & size);
      
      std::vector<double> const& ClusterXPos() const { return fClusterXPos; }
      std::vector<double> const& ClusterYPos() const { return fClusterYPos; }
      std::vector<double> const& ClusterZPos() const { return fClusterZPos; }
      std::vector<double> const& ClusterTime() const { return fClusterTime; }
      std::vector<int   > const& ClusterSize() const { return fClusterSize; }
      
    private:
      
      std::vector<double> fClusterXPos;  ///< x positions of each cluster drifted
      std::vector<double> fClusterYPos;  ///< y positions of each cluster drifted
      std::vector<double> fClusterZPos;  ///< z positions of each cluster drifted
      std::vector<double> fClusterTime;  ///< arrival time of each cluster drifted
      std::vector<int   > fClusterSize;  ///< size of each cluster drifted
      
    };
    
    class ElectronDriftAlg{
      
    public:
      
      ElectronDriftAlg(CLHEP::HepRandomEngine      & engine,
                       fhicl::ParameterSet    const& pset);
      virtual ~ElectronDriftAlg();
      
      // The user must call IonizationAndScintillation::Instance()->Reset(dep)
      // from within the implementation of this method in the derived class
      virtual void DriftElectronsToReadout(gar::sdp::EnergyDeposit       const& dep,
                                           gar::rosim::ElectronDriftInfo      & driftInfo) = 0;
      
    protected:
      
      double                  fDriftVelocity;         ///< electron drift velocity
      double                  fInverseVelocity;       ///< stored for computational convenience
      double                  fLifetimeCorrection;    ///< electron lifetime correction in negative ms
      double                  fLongitudinalDiffusion; ///< diffusion along the drift in microns/sqrt(cm)
      double                  fLongDiffConst;         ///< stored for computational convenience in sqrt(cm)
      double                  fTransverseDiffusion;   ///< diffusion transverse to the drift in microns/sqrt(cm)
      double                  fTransDiffConst;        ///< stored for computational convenience in sqrt(cm)
      double                  fFanoFactor;            ///< Fano factor
      CLHEP::HepRandomEngine& fEngine;                ///< random number engine
    };
    
  }
} // gar

#endif /* GAR_READOUTSIMULATION_ElectronDriftAlg_h */
