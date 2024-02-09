//
//  ElectronDriftStandardAlg.h
//
//  Default implementation of an electron drift algorithm
//
//  Created by Brian Rebel on 11/18/16.
//  Copyright © 2016 Brian Rebel. All rights reserved.
//

#ifndef ElectronDriftStandardAlgTOAD_h
#define ElectronDriftStandardAlgTOAD_h

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "Geometry/GeometryGAr.h"
#include "ReadoutSimulation/ElectronDriftAlg.h"


namespace gar {

  namespace sdp{
    class EnergyDeposit;
  }
  
  namespace rosim{
    
    class ElectronDriftStandardAlgToad : public ElectronDriftAlg {
      
    public:
      
      ElectronDriftStandardAlgToad(CLHEP::HepRandomEngine      & engine,
                               fhicl::ParameterSet    const& pset);
      virtual ~ElectronDriftStandardAlgToad();
      
      void DriftElectronsToReadout(gar::sdp::EnergyDeposit       const& dep,
                                   gar::rosim::ElectronDriftInfo      & driftInfo);
      
    private:
      
      int                                 fElectronsPerCluster; ///< Number of electrons to drift in a cluster
      size_t                              fMinClusters;         ///< Minimum number of clusters for diffusion integral
      gar::detinfo::ElecClock             fClock;               ///< electronics clock
      const gar::detinfo::DetectorClocks* fTime;                ///< electronics clock
      const gar::geo::GeometryCore*       fGeo;                 ///< Geometry
      
    };
    
  }
} // gar



#endif /* ElectronDriftStandardAlg_h */
