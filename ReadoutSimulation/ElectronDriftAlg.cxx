//
//  ElectronDriftAlg.cxx
//  garsoft-mrb
//
//  Created by Brian Rebel on 11/18/16.
//

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "ReadoutSimulation/ElectronDriftAlg.h"
#include "CoreUtils/ServiceUtil.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

namespace gar {
  namespace rosim {
   
    //--------------------------------------------------------------------------
    ElectronDriftInfo::ElectronDriftInfo()
    {
      return;
    }
    
    //--------------------------------------------------------------------------
    void ElectronDriftInfo::Reset(std::vector<double> & xPos,
                                  std::vector<double> & yPos,
                                  std::vector<double> & zPos,
                                  std::vector<double> & time,
                                  std::vector<int   > & size)
    {
      fClusterXPos.clear();
      fClusterYPos.clear();
      fClusterZPos.clear();
      fClusterTime.clear();
      fClusterSize.clear();
      
      fClusterXPos.swap(xPos);
      fClusterYPos.swap(yPos);
      fClusterZPos.swap(zPos);
      fClusterTime.swap(time);
      fClusterSize.swap(size);
      
      if(fClusterXPos.size() != fClusterYPos.size() ||
         fClusterXPos.size() != fClusterZPos.size() ||
         fClusterXPos.size() != fClusterTime.size() ||
         fClusterXPos.size() != fClusterSize.size() )
        throw cet::exception("ElectronDriftInfo")
        << "ElectronDriftInfo vector sizes are not consistent: "
        << "\n\tX:"
        << fClusterXPos.size()
        << "\n\tY:"
        << fClusterYPos.size()
        << "\n\tZ:"
        << fClusterZPos.size()
        << "\n\tT:"
        << fClusterTime.size()
        << "\n\tSize:"
        << fClusterSize.size();
        
      return;
    }

    
    //--------------------------------------------------------------------------
    // The parameter set is in case a derived algorithm needs to be configured
    // beyond what is available from either DetectorPropertiesService or
    // GArPropertiesService.
    ElectronDriftAlg::ElectronDriftAlg(CLHEP::HepRandomEngine      & engine,
                                       fhicl::ParameterSet    const& /*pset*/)
    : fEngine(engine)
    {
      // get the corrections and constants from the necessary places
      auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      auto garProp = gar::providerFrom<detinfo::GArPropertiesService>();
      
      fDriftVelocity         = detProp->DriftVelocity(detProp->Efield(),
                                                      detProp->Temperature());
      fLifetimeCorrection    = -1000. * detProp->ElectronLifetime();
      fLongitudinalDiffusion = garProp->LongitudinalDiffusion();
      fTransverseDiffusion   = garProp->TransverseDiffusion();
      fFanoFactor            = garProp->FanoFactor();
      
      fInverseVelocity       = 1. / fDriftVelocity;

      // fLongitudinalDiffusion and fTransverseDiffusion are in microns/sqrt(cm)

      fLongDiffConst         = fLongitudinalDiffusion * 1E-4;  // in sqrt(cm)
      fTransDiffConst        = fTransverseDiffusion * 1E-4;   // in sqrt(cm)
      
      return;
    }
    
    //--------------------------------------------------------------------------
    ElectronDriftAlg::~ElectronDriftAlg()
    {
    }

    
    
  }
} // rosim
