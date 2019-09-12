////////////////////////////////////////////////////////////////////////
/// \file  ISCalcuationSeparate.cxx
/// \brief Interface to algorithm class for calculating ionization electrons
///        and scintillation photons using separate algorithms for each
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorProperties.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "ReadoutSimulation/ISCalculationSeparate.h"
#include "CoreUtils/ServiceUtil.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

namespace gar {
  namespace rosim{
    
    //----------------------------------------------------------------------------
    ISCalculationSeparate::ISCalculationSeparate(CLHEP::HepRandomEngine&)
    : ISCalculation()
    {
    }
    
    //----------------------------------------------------------------------------
    ISCalculationSeparate::~ISCalculationSeparate()
    {
    }
    
    //----------------------------------------------------------------------------
    void ISCalculationSeparate::Initialize()
    {
      const detinfo::DetectorProperties* detprop = gar::providerFrom<detinfo::DetectorPropertiesService>();
      
      //double density       = detprop->Density(detprop->Temperature());
      fEfield              = detprop->Efield();
      fGeVToElectrons      = gar::detinfo::kGeVToElectrons;
      
      return;
    }
    
    //----------------------------------------------------------------------------
    // fNumIonElectrons returns a value that is not corrected for life time effects
    void ISCalculationSeparate::Reset()
    {
      fEnergyDeposit   = 0.;
      fNumScintPhotons = 0.;
      fNumIonElectrons = 0.;
      
      return;
    }
    
    //----------------------------------------------------------------------------
    // fNumIonElectrons returns a value that is not corrected for life time effects
    void ISCalculationSeparate::CalculateIonizationAndScintillation(const gar::sdp::EnergyDeposit* dep)
    {
      fEnergyDeposit = dep->Energy();
      
      fNumIonElectrons = fGeVToElectrons * fEnergyDeposit;
      
      LOG_DEBUG("ISCalculationSeparate")
      << " Electrons produced for " << fEnergyDeposit * 1.e3
      << " MeV deposited: "       << fNumIonElectrons;
      
      // Now do the scintillation
      // TODO: fix this for gaseous argon
      fNumScintPhotons = std::numeric_limits<float>::min();
      
      return;
    }
    
  }// namespace
} // gar
