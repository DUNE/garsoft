////////////////////////////////////////////////////////////////////////
//
//  GArProperties
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <cmath>
#include <iostream>

// GArSoft includes
#include "CoreUtils/ProviderUtil.h" // gar::IgnorableProviderConfigKeys()
#include "DetectorInfo/GArPropertiesStandard.h"

// ROOT includes
#include "TMath.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
gar::detinfo::GArPropertiesStandard::GArPropertiesStandard() : fIsConfigured(false) {}

//-----------------------------------------------
gar::detinfo::GArPropertiesStandard::GArPropertiesStandard(
  fhicl::ParameterSet const& pset,
  std::set<std::string> ignore_params /* = {} */)
  : GArPropertiesStandard()
{
  this->Configure(pset, ignore_params);
}

#if 0
//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Configure(fhicl::ParameterSet const& pset)
{
  this->SetRadiationLength  (pset.get< double >("RadiationLength" ));
  this->SetAtomicNumber     (pset.get< double >("AtomicNumber"));
  this->SetAtomicMass       (pset.get< double >("AtomicMass"));
  this->SetMeanExcitationEnergy (pset.get< double >("ExcitationEnergy"));
  fIsConfigured = true;
  return true;
}
#endif // 0

//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Configure(fhicl::ParameterSet const& pset,
                                                    std::set<std::string> ignore_params /* = {} */)
{
  std::set<std::string> ignorable_keys = gar::IgnorableProviderConfigKeys();
  ignorable_keys.insert(ignore_params.begin(), ignore_params.end());

  // validation happens here:
  fhicl::Table<Configuration_t> config_table{pset, gar::IgnorableProviderConfigKeys()};
  Configuration_t const& config = config_table();

  SetRadiationLength(config.RadiationLength());
  SetAtomicNumber(config.AtomicNumber());
  SetAtomicMass(config.AtomicMass());
  SetMeanExcitationEnergy(config.MeanExcitationEnergy());
  SetFanoFactor(config.FanoFactor());
  SetLongitudinalDiffusion(config.LongitudinalDiffusion());
  SetTransverseDiffusion(config.TransverseDiffusion());
  fIsConfigured = true;

  return true;
}

//------------------------------------------------
bool gar::detinfo::GArPropertiesStandard::Update(uint64_t ts)
{
  if (ts == 0) return false;

  return true;
}
