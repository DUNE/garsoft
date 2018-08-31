////////////////////////////////////////////////////////////////////////
//
//  ECALProperties
//
////////////////////////////////////////////////////////////////////////
// Framework includes

// C++ language includes
#include <cmath>
#include <iostream>

// GArSoft includes
#include "DetectorInfo/ECALPropertiesStandard.h"
#include "CoreUtils/ProviderUtil.h" // gar::IgnorableProviderConfigKeys()

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
gar::detinfo::ECALPropertiesStandard::ECALPropertiesStandard()
: fIsConfigured(false)
{
}

//-----------------------------------------------
gar::detinfo::ECALPropertiesStandard::ECALPropertiesStandard(fhicl::ParameterSet   const& pset,
                                                      std::set<std::string>        ignore_params /* = {} */)
: ECALPropertiesStandard()
{
  this->Configure(pset, ignore_params);
}

//------------------------------------------------
bool gar::detinfo::ECALPropertiesStandard::Configure(fhicl::ParameterSet   const& pset,
                                                    std::set<std::string>        ignore_params /* = {} */)
{
  std::set<std::string> ignorable_keys = gar::IgnorableProviderConfigKeys();
  ignorable_keys.insert(ignore_params.begin(), ignore_params.end());

  // validation happens here:
  fhicl::Table<Configuration_t> config_table { pset, gar::IgnorableProviderConfigKeys() };
  Configuration_t const& config = config_table();

  SetEffectivePixel      (config.EffectivePixel()      );
  SetLightYield         (config.LightYield()         );
  SetSiPMGain           (config.SiPMGain()           );
  fIsConfigured = true;

  return true;
}

//------------------------------------------------
bool gar::detinfo::ECALPropertiesStandard::Update(uint64_t ts)
{
  if (ts == 0) return false;

  return true;
}
