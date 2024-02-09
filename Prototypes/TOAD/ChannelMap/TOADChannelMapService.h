#ifndef TOADChannelMapService_H
#define TOADChannelMapService_H

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "TOADChannelMapSP.h"

namespace dune {
  class TOADChannelMapService;
}

class dune::TOADChannelMapService {
  public:

  TOADChannelMapService(fhicl::ParameterSet const& pset);
  TOADChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  dune::TOADChannelMapSP::HDChanInfo_t GetChanInfoFromTOADElements(
   unsigned int toad_index) const;

  dune::TOADChannelMapSP::HDChanInfo_t GetChanInfoFromOfflChan(unsigned int offlchan) const;

  private:

  dune::TOADChannelMapSP fHDChanMap;

};

DECLARE_ART_SERVICE(dune::TOADChannelMapService, LEGACY)

#endif
