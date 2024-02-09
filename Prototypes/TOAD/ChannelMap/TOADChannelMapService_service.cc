#include "TOADChannelMapService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

dune::TOADChannelMapService::TOADChannelMapService(fhicl::ParameterSet const& pset) {

  std::string channelMapFile = pset.get<std::string>("FileName");

  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(channelMapFile, fullname);

  if (fullname.empty()) {
    std::cout << "Input file " << channelMapFile << " not found" << std::endl;
    throw cet::exception("File not found");
  }
  else
    std::cout << "TOAD Channel Map: Building TPC wiremap from file " << channelMapFile << std::endl;

  fHDChanMap.ReadMapFromFile(fullname);
}

dune::TOADChannelMapService::TOADChannelMapService(fhicl::ParameterSet const& pset, art::ActivityRegistry&) : TOADChannelMapService(pset) {
}

dune::TOADChannelMapSP::HDChanInfo_t dune::TOADChannelMapService::GetChanInfoFromTOADElements(
    unsigned int toad_index ) const {

  return fHDChanMap.GetChanInfoFromTOADElements(toad_index);
}


dune::TOADChannelMapSP::HDChanInfo_t dune::TOADChannelMapService::GetChanInfoFromOfflChan(unsigned int offlineChannel) const {

  return fHDChanMap.GetChanInfoFromOfflChan(offlineChannel);

}


DEFINE_ART_SERVICE(dune::TOADChannelMapService)
