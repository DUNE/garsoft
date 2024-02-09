#include "TOADChannelMapSP.h"

#include <iostream>
#include <fstream>
#include <sstream>

// so far, nothing needs to be done in the constructor

dune::TOADChannelMapSP::TOADChannelMapSP()
{
}

void dune::TOADChannelMapSP::ReadMapFromFile(std::string &fullname)
{
  std::ifstream inFile(fullname, std::ios::in);
  std::string line;

  while (std::getline(inFile,line)) {
    std::stringstream linestream(line);
    HDChanInfo_t chanInfo;
    linestream 
      >> chanInfo.offlchan 
      >> chanInfo.padrow 
      >> chanInfo.pad
      >> chanInfo.connector 
      >> chanInfo.pin 
      >> chanInfo.alice_fec 
      >> chanInfo.alice_fec_chan 
      >> chanInfo.alice_fec_connector 
      >> chanInfo.toad_fec 
      >> chanInfo.toad_fec_connector 
      >> chanInfo.toad_fec_new 
      >> chanInfo.toad_agg
      >> chanInfo.toad_agg_row
      >> chanInfo.toad_agg_col
      >> chanInfo.toad_fec_chan
      >> chanInfo.toad_index; 

    chanInfo.valid = true;

    // fill maps.

    check_offline_channel(chanInfo.offlchan);

    DetToChanInfo[chanInfo.toad_index] = chanInfo;
    OfflToChanInfo[chanInfo.offlchan] = chanInfo;

  }
  inFile.close();

}

dune::TOADChannelMapSP::HDChanInfo_t dune::TOADChannelMapSP::GetChanInfoFromTOADElements(
    unsigned int toad_index) const {

  HDChanInfo_t badInfo = {};
  badInfo.valid = false;

  auto fm1 = DetToChanInfo.find(toad_index);
  if (fm1 == DetToChanInfo.end()) return badInfo;
  return fm1->second;
}


dune::TOADChannelMapSP::HDChanInfo_t dune::TOADChannelMapSP::GetChanInfoFromOfflChan(unsigned int offlineChannel) const {
  auto ci = OfflToChanInfo.find(offlineChannel);
  if (ci == OfflToChanInfo.end()) 
    {
      HDChanInfo_t badInfo = {};
      badInfo.valid = false;
      return badInfo;
    }
  return ci->second;
}
