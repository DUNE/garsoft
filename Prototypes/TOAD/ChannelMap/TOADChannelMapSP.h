#ifndef TOADChannelMapSP_H
#define TOADChannelMapSP_H

#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>

namespace dune {
  class TOADChannelMapSP;
}

class dune::TOADChannelMapSP {

public:

  typedef struct HDChanInfo {
    unsigned int offlchan;        // OROC Offline Channel Index for Garsoft (9984-19967)
    unsigned int padrow;           // Pad Row Number (0-95)
    unsigned int pad;          // Pad number (0-(Np-1)) Np = Pad row number
    unsigned int connector;             // Connector Number (1-468)
    unsigned int pin;            // Pin Number on Connector (0-22)
    unsigned int alice_fec;    // FEC number in ALICE Map (0-77)
    unsigned int alice_fec_chan;         // Channel number on FEC in ALICE Map (0-127)
    unsigned int alice_fec_connector;           // ALICE FEC Connector Number (0-5)
    unsigned int toad_fec;   // TOAD FEC number from ALICE Map (0-155)
    unsigned int toad_fec_connector;            // TOAD FEC Connector number (0-2)
    unsigned int toad_fec_new;            // TOAD FEC Final Slot Number (mapped to TOAD slot layout) (0-155)
    unsigned int toad_agg;        // TOAD Aggregator
    unsigned int toad_agg_row;    // TOAD Aggregator Row Connector Number (a=1, g=7) (1-7)
    unsigned int toad_agg_col; //TOAD Aggregator Column Connector (1-8)
    unsigned int toad_fec_chan; // TOAD FEC Channel Number (1-64)
    unsigned int toad_index; //TOAD Frame Channel Index
    bool valid;          // true if valid, false if not
  } HDChanInfo_t;

  TOADChannelMapSP();  // constructor

  // initialize:  read map from file

  void ReadMapFromFile(std::string &fullname);

  // TPC channel map accessors

  // Map instrumentation numbers to offline channel number.

  HDChanInfo_t GetChanInfoFromTOADElements(unsigned int toad_index) const;

  HDChanInfo_t GetChanInfoFromOfflChan(unsigned int offlchan) const;

private:

  const unsigned int fNChans = 9984;

  // map of toad_index to channel info struct

  std::unordered_map<unsigned int,   // toad_index
    HDChanInfo_t > DetToChanInfo;

  // map of chan info indexed by offline channel number

  std::unordered_map<unsigned int, HDChanInfo_t> OfflToChanInfo;

  //-----------------------------------------------

  void check_offline_channel(unsigned int offlineChannel) const
  {
  if (offlineChannel >= (fNChans*2) && offlineChannel <= 9983)
    {      
      //throw std::range_error("TOADChannelMapSP offline Channel out of range"); 
    }
  };

};


#endif
