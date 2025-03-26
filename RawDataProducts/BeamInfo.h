/////////////////////////////////////////////////////////////////////
// $Id: BeamInfo.h,v 1.3 2010/03/26 19:36:42 brebel Exp $
//  information about the neutrino beam
//  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////
#ifndef RAWDATA_BEAMINFO_H
#define RAWDATA_BEAMINFO_H

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace gar {
  namespace raw {

    class BeamInfo {

    public:
      BeamInfo();
      ~BeamInfo();

    private:
      double tor101;
      double tortgt;
      double trtgtd;
      long long int t_ms;
      uint8_t fRecordType;
      uint32_t fSeconds; // GPS clock. Since Jan 1, 2012.
      uint16_t fMilliSeconds;
      uint16_t fNumberOfDevices;
      std::map<std::string, std::vector<double>> fDataMap;

#ifndef __GCCXML__

    public:
      BeamInfo(double tor101, double tortgt, double trtgtd, long long int t_ms);

      double get_tor101() const;
      double get_tortgt() const;
      double get_trtgtd() const;
      long long int get_t_ms() const;

      void SetTOR101(double val);
      void SetTORTGT(double val);
      void SetTRTGTD(double val);
      void SetT_MS(long long int val);

      void SetRecordType(uint8_t val) { fRecordType = val; };
      void SetSeconds(uint32_t val) { fSeconds = val; };
      void SetMilliSeconds(uint16_t val) { fMilliSeconds = val; };
      void SetNumberOfDevices(uint16_t val) { fNumberOfDevices = val; };

      void Set(std::string device_name, double val);
      void Set(std::string device_name, std::vector<double> val);

      uint8_t GetRecordType() const { return fRecordType; };
      uint32_t GetSeconds() const { return fSeconds; };
      uint16_t GetMilliSeconds() const { return fMilliSeconds; };
      uint16_t GetNumberOfDevices() const { return fNumberOfDevices; };

      std::vector<double> Get(std::string device_name);
      std::map<std::string, std::vector<double>> GetDataMap() const { return fDataMap; };

      friend std::ostream& operator<<(std::ostream&, const BeamInfo&);

#endif
    };

  }
} // gar

////////////////////////////////////////////////////////////////////////
#endif // RAWDATA_BEAMINFO_H
