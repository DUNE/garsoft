/////////////////////////////////////////////////////////////////////
// \version $Id: BeamInfo.cxx,v 1.3 2010/03/26 19:36:42 brebel Exp $
// \file    BeamInfo class
// \author  kinga.partyka@yale.edu
////////////////////////////////////////////////////////////////////

#include "RawDataProducts/BeamInfo.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <time.h>

namespace gar {
  namespace raw {

    //-------------------------------------------------------------------------
    BeamInfo::BeamInfo()
      : tor101(0.)
      , tortgt(0.)
      , trtgtd(0.)
      , t_ms(0)
      , fRecordType(0)
      , fSeconds(0)
      , fMilliSeconds(0)
      , fNumberOfDevices(0)
    {}

    //-------------------------------------------------------------------------
    BeamInfo::~BeamInfo() {}

    //-------------------------------------------------------------------------
    BeamInfo::BeamInfo(double vtor101, double vtortgt, double vtrtgtd, long long int vt_ms)
      : tor101(vtor101)
      , tortgt(vtortgt)
      , trtgtd(vtrtgtd)
      , t_ms(vt_ms)
      , fRecordType(0)
      , fSeconds(0)
      , fMilliSeconds(0)
      , fNumberOfDevices(0)
    {}

    //-------------------------------------------------------------------------
    double BeamInfo::get_tor101() const
    {
      return tor101;
    }

    //-------------------------------------------------------------------------
    double BeamInfo::get_tortgt() const
    {
      return tortgt;
    }

    //-------------------------------------------------------------------------
    double BeamInfo::get_trtgtd() const
    {
      return trtgtd;
    }

    //-------------------------------------------------------------------------
    long long int BeamInfo::get_t_ms() const
    {
      return t_ms;
    }

    //-------------------------------------------------------------------------
    void BeamInfo::SetTOR101(double val)
    {
      tor101 = val;
      Set("E:TOR101", val);
    }

    //-------------------------------------------------------------------------
    void BeamInfo::SetTORTGT(double val)
    {
      tortgt = val;
      Set("E:TORTGT", val);
    }

    //-------------------------------------------------------------------------
    void BeamInfo::SetTRTGTD(double val)
    {
      trtgtd = val;
      Set("E:TRTGTD", val);
    }

    //-------------------------------------------------------------------------
    void BeamInfo::SetT_MS(long long int val)
    {
      t_ms = val;
    }

    //-------------------------------------------------------------------------
    void BeamInfo::Set(std::string device, double val)
    {
      std::vector<double> vec;
      vec.push_back(val);
      Set(device, vec);
    }

    //-------------------------------------------------------------------------
    void BeamInfo::Set(std::string device, std::vector<double> val)
    {
      if (fDataMap.find(device) != fDataMap.end()) {
        //device already listed
        return;
      };
      std::pair<std::string, std::vector<double>> p(device, val);
      fDataMap.insert(p);
    }

    //-------------------------------------------------------------------------
    //output operator
    std::ostream& operator<<(std::ostream& os, const raw::BeamInfo& o)
    {

      os << "Record type:" << o.GetRecordType() << std::endl;
      ;
      os << "Timestamp: " << o.GetSeconds() << "\t" << o.GetMilliSeconds() << std::endl;
      os << "Number of Devices: " << o.GetNumberOfDevices() << std::endl;

      const std::map<std::string, std::vector<double>> dm = o.GetDataMap();
      std::map<std::string, std::vector<double>>::const_iterator it = dm.begin();
      while (it != dm.end()) {
        os << it->first << ": ";
        for (size_t i = 0; i < it->second.size(); i++)
          os << it->second[i] << ", ";
        os << std::endl;
        ++it;
      }

      return os;
    }
  } // namespace
} // gar
