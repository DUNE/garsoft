////////////////////////////////////////////////////////////////////////
// \file GArProperties.h
//
// \brief pure virtual base interface for run history
//
// \author jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_IRUNHISTORY_H
#define DETINFO_IRUNHISTORY_H

///General GArSoft Utilities
namespace gar {
  namespace detinfo {

    enum RunType_t {
      kUnknownRunType = 0,
      kProductionRun,
      kCommissioningRun,
      kTestRun,
      kPedestalRun,
      kCalibrationRun,
      kNRunType
    };

    class SubRun {
    public:
      virtual ~SubRun() = default;

      virtual uint64_t TStart() const = 0;
    };

    class RunHistory {
    public:
      virtual ~RunHistory() = default;

      virtual bool Update(uint64_t) = 0;

      virtual int RunNumber() const = 0;
      virtual int NSubruns() const = 0;
      virtual int RunType() const = 0;
      virtual std::string RunTypeAsString() const = 0;
      virtual uint64_t TStart() const = 0;
      virtual uint64_t TStop() const = 0;
      virtual uint64_t Duration() const = 0;

    }; // class RunHistory
  }    //namespace detinfo
} // gar
#endif // RunHistory_H
