////////////////////////////////////////////////////////////////////////
// RunHistoryStandard.h
//
//  Detector info class for run history
//
// jpaley@fnal.gov
//
////////////////////////////////////////////////////////////////////////
#ifndef DETINFO_RUNHISTORY_H
#define DETINFO_RUNHISTORY_H

#include <map>
#include <string>
#include <vector>

#include "DetectorInfo/RunHistory.h"
#include "fhiclcpp/ParameterSet.h"

///General GArSoft Utilities
namespace gar {
  namespace detinfo {

    class SubRunStandard : public SubRun {
    public:
      SubRunStandard() : fTStart(0){};
      SubRunStandard(SubRunStandard const&) = delete;
      virtual ~SubRunStandard(){};

      virtual uint64_t TStart() const override { return fTStart; }
      void SetTStart(uint64_t t) { fTStart = t; }

    private:
      uint64_t fTStart;
    };

    class RunHistoryStandard : public RunHistory {
    public:
      RunHistoryStandard();
      RunHistoryStandard(int runnum);
      RunHistoryStandard(RunHistoryStandard const&) = delete;
      virtual ~RunHistoryStandard();

      virtual bool Update(uint64_t ts = 0) override;

      virtual int RunNumber() const override { return fRun; }
      virtual int NSubruns() const override { return fNSubruns; }
      virtual int RunType() const override { return fRunType; }
      virtual std::string RunTypeAsString() const override;
      virtual uint64_t TStart() const override { return fTStart; }
      virtual uint64_t TStop() const override { return fTStop; }
      virtual uint64_t Duration() const override { return fTStop - fTStart; }

      std::vector<std::string> Shifters() { return fShifter; }

      void SetNSubruns(int nsr) { fNSubruns = nsr; }
      void SetRunType(int rt) { fRunType = rt; }
      void SetDetId(int id) { fDetId = id; }
      void SetTStart(uint64_t t) { fTStart = t; }
      void SetTStop(uint64_t t) { fTStop = t; }
      void AddShifter(std::string sh) { fShifter.push_back(sh); }
      void SetShifters(std::vector<std::string> sh) { fShifter = sh; }
      void SetDetName(std::string dn) { fDetName = dn; }

    private:
    protected:
      int fRun;
      int fNSubruns;
      int fRunType;
      int fDetId;

      uint64_t fTStart;
      uint64_t fTStop;

      std::vector<std::string> fShifter;
      std::string fDetName;

      std::vector<SubRunStandard> fSubrun;

    }; // class RunHistoryStandard
  }    //namespace detinfo
} // gar

#endif // DETINFO_RUNHISTORY_H
