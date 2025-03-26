#include "messagefacility/MessageLogger/MessageLogger.h"

#include "DetectorInfo/DetectorClocksStandardGAr.h"

namespace gar {

  //-------------------------------------------------------------------------
  detinfo::DetectorClocksStandardGAr::DetectorClocksStandardGAr()
    : fConfigName(detinfo::kInheritConfigTypeMax, "")
    , fConfigValue(detinfo::kInheritConfigTypeMax, 0)
    , fTrigModuleName("")
    , fG4RefTime(detinfo::kDEFAULT_MC_CLOCK_T0)
    , fFramePeriod(detinfo::kDEFAULT_FRAME_PERIOD)
    , fTPCClock(0, detinfo::kDEFAULT_FRAME_PERIOD, detinfo::kDEFAULT_FREQUENCY_TPC)
    , fTriggerClock(0, detinfo::kDEFAULT_FRAME_PERIOD, detinfo::kDEFAULT_FREQUENCY_TRIGGER)
    , fExternalClock(0, detinfo::kDEFAULT_FRAME_PERIOD, detinfo::kDEFAULT_FREQUENCY_EXTERNAL)
    , fTriggerOffsetTPC(detinfo::kDEFAULT_TRIG_OFFSET_TPC)
    , fTriggerTime(0)
    , fBeamGateTime(0)
  {

    fConfigName.at(detinfo::kG4RefTime) = "G4RefTime";
    fConfigName.at(detinfo::kTriggerOffsetTPC) = "TriggerOffsetTPC";
    fConfigName.at(detinfo::kFramePeriod) = "FramePeriod";
    fConfigName.at(detinfo::kClockSpeedTPC) = "ClockSpeedTPC";
    fConfigName.at(detinfo::kClockSpeedTrigger) = "ClockSpeedTrigger";
    fConfigName.at(detinfo::kClockSpeedExternal) = "ClockSpeedExternal";
    fConfigName.at(detinfo::kDefaultTrigTime) = "DefaultTrigTime";
    fConfigName.at(detinfo::kDefaultBeamTime) = "DefaultBeamTime";
    fConfigName.at(detinfo::kDefaultSpillLength) = "DefaultSpillLength";

    fInheritClockConfig = false;
  }

  //-------------------------------------------------------------------------
  detinfo::DetectorClocksStandardGAr::DetectorClocksStandardGAr(fhicl::ParameterSet const& pset)
    : DetectorClocksStandardGAr()
  {
    // In a constructor, the version of virtual method that is called
    // is always the one specific of the class being constructed
    // (the one mentioned in the name of the constructor itself).
    // For clarity, we explicitly show that:
    DetectorClocksStandardGAr::Configure(pset);
  }

  //------------------------------------------------------------------
  bool detinfo::DetectorClocksStandardGAr::Update(uint64_t)
  {
    return true;
  }

  //------------------------------------------------------------------
  bool detinfo::DetectorClocksStandardGAr::Configure(fhicl::ParameterSet const& pset)
  {

    // Read fcl parameters
    fTrigModuleName = pset.get<std::string>("TrigModuleName");
    fInheritClockConfig = pset.get<bool>("InheritClockConfig");
    fConfigValue.at(kG4RefTime) = pset.get<double>(fConfigName.at(kG4RefTime).c_str());
    fConfigValue.at(kFramePeriod) = pset.get<double>(fConfigName.at(kFramePeriod).c_str());
    fConfigValue.at(kTriggerOffsetTPC) =
      pset.get<double>(fConfigName.at(kTriggerOffsetTPC).c_str());
    fConfigValue.at(kClockSpeedTPC) = pset.get<double>(fConfigName.at(kClockSpeedTPC).c_str());
    fConfigValue.at(kClockSpeedTrigger) =
      pset.get<double>(fConfigName.at(kClockSpeedTrigger).c_str());
    fConfigValue.at(kClockSpeedExternal) =
      pset.get<double>(fConfigName.at(kClockSpeedExternal).c_str());
    fConfigValue.at(kDefaultTrigTime) = pset.get<double>(fConfigName.at(kDefaultTrigTime).c_str());
    fConfigValue.at(kDefaultBeamTime) = pset.get<double>(fConfigName.at(kDefaultBeamTime).c_str());
    fConfigValue.at(kDefaultSpillLength) =
      pset.get<double>(fConfigName.at(kDefaultSpillLength).c_str());

    // Save fcl parameters in a container to check for inheritance
    fBeamGateTime = fTriggerTime = 0;

    ApplyParams();

    return true;
  }

  //-----------------------------------
  void detinfo::DetectorClocksStandardGAr::ApplyParams()
  //-----------------------------------
  {

    fG4RefTime = fConfigValue.at(kG4RefTime);
    fFramePeriod = fConfigValue.at(kFramePeriod);
    fTriggerOffsetTPC = fConfigValue.at(kTriggerOffsetTPC);

    fTPCClock = detinfo::ElecClock(fTriggerTime, fFramePeriod, fConfigValue.at(kClockSpeedTPC));
    fTriggerClock =
      detinfo::ElecClock(fTriggerTime, fFramePeriod, fConfigValue.at(kClockSpeedTrigger));
  }

  //------------------------------------------------------------------------
  bool detinfo::DetectorClocksStandardGAr::IsRightConfig(const fhicl::ParameterSet& ps) const
  //------------------------------------------------------------------------
  {
    std::string s;
    double d;

    bool result = !ps.get_if_present("module_label", s);
    for (size_t i = 0; result && i < kInheritConfigTypeMax; ++i)

      result = result && ps.get_if_present(fConfigName.at(i).c_str(), d);

    return result;
  }

  //-----------------------------------------
  void detinfo::DetectorClocksStandardGAr::debugReport() const
  //-----------------------------------------
  {
    MF_LOG_VERBATIM("DetectorClocksStandardGAr") << "fConfigValues contents: ";

    for (size_t i = 0; i < kInheritConfigTypeMax; ++i)
      MF_LOG_VERBATIM("DetectorClocksStandardGAr")
        << "    " << fConfigName.at(i).c_str() << " ... " << fConfigValue.at(i);

    MF_LOG_VERBATIM("DetectorClocksStandardGAr")
      << "Trigger  time @ " << fTriggerTime << "\nBeamGate time @ " << fBeamGateTime
      << "\nTrigOffsetTPC @ " << TriggerOffsetTPC() << "\nG4RefTime     @ " << fG4RefTime
      << "\nTPC     Freq. @ " << fTPCClock.Frequency() << "\nTrigger Freq. @ "
      << fTriggerClock.Frequency() << "\nExternal Freq. @ " << fExternalClock.Frequency()
      << "\nTPC start tick [tdc]             : " << TPCTick2TDC(0)
      << "\nTPC start tick from trigger [ns] : " << TPCTick2TrigTime(0)
      << "\nTPC start tick from beam    [ns] : " << TPCTick2BeamTime(0)
      << "\nTPC tdc=0 in tick     : " << TPCTDC2Tick(0)
      << "\nTPC G4 time 0 in tick : " << TPCG4Time2Tick(0);
  }

}
