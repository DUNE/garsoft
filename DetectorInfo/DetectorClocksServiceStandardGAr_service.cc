#include "DetectorInfo/DetectorClocksServiceStandardGAr.h"

#include "TFile.h"
#include "art_root_io/RootDB/SQLite3Wrapper.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------------------------------------------------
gar::detinfo::DetectorClocksServiceStandardGAr::DetectorClocksServiceStandardGAr(fhicl::ParameterSet   const& pset,
                                                                           ::art::ActivityRegistry      &reg)
: fClocks(std::make_unique<detinfo::DetectorClocksStandardGAr>(pset))
{
  
  reg.sPreProcessEvent.watch (this, &DetectorClocksServiceStandardGAr::preProcessEvent);
  reg.sPostOpenFile.watch    (this, &DetectorClocksServiceStandardGAr::postOpenFile);
  reg.sPreBeginRun.watch     (this, &DetectorClocksServiceStandardGAr::preBeginRun);
  
}

//------------------------------------------------------------------
void gar::detinfo::DetectorClocksServiceStandardGAr::reconfigure(fhicl::ParameterSet const& pset)
//------------------------------------------------------------------
{
  fClocks->Configure(pset);
  
}

//------------------------------------------------------------
void gar::detinfo::DetectorClocksServiceStandardGAr::preProcessEvent(::art::Event const& evt, art::ScheduleContext)
//------------------------------------------------------------
{
  auto trig_handle = evt.getHandle<std::vector<raw::Trigger> >(fClocks->TrigModuleName());

  std::vector<std::string> cfgNames = fClocks->ConfigNames();
  std::vector<double> cfgValues = fClocks->ConfigValues();

  fClocks->SetSpillLength( cfgValues.at(detinfo::kDefaultSpillLength) );
  
  if(!trig_handle.isValid() || trig_handle->empty()) {
      // Trigger simulation has not run yet!
    fClocks->SetTriggerTime(cfgValues.at(detinfo::kDefaultTrigTime),
                            cfgValues.at(detinfo::kDefaultBeamTime) );
    return;
  }
  
  if(trig_handle->size()>1)
    
    throw cet::exception("DetectorClocksServiceStandardGAr::preProcessEvent")
    << "Found "
    << trig_handle->size()
    << " triggers (only 1 trigger/event supported)\n";
  
  const ::art::Ptr<raw::Trigger> trig_ptr(trig_handle,0);
  
  fClocks->SetTriggerTime(trig_ptr->TriggerTime(),
                          trig_ptr->BeamGateTime() );

  return;
}

//------------------------------------------------------
void gar::detinfo::DetectorClocksServiceStandardGAr::preBeginRun(::art::Run const& ) // run)
//------------------------------------------------------
{
  // run number is unsigned so clang says this statement cannot happen
  //if (run.id().run() < 0) return;
  
  fClocks->ApplyParams();
}


//---------------------------------------------------------------
void gar::detinfo::DetectorClocksServiceStandardGAr::postOpenFile(std::string const& filename)
//---------------------------------------------------------------
{
  
  // Method inheriting from DetectorProperties
  
  if(!fClocks->InheritClockConfig()) return;
  
  // The only way to access art service metadata from the input file
  // is to open it as a separate TFile object.  Do that now.
  
  if(!filename.empty()) {
    
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if(file != 0 && !file->IsZombie() && file->IsOpen()) {
      
      std::vector<std::string> cfgName = fClocks->ConfigNames();
      std::vector<double> cfgValue = fClocks->ConfigValues();
      
      // Open the sqlite datatabase.
      
      ::art::SQLite3Wrapper sqliteDB(file, "RootFileDB");
      
      // Loop over all stored ParameterSets.
      
      std::vector<size_t> config_count(detinfo::kInheritConfigTypeMax,0);
      std::vector<double> config_value(detinfo::kInheritConfigTypeMax,0);
      
      sqlite3_stmt * stmt = 0;
      sqlite3_prepare_v2(sqliteDB, "SELECT PSetBlob from ParameterSets;", -1, &stmt, NULL);
      
      while (sqlite3_step(stmt) == SQLITE_ROW) {
        
        //fhicl::ParameterSet ps;
        //fhicl::make_ParameterSet(reinterpret_cast<char const *>(sqlite3_column_text(stmt, 0)), ps);
	auto ps = fhicl::ParameterSet::make(reinterpret_cast<char const *>(sqlite3_column_text(stmt, 0)));

        if(!fClocks->IsRightConfig(ps)) continue;
        
        for(size_t i=0; i<detinfo::kInheritConfigTypeMax; ++i) {
          
          double value_from_file = ps.get<double>(cfgName.at(i).c_str());
          
          if(!(config_count.at(i)))
            
            config_value.at(i) = value_from_file;
          
          else if(config_value.at(i) != value_from_file)
            
            throw cet::exception(__FUNCTION__)
            << Form("\033[95mFound historical value disagreement for %s ... %g != %g",
                    cfgName.at(i).c_str(),
                    config_value.at(i),
                    value_from_file)
            << "\033[00m";
          config_count.at(i) +=1;
          
        }
        
      }
      
      // Override parameters
      
      for(size_t i=0; i<detinfo::kInheritConfigTypeMax; ++i)
        
        if(config_count.at(i) && cfgValue.at(i) != config_value.at(i)) {
          
          MF_LOG_INFO("DetectorClocksServiceStandardGAr")
          << Form("\033[93mOverriding configuration parameter %s ... %g (fcl) => %g (data file)\033[00m",
                  cfgName.at(i).c_str(),
                  cfgValue.at(i),
                  config_value.at(i));
          
          fClocks->SetConfigValue(i,config_value.at(i));
          
        }
    }
    
    // Close file.
    if(file != 0) {
      if(file->IsOpen())
        file->Close();
      delete file;
    }
  }
  
  // Reset parameters
  fClocks->ApplyParams();
  
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(gar::detinfo::DetectorClocksServiceStandardGAr, gar::detinfo::DetectorClocksServiceGAr)

