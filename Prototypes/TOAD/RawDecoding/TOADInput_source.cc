#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "TOADInput.h"
#include "RawDataProducts/RawDigit.h"
#include "TOADFrameOverlay.hpp"
#include "Prototypes/TOAD/ChannelMap/TOADChannelMapService.h"

// maybe make our own timestamp data product if need be.
// #include "lardataobj/RawData/RDTimeStamp.h"

gar::TOADInputDetail::TOADInputDetail(
				      fhicl::ParameterSet const & ps,
				      art::ProductRegistryHelper & rh,
				      art::SourceHelper const & sh) 
  : pretend_module_name(ps.get<std::string>("raw_data_label", "daq")),
    fLogLevel(ps.get<int>("LogLevel", 0)),
    fClockFreqMHz(ps.get<double>("ClockFrequencyMHz", 50.0)),
    fHandleSequenceOption(ps.get<std::string>("HandleSequenceOption","ignore")),
    fTrnScale(ps.get<unsigned int>("TrnScale",10000)),
    pmaker(sh) {
    rh.reconstitutes<std::vector<gar::raw::RawDigit>, art::InEvent>(pretend_module_name);
		   // don't link to lardataobj just to get a RDTimeStamp definition.  If we need one, maybe make a GAr one
		   // rh.reconstitutes<raw::RDTimeStamp, art::InEvent>(pretend_module_name, "trigger");
}

  void gar::TOADInputDetail::readFile(
				      std::string const & filename, art::FileBlock*& fb) {

    // open the input file with the dunedaq class

    fRawDataFilePtr = std::make_unique<dunedaq::hdf5libs::HDF5RawDataFile>(filename);

    fUnprocessedEventRecordIDs = fRawDataFilePtr->get_all_trigger_record_ids();
    fLastEvent = 0;

    uint32_t run_number = fRawDataFilePtr->get_attribute<uint32_t>("run_number");
    MF_LOG_INFO("HDF5")
      << "HDF5 opened HDF file with run number " <<
      run_number  << " and " <<
      fUnprocessedEventRecordIDs.size() << " events";
    if (fLogLevel > 0)
      {
	for (const auto & e : fUnprocessedEventRecordIDs) MF_LOG_INFO("HDF5") << e.first << " " << e.second;
      }

    fb = new art::FileBlock(art::FileFormatVersion(1, "RawEvent2011"),
			    filename); 
  }

bool gar::TOADInputDetail::readNext(art::RunPrincipal const* const inR,
				    art::SubRunPrincipal const* const inSR,
				    art::RunPrincipal*& outR,
				    art::SubRunPrincipal*& outSR,
				    art::EventPrincipal*& outE) 
{
  
  // Establish default 'results'
  outR = 0;
  outSR = 0;
  outE = 0;

  if (fUnprocessedEventRecordIDs.empty()) {return false;}

  // take the first unprocessed event off the set
  //MF_LOG_INFO("HDF5") << "Here 1";
  auto nextEventRecordID_i = fUnprocessedEventRecordIDs.cbegin();
  auto nextEventRecordID = *nextEventRecordID_i;
  fUnprocessedEventRecordIDs.erase(nextEventRecordID_i);
  //MF_LOG_INFO("HDF5") << "Here 2";
 
  uint32_t run_id = fRawDataFilePtr->get_attribute<uint32_t>("run_number");

  // get trigger record header pointer

  auto trh = fRawDataFilePtr->get_trh_ptr(nextEventRecordID);

  //check that the run number in the trigger record header agrees with that in the file attribute
  //MF_LOG_INFO("HDF5") << "Here 3";

  auto run_id2 = trh->get_run_number();
  if (run_id != run_id2)
    {
      throw cet::exception("TOADInput_source.cc") << " Inconsistent run IDs in file attribute and trigger record header " << run_id << " " << run_id2;
    }

  // check that the trigger record number in the header agrees with what we are expecting

  auto trhtn = trh->get_trigger_number();
  if (nextEventRecordID.first != trhtn)
    {
      throw cet::exception("TOADInput_source.cc") << " Inconsistent trigger record numbers: " << nextEventRecordID.first << " " << trhtn;
    }
  
  // trigTimeStamp is NOT time but Clock-tick since the epoch.
  //MF_LOG_INFO("HDF5") << "Here 4";

  uint64_t trigTimeStamp = trh->get_trigger_timestamp();
  if (fLogLevel > 0)
    {
      std::cout << "TOADInput_source: Trigger Time Stamp: " << trigTimeStamp << std::endl;
    }

  // this reformats the timestamp -- figure out what you need for TOAD
  // uint64_t getTrigTime = formatTrigTimeStampWithFrequency(trigTimeStamp, fClockFreqMHz);
  // currently just a straight copy
  uint64_t getTrigTime = trigTimeStamp;
  if (fLogLevel > 0)
    {
      std::cout << "TOADInput_source: getTrigTime: " << getTrigTime << std::endl;
    }
  //MF_LOG_INFO("HDF5") << "Here 5";

  art::Timestamp artTrigStamp (getTrigTime);
  if (fLogLevel > 0)
    {
      std::cout << "TOADInput_source: artTrigStamp: " <<   artTrigStamp.value()<< std::endl;
    }

  //std::unique_ptr<raw::RDTimeStamp> rd_timestamp(
  //                                               new raw::RDTimeStamp(trigTimeStamp));

  // make new run if inR is 0 or if the run has changed
  if (inR == 0 || inR->run() != run_id) {
    outR = pmaker.makeRunPrincipal(run_id,artTrigStamp);
  }

  // make new subrun if inSR is 0 or if the subrun has changed
  art::SubRunID subrun_check(run_id, 1);
  if (inSR == 0 || subrun_check != inSR->subRunID()) {
    outSR = pmaker.makeSubRunPrincipal(run_id, 1, artTrigStamp);
  }

  // this is the trigger record ID:
  //MF_LOG_INFO("HDF5") << "Here 6";

  int event = nextEventRecordID.first;
  int sequenceno = nextEventRecordID.second;
  if ( fHandleSequenceOption == "increment" )
    {
      event = fLastEvent + 1;
      if (fLogLevel > 0)
	{
	  std::cout << "TOADInput_source: assigning event number via incrment: " << event << " orig: " << nextEventRecordID.first << " " << nextEventRecordID.second;
	}
    }
  else if ( fHandleSequenceOption == "shiftadd" )
    {
      event = event*fTrnScale + sequenceno;
	if (fLogLevel > 0)
	{
	  std::cout << "TOADInput_source: assigning event number via shfitadd: " << event << " orig: " << nextEventRecordID.first << " " << nextEventRecordID.second;
	}
    }
  else if ( fHandleSequenceOption == "ignore" )
    {
      if (fLogLevel > 0)
	{
	  std::cout << "TOADInput_source: assigning event number ignoring seqID: " << event << " orig: " << nextEventRecordID.first << " " << nextEventRecordID.second;
	}
    }
  else
    {
      throw cet::exception("TOADInput_source.cc") << "Ununderstood HandleSequenceOption: " << fHandleSequenceOption << " use: ignore, increment or shiftadd";
    }
  fLastEvent = event;

  outE = pmaker.makeEventPrincipal(run_id, 1, event, artTrigStamp);
  if (fLogLevel > 0)
    {
      std::cout << "TOADInput_source: Event Time Stamp :" << outE->time().value() << std::endl;
    }

  //MF_LOG_INFO("HDF5") << "Here 7";

  std::unique_ptr<std::vector<gar::raw::RawDigit>> rdcol = std::make_unique<std::vector<gar::raw::RawDigit>>();
  //MF_LOG_INFO("HDF5") << "Here 8";

  // put code here to unpack the data and put it in rdcol.
  // see the example at duneprototypes/Protodune/hd/RawDecoding/PDHDDataInterfaceWIB3_tool.cc in the method getFragmentsForEvent
  // instead of retrieving the raw file pointer from a service, it's fRawDataFilePtr as a private member of this class
  // use it instead of rf
  
  art::ServiceHandle<dune::TOADChannelMapService> channelMap;
  dunedaq::hdf5libs::HDF5RawDataFile::record_id_t rid = std::make_pair(event, sequenceno);
  auto sourceids = fRawDataFilePtr->get_source_ids(rid);
  //MF_LOG_INFO("HDF5") << "Here 9";
    
  for (const auto &source_id : sourceids)  
    {
      // only want detector readout data (i.e. not trigger info)
      if (source_id.subsystem != dunedaq::daqdataformats::SourceID::Subsystem::kDetectorReadout) continue;

      // look through the geo IDs and see if we are in the right crate
      auto gids = fRawDataFilePtr->get_geo_ids_for_source_id(rid, source_id);
      for (const auto &gid : gids)
        {
          if (fLogLevel > 1)
            {
              std::cout << "TOADInputSource Geoid: " << std::hex << gid << std::dec << std::endl;
            }
          uint16_t detid = 0xffff & gid;
          dunedaq::detdataformats::DetID::Subdetector detidenum = static_cast<dunedaq::detdataformats::DetID::Subdetector>(detid);
          auto subdetector_string = dunedaq::detdataformats::DetID::subdetector_to_string(detidenum);
	  if (fLogLevel > 1)
	    {
	      std::cout << "TOADInputSource subdetector string: " << subdetector_string << std::endl;
	      //std::cout << "TOADInputSource looking for subdet: " << " Put string here.  ND_LAr?" << std::endl;
	      std::cout << "TOADInputSource looking for subdet: " << "ND_GAr" << std::endl;
	    }
	  // I stopped copying code here ... Need to check the subdetector string and unpack the data.
          //MF_LOG_INFO("HDF5") << "Here 10";

	  if(subdetector_string == "ND_GAr")
            {
              //MF_LOG_INFO("HDF5") << "Here 11";

              auto frag = fRawDataFilePtr->get_frag_ptr(rid, source_id);
              //auto frag_hdr = frag->get_header();
              //auto frag_ts = frag->get_trigger_timestamp();
              size_t fhs = sizeof(dunedaq::daqdataformats::FragmentHeader);
              auto frag_size = frag->get_size() - fhs;
              MF_LOG_INFO("HDF5") << "frag_size "<<(int)frag_size<<" fhs "<<(int)fhs;
              //int frag_start = 0;
	      //int frag_count = 0;
              int frag_interval=0;
              gar::raw::ADCvector_t adc_vector;
              //MF_LOG_INFO("HDF5") << "Here 12";
              //int i_test = 0;
              ULong64_t time;
	      ULong64_t t0 =0;
	      //short t_test = 1;
//              for(frag_start = 0; frag_start<=(int)(frag_size); frag_start =+ frag_count) // frag_start<=number ; frag_start+=frag+interval
	      while(frag_interval<(int)(frag_size))
                {
                auto fdp = static_cast<void*>(static_cast<char*>(frag->get_data()) + frag_interval);
                auto toad_frame = reinterpret_cast<dunedaq::detdataformats::toad::TOADFrameOverlay*>(static_cast<char*>(fdp));
                size_t nsamples = toad_frame->n_samples;
                //MF_LOG_INFO("HDF5") << "frag_start "<< frag_start<<" n_samples "<< nsamples;
		adc_vector.push_back(nsamples);
		adc_vector.push_back(1);
		adc_vector.push_back(4);
		adc_vector.push_back(nsamples);
                for (size_t i_samp=0; i_samp<nsamples; i_samp++)
                {
                  adc_vector.push_back(toad_frame->get_samples(i_samp));
                  //MF_LOG_INFO("HDF5") << "samples"<<(int)(toad_frame->get_samples(i_samp));
                }
                auto chaninfo = channelMap->GetChanInfoFromTOADElements (toad_frame->fec);
                unsigned int offline_chan = chaninfo.offlchan;
                MF_LOG_INFO("HDF5") << "channel"<< (int)(toad_frame->fec) << "offlchan"<< offline_chan;
                //unsigned int offline_chan = 0;
                gar::raw::Compress_t compress = gar::raw::kZeroSuppression;
		if(adc_vector[0] == 0){
		  continue;
		}
		
		if(frag_interval == 0){
		  t0 = (ULong64_t)(toad_frame->tstmp);
		}
		time = (ULong64_t)(toad_frame->tstmp)-(ULong64_t)t0+1;
                MF_LOG_INFO("HDF5") << "timestamps "<< time << "t0"<<t0;
		//adc_vector[2] = 1;
		//t_test++;
                //t0 = 0;
		//time = t0;
                gar::raw::RawDigit rd(offline_chan, adc_vector.size(), adc_vector, compress, time);
                adc_vector.clear();
		(*rdcol).push_back(rd);
                //frag_start = frag_start + (int)(toad_frame ->n_bytes);
                frag_interval = frag_interval + (int)(toad_frame ->n_bytes);
		//MF_LOG_INFO("HDF5") << "frag_interval, frag_start_sum "<<frag_interval<<" "<<(frag_start+frag_interval);
                //MF_LOG_INFO("HDF5") << "Here 13 "<< frag_interval;
                //i_test++;

                }
            }
	}
    }
  MF_LOG_INFO("HDF5") << "Here 14";

  put_product_in_principal(std::move(rdcol), *outE, pretend_module_name,
                           "");
  MF_LOG_INFO("HDF5") << "Here 15";

  //put_product_in_principal(std::move(rd_timestamp), *outE, pretend_module_name,
  //                         "trigger");

  return true;
}

//typedef for shorthand
namespace gar {
  using TOADInputSource = art::Source<TOADInputDetail>;
}


DEFINE_ART_INPUT_SOURCE(gar::TOADInputSource)
