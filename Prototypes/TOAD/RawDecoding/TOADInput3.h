// Input source for the second HDF5 file format version from the DAQ people
// using the DUNE-DAQ HDF5RawInputFile class

#ifndef TOADInput3_h
#define TOADInput3_h
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/Core/fwd.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "dunedaqhdf5utils3/HDF5RawDataFile.hpp"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  //Forward declare the class
  class TOADInput3Detail;
}

class gar::TOADInput3Detail {
public:
  TOADInput3Detail(fhicl::ParameterSet const& ps,
                   art::ProductRegistryHelper& rh,
                   art::SourceHelper const& sh);

  void readFile(std::string const& filename, art::FileBlock*& fb);

  bool readNext(art::RunPrincipal const* const inR,
                art::SubRunPrincipal const* const inSR,
                art::RunPrincipal*& outR,
                art::SubRunPrincipal*& outSR,
                art::EventPrincipal*& outE);

  void closeCurrentFile() { fRawDataFilePtr.reset(); };

private:
  dunedaq::hdf5libs::HDF5RawDataFile::record_id_set fUnprocessedEventRecordIDs;
  std::string pretend_module_name;
  int fLogLevel;
  double fClockFreqMHz; // clock frequency in MHz -- used to unpack trigger timestamps for the event
  std::string fHandleSequenceOption; // to steer what to do with trigger record sequence numbers
  unsigned int fTrnScale; // in case we are doing shiftadd, this the scale factor on trig number
  art::SourceHelper const& pmaker;

  std::unique_ptr<dunedaq::hdf5libs::HDF5RawDataFile> fRawDataFilePtr;
  int fLastEvent;
};
#endif
