#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TVectorT.h"
#include "ReconstructionDataProducts/Hit.h"

using namespace art;
using namespace std;

//  ROOT script using gallery to dump garsoft hits
// Tom Junk, Fermilab, Dec. 2018

// arguments:  filename -- input file, larsoft formatted
// ievstart:  file index of starting event (staring at zero)
// nev:       number of events to dump
// tickmin, tickmax -- to truncate the part of the event run on.  Set to big and small numbers for no truncation.
// inputtag: use "RecoProc" for a standard MC job made with recojob.fcl


void hitdump(std::string const& filename="reco.root", 
             size_t ievstart=0, 
	     size_t nev=1,
             std::string const& inputtag="hit")
{


  size_t evcounter=0;

  InputTag hit_tag{ inputtag };
  // Create a vector of length 1, containing the input filename.
  vector<string> filenames(1, filename);

  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    auto const& evaux = ev.eventAuxiliary();
    if (evcounter >= ievstart && evcounter < ievstart+nev )
      {
	auto const& hits = *ev.getValidHandle<vector<gar::rec::Hit>>(hit_tag);
	if (!hits.empty())
	  {
	    std::cout << "Dumping hits for Run: " << evaux.run() << " " << " Subrun: " << evaux.subRun() << " Event: " << evaux.event() << std::endl;
	    const size_t nhits = hits.size();
	    size_t nchans=0;
	    for (size_t i=0; i<nhits; ++i)
	      {
		std::cout << i << " " << hits[i].Position()[0] << " " << hits[i].Position()[1] << " " << hits[i].Position()[2] << " " << 
		  hits[i].Signal() << " " << hits[i].RMS() << std::endl; 
	      }

	  }
      }
    ++evcounter;
  }
}
