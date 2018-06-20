////////////////////////////////////////////////////////////////////////
// Class:       CompressedHitFinder
// Plugin Type: producer (art v2_11_02)
// File:        CompressedHitFinder_module.cc
//
// Generated at Wed Jun 13 15:27:13 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Report the compressed raw digit blocks as hits -- just a threshold
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"
#include "ReconstructionDataProducts/Hit.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorClocksService.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include <memory>

#include "TMath.h"


namespace gar {
  namespace rec {

    class CompressedHitFinder : public art::EDProducer {
    public:
      explicit CompressedHitFinder(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      CompressedHitFinder(CompressedHitFinder const &) = delete;
      CompressedHitFinder(CompressedHitFinder &&) = delete;
      CompressedHitFinder & operator = (CompressedHitFinder const &) = delete;
      CompressedHitFinder & operator = (CompressedHitFinder &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      // Declare member data here.

      int fADCThreshold;   ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
      int fTicksBefore;    ///< zero-suppression ticks before
      int fTicksAfter;     ///< zero-suppression ticks after
      std::string fRawDigitLabel;  ///< label to find the right raw digits
      const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
      const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock
    };


    CompressedHitFinder::CompressedHitFinder(fhicl::ParameterSet const & p)
    // :
    {
      fADCThreshold = p.get<int>("ADCThreshold",5);
      fTicksBefore  = p.get<int>("TicksBefore",5);
      fTicksAfter   = p.get<int>("TicksAfter",5);
      fRawDigitLabel = p.get<std::string>("RawDigitLabel","daq");
      fTime    = gar::providerFrom<detinfo::DetectorClocksService>();
      fGeo     = gar::providerFrom<geo::Geometry>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      produces< std::vector<rec::Hit> >();
    }

    void CompressedHitFinder::produce(art::Event & e)
    {
      // create an emtpy output hit collection -- to add to.
      std::unique_ptr<std::vector<Hit> > hitCol (new std::vector<Hit> );
      // the input raw digits
      auto rdCol = e.getValidHandle< std::vector<raw::RawDigit> >(fRawDigitLabel);

      for (size_t ird = 0; ird < rdCol->size(); ++ ird)
	{
	  auto const& rd = (*rdCol)[ird];
	  auto channel = rd.Channel();
	  raw::ADCvector_t adc = rd.ADCs();
	  if (rd.Compression() == raw::kNone)
	    {
	      raw::ZeroSuppression(adc,fADCThreshold,fTicksBefore,fTicksAfter);
	    }
	  else if (rd.Compression() == raw::kZeroSuppression)
	    {
	      // we already have what we want
	    }
	  else
	    {
	      LOG_WARNING("CompressedHitFinder") << " Ununderstood compression mode: " << rd.Compression() << " Not making hits.";
	      e.put(std::move(hitCol));
	      return;
	    }
	  // use the format of the compressed raw digits -- a table of contents of the number of blocks, then all the block sizes, and then all the
	  // block start locations
	  if (adc.size() < 2)
	    {
	      e.put(std::move(hitCol));
	      return;
	    }

	  // walk through the zero-suppressed raw digits and call each block a hit

	  int nblocks = adc[1];
	  int zerosuppressedindex = nblocks*2 + 2;
	  float pos[3] = {0,0,0};
          fGeo->ChannelToPosition(channel, pos);
          float chanposx = pos[0];

	  for (int i=0; i<nblocks; ++i)
	    {
	      float hitSig = 0;
	      float hitTime = 0;
	      float hitSumSq = 0;
	      float hitRMS = 0;
	      unsigned int begT = adc[2+i];
	      int blocksize = adc[2+nblocks+i];
	      if (blocksize<1)
		{
		  throw cet::exception("CompressedHitFinder") << "Negative or zero block size in compressed data.";
		}
	      unsigned int endT = begT + blocksize;
	      for(int j = 0; j < blocksize; ++j)  // loop over time samples in each block
		{
		  int t = adc[2+i]+j;
		  int a = adc[zerosuppressedindex];
		  zerosuppressedindex++;
		  hitSig   += a;
		  hitTime  += a*t;
		  hitSumSq += a*t*t;
		}
	      if (hitSig > 0)  // otherwise leave the values at zero
		{
		  hitTime /= hitSig;
		  hitRMS = TMath::Sqrt(hitSumSq/hitSig - hitTime*hitTime);
		}
	      else
		{
		  hitTime = 0.5*(begT + endT);
		  hitRMS = 0;
		}
	      float driftdistance = fDetProp->DriftVelocity() * fTime->TPCTick2Time(hitTime);
	      if (chanposx < 0)
		{
		  pos[0] = chanposx + driftdistance;
		}
	      else
		{
		  pos[0] = chanposx - driftdistance;
		}

	      if (hitSig < 0)
		{
		  LOG_WARNING("CompressedHitFinder") << "Negative Signal in hit finder" << std::endl;
		}
	      hitCol->emplace_back(channel,
				   hitSig,
				   pos,
				   begT,
				   endT,
				   hitTime,
				   hitRMS);
	    }
	}
      e.put(std::move(hitCol));

    }

    DEFINE_ART_MODULE(CompressedHitFinder)

  } // namespace rec
} // namespace gar