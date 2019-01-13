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

      int fADCThreshold;        ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
      int fTicksBefore;         ///< zero-suppression ticks before
      int fTicksAfter;          ///< zero-suppression ticks after
      int fClusterHits;         ///< hit clustering algorithm number
      float fHitClusterDx;      ///< range in cm to look for hits to cluster in x
      float fHitClusterDyDz;    ///< range in cm to look for hits to cluster in y and z
      float fMinRMS;            ///< minimum RMS to report in case only one tick has signal on it
      float fHitMaxLen;         ///< maximum length of a hit in X, in cm
      float fHitFracADCNewHit;  ///< threshold below which if the ADC falls below multiplied by adcmax, to start a new hit
      float fHitFracADCRise;    ///< fraction ADC must rise back from the minimum to start a new hit

      std::string fRawDigitLabel;  ///< label to find the right raw digits
      const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
      const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock
    };


    CompressedHitFinder::CompressedHitFinder(fhicl::ParameterSet const & p)
    // :
    {
      fADCThreshold     = p.get<int>("ADCThreshold",5);
      fTicksBefore      = p.get<int>("TicksBefore",5);
      fTicksAfter       = p.get<int>("TicksAfter",5);
      fMinRMS           = p.get<float>("MinRMS",3);
      fRawDigitLabel    = p.get<std::string>("RawDigitLabel","daq");
      fClusterHits      = p.get<int>("ClusterHits",1);
      fHitClusterDx     = p.get<float>("HitClusterDx",1.0);
      fHitClusterDyDz   = p.get<float>("HitClusterDyDz",5.0);
      fHitMaxLen        = p.get<float>("HitMaxLen",1.0);
      fHitFracADCNewHit = p.get<float>("HitFracADCNewHit",0.5);
      fHitFracADCRise   = p.get<float>("HitFracADCRise",1.3);

      fTime    = gar::providerFrom<detinfo::DetectorClocksService>();
      fGeo     = gar::providerFrom<geo::Geometry>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();
      produces< std::vector<rec::Hit> >();
    }

    void CompressedHitFinder::produce(art::Event & e)
    {
      // create an emtpy output hit collection -- to add to.
      std::unique_ptr<std::vector<Hit> > hitCol (new std::vector<Hit> );
      std::unique_ptr<std::vector<Hit> > hitClusterCol (new std::vector<Hit> );

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
	      //LOG_WARNING("CompressedHitFinder") << " adc vector size < 2, skipping channel";
	      continue;
	    }

	  // walk through the zero-suppressed raw digits and call each block a hit

	  int nblocks = adc[1];
	  int zerosuppressedindex = nblocks*2 + 2;
	  float pos[3] = {0,0,0};
	  fGeo->ChannelToPosition(channel, pos);
	  float chanposx = pos[0];

	  for (int i=0; i<nblocks; ++i)
	    {
	      double hitSig = 0;
	      double hitTime = 0;
	      double hitSumSq = 0;
	      double hitRMS = 0;
	      unsigned int begT = adc[2+i];
	      int blocksize = adc[2+nblocks+i];
	      if (blocksize<1)
		{
		  throw cet::exception("CompressedHitFinder") << "Negative or zero block size in compressed data.";
		}
	      unsigned int endT = begT + blocksize;  // will update this as need be.

	      // divide the hits up into smaller hits 

	      int adcmax = 0;
	      int jhl = 0;

	      float distonetick = fDetProp->DriftVelocity() * (fTime->TPCTick2Time(1) - fTime->TPCTick2Time(0)) ;

	      for(int j = 0; j < blocksize; ++j)  // loop over time samples in each block
		{
		  ULong64_t t = adc[2+i]+j;
		  int a = adc[zerosuppressedindex];
		  zerosuppressedindex++;

		  ++jhl;
		  endT = begT + jhl;

		  hitSig   += a;
		  hitTime  += a*t;
		  hitSumSq += a*t*t;
		  //std::cout << "  In hit calc: " << t << " " << a << " " << hitSig << " " << hitTime << " " << hitSumSq << std::endl;

		  // make a new hit if the ADC value drops below a fraction of the max value or if we have exceeded the length limit
		  // add a check to make sure the ADC value goes back up after the dip, otherwise keep adding to the same hit.

		  bool splithit = jhl*distonetick > fHitMaxLen;
		  if (! splithit)  // only do this calc if we need to
		    {
                      if (a < adcmax*fHitFracADCNewHit)
			{
			  int zsi2 = zerosuppressedindex;
			  for (int k=j+1; k<blocksize; ++k)
			    {
			      int a2 = adc[zsi2];
			      zsi2++;
			      if (a2 > a*fHitFracADCRise)
				{
				  splithit = true;
				  break;
				}
			    }
			} 
		    }

		  if ( splithit )
		    {
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
		      if (hitRMS == 0) hitRMS = fMinRMS;
		      //std::cout << " hit RMS calc: " << hitSumSq << " " << hitSig << " " << hitTime << " " << hitRMS << std::endl;

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
		      if (hitSig>0)
			{
			  hitCol->emplace_back(channel,
					       hitSig,
					       pos,
					       begT,
					       endT,
					       hitTime,
					       hitRMS);
			}
		      jhl = 0;  // reset accumulators so we can make the next hit
	              hitSig = 0;
	              hitTime = 0;
	              hitSumSq = 0;
		      begT = endT + 1;
		      endT = begT;
		      adcmax = 0;
		    }
		  adcmax = TMath::Max(adcmax,a);
		}
	    }
        }

      // cluster hits if requested

      if (fClusterHits == 0)
        {
          e.put(std::move(hitCol));
          return;
        }
      else if (fClusterHits == 1)  // first hit-clustering algorithm. Start with hits in HitCol and make hitClusterCol
        {
          size_t nhits = hitCol->size();
          std::vector<float> hptmp(nhits);
          std::vector<int> hsi(nhits);
          std::vector<int> used(nhits,0);

          // sort hits by position in x

          for (size_t ihit=0;ihit<nhits;++ihit)
	    {
	      hptmp[ihit] = hitCol->at(ihit).Position()[0];
	    }
          TMath::Sort( (int) nhits, hptmp.data(), hsi.data() );

          // loop through hits in x, clustering as we go, marking hits as used.

          for (size_t ihitx=0; ihitx< nhits; ++ihitx)
	    {
	      size_t ihit = hsi[ihitx];  // unwound index to hit array
	      if (used[ihit]) continue;

	      const float *xyz = hitCol->at(ihit).Position();

	      // start a cluster with just this hit
	      double cpos[3] = {xyz[0], xyz[1], xyz[2]};
	      double csig = hitCol->at(ihit).Signal();
	      double cstime = hitCol->at(ihit).StartTime();
	      double cetime = hitCol->at(ihit).EndTime();
	      double crms = hitCol->at(ihit).RMS();
	      double ctime = hitCol->at(ihit).Time();

	      double xyzlow[3] =
		{
		  xyz[0] - fHitClusterDx,
		  xyz[1] - fHitClusterDyDz,
		  xyz[2] - fHitClusterDyDz
		};
	      if (xyzlow[0]<0 && xyz[0] >= 0) xyzlow[0] = 0;  // don't cluster across the cathode

	      double xyzhigh[3] =
		{
		  xyz[0] + fHitClusterDx,
		  xyz[1] + fHitClusterDyDz,
		  xyz[2] + fHitClusterDyDz
		};
	      if (xyzhigh[0]>0 && xyz[0] <= 0) xyzhigh[0] = 0;  // don't cluster across the cathode

	      for (size_t ix = ihitx; ix<nhits; ++ix) // look for candidate hits to cluster in with this one
		{
		  size_t ihc = hsi[ix];  // candidate hit to add to the cluster if it's in range
		  const float *xyz2 = hitCol->at(ihc).Position();
		  if (xyz2[0] > xyzhigh[0] || xyz2[0] < xyzlow[0]) break;
		  if (xyz2[1] < xyzhigh[1] && xyz2[1] > xyzlow[1] && xyz2[2] < xyzhigh[2] && xyz2[2] > xyzlow[2] && (used[ihc] == 0))
		    {
		      // add hit to cluster
		      used[ihc] = 1;
		      double signal = hitCol->at(ihit).Signal();
		      double totsig = csig + signal;
		      if (totsig > 0)
			{
			  for (size_t idim=0; idim<3; ++idim)
			    {
			      cpos[idim] = ( cpos[idim]*csig + xyz2[idim]*signal) / totsig;
			    }
			  cstime = TMath::Min(cstime, (double) hitCol->at(ihit).StartTime());
			  cetime = TMath::Max(cetime, (double) hitCol->at(ihit).EndTime());

			  double htime = hitCol->at(ihit).Time();
			  double hrms = hitCol->at(ihit).RMS();
			  crms = TMath::Sqrt( (csig*crms*crms + signal*hrms*hrms)/totsig  +
					      (csig*signal)*TMath::Sq(htime-ctime)/TMath::Sq(totsig));
			  ctime = (ctime*csig + htime*signal) / totsig;
			  csig = totsig;
			}
		    }
		}

	      float fcpos[3] = {0,0,0};
	      for (int i=0;i<3;++i)
		{
		  fcpos[i] = cpos[i];
		}
	      hitClusterCol->emplace_back(0,
					  csig,
					  fcpos,
					  cstime,
					  cetime,
					  ctime,
					  crms);

            }

	  e.put(std::move(hitClusterCol));
	}
    }

    DEFINE_ART_MODULE(CompressedHitFinder)

  } // namespace rec
} // namespace gar
