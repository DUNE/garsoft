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
#include "canvas/Persistency/Common/Assns.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art_root_io/TFileService.h"

#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"
#include "ReconstructionDataProducts/Hit.h"
#include "Geometry/GeometryGAr.h"
#include "DetectorInfo/DetectorClocksServiceGAr.h"
#include "DetectorInfo/DetectorPropertiesService.h"

#include <memory>

#include "TMath.h"
#include "TH1D.h"

namespace gar {
  namespace rec {

    class CompressedHitFinderToad : public art::EDProducer {
    public:
      explicit CompressedHitFinderToad(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      CompressedHitFinderToad(CompressedHitFinderToad const &) = delete;
      CompressedHitFinderToad(CompressedHitFinderToad &&) = delete;
      CompressedHitFinderToad & operator = (CompressedHitFinderToad const &) = delete;
      CompressedHitFinderToad & operator = (CompressedHitFinderToad &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      // Declare member data here.

      int fADCThreshold;        ///< zero-suppression threshold (in case the raw digits need to be zero-suppressed)
      int fTicksBefore;         ///< zero-suppression ticks before.  Only used if raw waveforms are not zero-suppressed and we must do it
      int fTicksAfter;          ///< zero-suppression ticks after.  Only used if raw waveforms are not zero-suppressed and we must do it

      float fMinRMS;            ///< minimum RMS to report in case only one tick has signal on it
      float fHitMaxLen;         ///< maximum length of a hit in X, in cm
      float fHitFracADCNewHit;  ///< threshold below which if the ADC falls below multiplied by adcmax, to start a new hit
      float fHitFracADCRise;    ///< fraction ADC must rise back from the minimum to start a new hit

      std::string fRawDigitLabel;  ///< label to find the right raw digits
      const detinfo::DetectorProperties*  fDetProp;      ///< detector properties
      const geo::GeometryCore*            fGeo;          ///< pointer to the geometry
      const gar::detinfo::DetectorClocks* fTime;         ///< electronics clock

      TH1D *fAvgPulseLongHist;
      TH1D *fAvgPulseShortHist;
    };


    CompressedHitFinderToad::CompressedHitFinderToad(fhicl::ParameterSet const & p)
    : EDProducer{p} 
    {
      fADCThreshold     = p.get<int>("ADCThreshold",5);
      fTicksBefore      = p.get<int>("TicksBefore",5);
      fTicksAfter       = p.get<int>("TicksAfter",5);
      fMinRMS           = p.get<float>("MinRMS",3);
      fRawDigitLabel    = p.get<std::string>("RawDigitLabel","daq");
      fHitMaxLen        = p.get<float>("HitMaxLen",1.0);
      fHitFracADCNewHit = p.get<float>("HitFracADCNewHit",0.5);
      fHitFracADCRise   = p.get<float>("HitFracADCRise",1.3);

      fTime    = gar::providerFrom<detinfo::DetectorClocksServiceGAr>();
      fGeo     = gar::providerFrom<geo::GeometryGAr>();
      fDetProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

      consumes< std::vector<raw::RawDigit> >(fRawDigitLabel);
      produces< std::vector<rec::Hit> >();
      produces< art::Assns<rec::Hit, raw::RawDigit> >();

      art::ServiceHandle<art::TFileService> tfs;
      fAvgPulseLongHist = tfs->make<TH1D>("garAvgPulseLongHist","Pulse ADC Average;tick;ADC Sum",1000,-0.5,999.5);
      fAvgPulseShortHist = tfs->make<TH1D>("garAvgPulseShortHist","Pulse ADC Average;tick;ADC Sum",200,-0.5,199.5);
    }

    void CompressedHitFinderToad::produce(art::Event & e)
    {
      // create an emtpy output hit collection -- to add to.
      std::unique_ptr<std::vector<Hit> > hitCol (new std::vector<Hit> );

      // create art::Assns between (multiple) rec::Hit and a raw::RawDigit
      std::unique_ptr< art::Assns<rec::Hit, raw::RawDigit> >
        hitDigAssns(new::art::Assns<rec::Hit, raw::RawDigit>);

      auto const hitPtrMaker = art::PtrMaker<rec::Hit>(e);

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
              MF_LOG_WARNING("CompressedHitFinder") << " Ununderstood compression mode: " << rd.Compression() << " Not making hits.";
              e.put(std::move(hitCol));
              e.put(std::move(hitDigAssns));
              return;
            }
          // use the format of the compressed raw digits -- a table of contents of the number of blocks, then all the block sizes, and then all the
          // block start locations
          if (adc.size() < 2)
            {
              //MF_LOG_WARNING("CompressedHitFinder") << " adc vector size < 2, skipping channel";
              continue;
            }

          // walk through the zero-suppressed raw digits and call each block a hit

          int nBlocks = adc[1];
          int zerosuppressedindex = nBlocks*2 + 2;
          float pos[3] = {0,0,0};
          fGeo->ChannelToPosition(channel, pos);
          float chanposx = pos[0];

          int t0adcblockhist = 0;

          for (int iBlock=0; iBlock<nBlocks; ++iBlock)
            {
              double hitSig = 0;
              double hitTime = 0;
              double hitSumSq = 0;
              double hitRMS = 0;
              unsigned int begT = adc[2+iBlock];
              int blocksize = adc[2+nBlocks+iBlock];
              if (blocksize<1)
                {
                  throw cet::exception("CompressedHitFinder") << "Negative or zero block size in compressed data.";
                }
              unsigned int endT = begT + blocksize;  // will update this as need be.

              // divide the hits up into smaller hits 

              int adcmax = 0;
              int jhl = 0;

              double distonetick = fDetProp->DriftVelocity() * (fTime->TPCTick2Time(1) - fTime->TPCTick2Time(0)) ;

              for(int jInBlock = 0; jInBlock < blocksize; ++jInBlock)  // loop over time samples in each block
                {
                  ULong64_t t = adc[2+iBlock]+jInBlock;
                  int a = adc[zerosuppressedindex];
                  zerosuppressedindex++;

                  fAvgPulseShortHist->Fill(jInBlock,a);
                  if (t0adcblockhist == 0)
                    {
                      t0adcblockhist = t;
                    }
                  int deltat = t - t0adcblockhist;
                  if (deltat < fAvgPulseLongHist->GetNbinsX())
                    {
                      fAvgPulseLongHist->Fill(deltat,a);
                    }
                  else
                    {
                      t0adcblockhist = 0;
                    }
                  endT = begT + jhl;
                  ++jhl;

                  hitSig   += a;
                  hitTime  += a*t;
                  hitSumSq += a*t*t;
                  //std::cout << "  In hit calc: " << t << " " << a << " " << hitSig << " " << hitTime << " " << hitSumSq << std::endl;

                  // make a new hit if the ADC value drops below a fraction of the max value or if we have exceeded the length limit
                  // add a check to make sure the ADC value goes back up after the dip, otherwise keep adding to the same hit.

                  bool splithit = jhl*distonetick > fHitMaxLen;
                  if (!splithit)  // only do this calc if we need to
                    {
                      if (a < adcmax*fHitFracADCNewHit)
                        {
                          int zsi2 = zerosuppressedindex;
                          for (int kInBlock=jInBlock+1; kInBlock<blocksize; ++kInBlock)
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
                  adcmax = TMath::Max(adcmax,a);

		  // create a hit if we have split a hit or if we are out of waveforms

                  if ( splithit || jInBlock == (blocksize - 1))
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

                      double driftdistance = fDetProp->DriftVelocity() * fTime->TPCTick2Time(hitTime);
                      if (chanposx < fGeo->TPCXCent())
                        {
                          pos[0] = chanposx + driftdistance;
                        }
                      else
                        {
                          pos[0] = chanposx + driftdistance;
                        }

		      // noise hits are plentiful
                      //if (hitSig < 0)
                      //  {
                      //    MF_LOG_WARNING("CompressedHitFinder") << "Negative Signal in hit finder" << std::endl;
                      //  }
                      if (hitSig>0)
                        {
                          rec::Hit newHit(channel,hitSig,pos,begT,endT,hitTime,hitRMS);
                          hitCol->emplace_back(newHit);

                          // Make art::Assn<Hit,RawDigit>
                          art::Ptr<rec::Hit>      hitPtr = hitPtrMaker(hitCol->size()-1);
                          art::Ptr<raw::RawDigit> digPtr = art::Ptr<raw::RawDigit>(rdCol,ird);
                          hitDigAssns->addSingle(hitPtr,digPtr);
                        }
                      jhl = 0;  // reset accumulators so we can make the next hit
                      hitSig = 0;
                      hitTime = 0;
                      hitSumSq = 0;
                      begT = endT + 1;
                      endT = begT;
                      adcmax = 0;
                    }
                }
            }
        }

      // cluster hits if requested

      e.put(std::move(hitCol));
      e.put(std::move(hitDigAssns));
      return;

    }

    DEFINE_ART_MODULE(CompressedHitFinderToad)

  } // namespace rec
} // namespace gar
