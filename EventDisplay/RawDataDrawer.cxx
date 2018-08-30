/**
* @file    RawDataDrawer.cxx
* @brief   Class to aid in the rendering of RawData objects
* @author  trj@fnal.gov
*
* This class prepares the rendering of the raw digits content for the 3D view
* of the display.
*/
#include <cmath>       // std::abs(), ...
#include <utility>     // std::pair<>, std::move()
#include <memory>      // std::unique_ptr()
#include <tuple>
#include <limits>      // std::numeric_limits<>
#include <type_traits> // std::add_const_t<>, ...
#include <typeinfo>    // to use typeid()
#include <algorithm>   // std::fill(), std::find_if(), ...
#include <cstddef>     // std::ptrdiff_t

#include "TH1F.h"
#include "TMarker3DBox.h"
#include "TFrame.h"
#include "TVirtualPad.h"

#include "EventDisplay/ChangeTrackers.h"
#include "EventDisplay/RawDataDrawer.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "RawDataProducts/RawDigit.h"
#include "RawDataProducts/raw.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/GArPropertiesService.h"
#include "DetectorInfo/ECALPropertiesService.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "DetectorInfo/DetectorClocksService.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/demangle.h"

#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TBox.h"
#include "TText.h"
#include "TMath.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
    cet::exception const& e)
    {
      mf::LogWarning("RecoBaseDrawer") << "RecoBaseDrawer::" << fcn
      << " failed with message:\n"
      << e;
    }
  }


  namespace gar {
    namespace evd {

      //......................................................................
      RawDataDrawer::RawDataDrawer()
      {
        fEventTQHist = (TH1F*) new TH1F("eventTQHist","ADC Sum vs Time, All Channels;time (ticks);ADC Sum",100,0,10000);
        fDigitTQHist = (TH1F*) new TH1F("digitTQHist","ADC Sum vs Time, One Channel;time (ticks);ADC Sum",100,0,10000);
      }

      //......................................................................
      RawDataDrawer::~RawDataDrawer()
      {
      }

      //......................................................................
      void RawDataDrawer::FillQHisto(gar::raw::RawDigit const& dig,
        TH1F*                     histo)
        {

          gar::raw::ADCvector_t uncompressed;
          uncompressed.resize(dig.Samples());
          short ped = dig.Pedestal();
          gar::raw::Uncompress(dig.ADCs(),uncompressed,ped,dig.Compression());
          for(size_t t = 0; t < uncompressed.size(); ++t)
          histo->Fill(t, 1. * uncompressed[t] - 1. * dig.Pedestal());

          return;
        }//end loop over raw hits

        // to do -- need to get xyz position of a channel from the channel map -- it's not just a rectangular grid of pixels

        //......................................................................
        void RawDataDrawer::RawDigit3D(const art::Event& evt,
          evdb::View3D*     view)
          {
            // Check if we're supposed to draw raw hits at all
            art::ServiceHandle<evd::RawDrawingOptions>   rawopt;
            art::ServiceHandle<evd::ColorDrawingOptions> cdopt;
            art::ServiceHandle<geo::Geometry>            geom;
            auto fTime = gar::providerFrom<detinfo::DetectorClocksService>();
            auto detProp = gar::providerFrom<detinfo::DetectorPropertiesService>();

            double driftVelocity = detProp->DriftVelocity(detProp->Efield(),
            detProp->Temperature());

            if(rawopt->fDrawRawOrReco > 0) return;

            std::vector<const raw::RawDigit*> rawhits;
            this->GetRawDigits(evt, rawhits);

            fEventTQHist->Reset();
            fDigitTQHist->Reset();

            //auto const& colorSet = cdopt->RawQ();

            gar::raw::ADCvector_t uncompressed;
            for(auto dig : rawhits){

              this->FillQHisto(*dig, fEventTQHist);
              if(dig->Channel() == rawopt->fChannel)
              this->FillQHisto(*dig, fDigitTQHist);

              // draw the digit as time, y, z with the latter coordinates
              // given by the channel map.
              float xyz[3];

              geom->ChannelToPosition((unsigned int)dig->Channel(),
              xyz);

              //std::cout << "RDGeo: " << dig->Channel() << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
              double chanposx = xyz[0];

              // draw hit channels on the readout planes in gray.  No timing or x information.
              //int c = kGray;
              //int s=1;
              //int w=1;
              //TPolyLine3D& rdpos = view->AddPolyLine3D(2,c,w,s);
              //rdpos.SetPoint(0,xyz[0],xyz[1],xyz[2]);
              //rdpos.SetPoint(1,xyz[0]+1,xyz[1],xyz[2]);
              //// std::cout << "adding a raw digit: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;

              // loop over the ADC values for this digit and draw the value if it
              // is above the threshold
              //TMarker3DBox *box = nullptr;

              uncompressed.resize(dig->Samples());
              short ped = dig->Pedestal();
              gar::raw::Uncompress(dig->ADCs(),uncompressed,ped,dig->Compression());

              for(size_t t = 0; t < dig->Samples(); ++t){
                if(uncompressed[t] < rawopt->fMinSignal) continue;

                double driftdistance = fTime->TPCTick2Time(t) * driftVelocity;

                if (chanposx < 0)
                {
                  xyz[0] = chanposx + driftdistance;
                }
                else
                {
                  xyz[0] = chanposx - driftdistance;
                }

                // somehow these boxes don't work -- draw short lines instead
                //box = &(view->AddMarker3DBox(xyz[0],
                //                            xyz[1],
                //                            xyz[2],
                //                            0.5,    // the extent is 1/2 tick
                //			       0.5 * 0.3, // to fix  geom->ChannelPitch(),
                //                         0.5 * 0.3  // to fix geom->ChannelPitch()
                //			       ));
                //box->SetFillStyle(1001);
                ////auto tmpcolor = colorSet.GetColor(uncompressed[t]);
                //box->SetFillColor(kGreen);
                //box->SetBit(kCannotPick);

                //std::cout << "adding a raw digit: " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;

                // try drawing with polylines, little red segments


                int c2 = kGreen;
                int s2=1;
                int w2=1;

                TPolyLine3D& rdpos = view->AddPolyLine3D(2,c2,w2,s2);
                rdpos.SetPoint(0,xyz[0]-0.5,xyz[1],xyz[2]);
                rdpos.SetPoint(1,xyz[0]+0.5,xyz[1],xyz[2]);

              } // end loop over ADC values for the digit
            }//end loop over raw digits

            return;
          }

          //......................................................................
          void RawDataDrawer::GetRawDigits(art::Event                        const& evt,
            std::vector<const raw::RawDigit*>      & digits)
            {
              digits.clear();

              art::ServiceHandle<evd::RawDrawingOptions> rawopt;

              std::vector<const raw::RawDigit*> temp;

              try{
                evt.getView(rawopt->fRawDataLabel, temp);
                for(size_t t = 0; t < temp.size(); ++t){
                  digits.push_back(temp[t]);
                }
              }
              catch(cet::exception& e){
                writeErrMsg("GetDigits", e);
              }

              return;
            } // RawDataDrawer::GetRawDigits()

            } // namespace evd
          }

          ////////////////////////////////////////////////////////////////////////
