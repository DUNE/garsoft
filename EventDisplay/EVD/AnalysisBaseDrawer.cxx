/// \file    AnalysisBaseDrawer.cxx
/// \brief   Class to aid in the rendering of AnalysisBase objects
/// \author  msoderbe@syr.edu
#include "TMarker.h"
#include "TBox.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TVector3.h"
#include "TLatex.h"

#include "nuevdb/EventDisplayBase/View2D.h"
#include "EventDisplay/EVD/eventdisplay.h"
#include "EventDisplay/EVD/Style.h"
#include "EventDisplay/EVD/AnalysisBaseDrawer.h"
#include "EventDisplay/EVD/RecoDrawingOptions.h"
#include "EventDisplay/EVD/AnalysisDrawingOptions.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Hit.h"
#include "Utilities/AssociationUtil.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <math.h>

/* unused function
namespace {
   // Utility function to make uniform error messages.
   void writeErrMsg(const char* fcn,
                    cet::exception const& e)
   {
      mf::LogWarning("AnalysisBaseDrawer") << "AnalysisBaseDrawer::" << fcn
                                           << " failed with message:\n"
                                           << e;
   }
}
*/

namespace gar {
namespace evd{

   //......................................................................
   AnalysisBaseDrawer::AnalysisBaseDrawer()
   {

   }

   //......................................................................
   AnalysisBaseDrawer::~AnalysisBaseDrawer()
   {

   }

   //......................................................................
   void AnalysisBaseDrawer::DrawDeDx(const art::Event& evt,
                                     evdb::View2D* ) // view)
   {
      art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
      art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;

      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {

         //Get Track collection
         std::string which = recoOpt->fTrackLabels[imod];
	 auto trackListHandle = evt.getValidHandle<std::vector<rec::Track> >(which);
         std::vector<art::Ptr<rec::Track> > tracklist;
         art::fill_ptr_vector(tracklist, trackListHandle);

      }
   }

  //......................................................................
   void AnalysisBaseDrawer::DrawKineticEnergy(const art::Event& evt,
                                              evdb::View2D* view)
   {
     art::ServiceHandle<evd::RecoDrawingOptions>     recoOpt;
     art::ServiceHandle<evd::AnalysisDrawingOptions> anaOpt;

     //add some legend-like labels with appropriate grayscale
     TLatex& proton_tex = view->AddLatex(2.0, 180.0, "proton");
     TLatex& kaon_tex   = view->AddLatex(2.0, 165.0, "kaon"  );
     TLatex& pion_tex   = view->AddLatex(2.0, 150.0, "pion"  );
     TLatex& muon_tex   = view->AddLatex(2.0, 135.0, "muon"  );
     proton_tex.SetTextColor(kBlack);
     kaon_tex  .SetTextColor(kGray+2);
     pion_tex  .SetTextColor(kGray+1);
     muon_tex  .SetTextColor(kGray);
     proton_tex.SetTextSize (0.075);
     kaon_tex  .SetTextSize (0.075);
     pion_tex  .SetTextSize (0.075);
     muon_tex  .SetTextSize (0.075);

       //int color = 1;

     //now get the actual data
     for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
       //Get Track collection
       std::string which = recoOpt->fTrackLabels[imod];
       auto trackListHandle = evt.getValidHandle<std::vector<rec::Track> >(which);
       std::vector<art::Ptr<rec::Track> > tracklist;
       art::fill_ptr_vector(tracklist, trackListHandle);

       //Loop over Tracks
       //for(size_t trkIter = 0; trkIter < tracklist.size(); ++trkIter){
       //  color = trkIter%evd::kNCOLS;
       //}
     }

   }
}
}// namespace
////////////////////////////////////////////////////////////////////////
