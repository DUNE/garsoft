/// \file    RecoBaseDrawer.cxx
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  brebel@fnal.gov
/// \version $Id: RecoBaseDrawer.cxx,v 1.3 2010/11/11 22:47:19 p-novaart Exp $
#include <cmath>
#include <map>
#include <stdint.h>

#include "TMath.h"
#include "TMarker.h"
#include "TBox.h"
#include "TH1.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TVector3.h"
#include "TText.h"
#include "TColor.h"

#include "nutools/EventDisplayBase/View2D.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "EventDisplay/eventdisplay.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "EventDisplay/Style.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Shower.h"
#include "Geometry/Geometry.h"
#include "DetectorInfo/DetectorPropertiesService.h"
#include "Utilities/AssociationUtil.h"

#include "cetlib_except/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/FindMany.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

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
namespace evd{

  //......................................................................
  RecoBaseDrawer::RecoBaseDrawer()
  {
  }
  
  //......................................................................
  RecoBaseDrawer::~RecoBaseDrawer()
  {
    
  }
  
  //......................................................................
  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void RecoBaseDrawer::Hit3D(const art::Event& evt,
                             evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;

    if(recoOpt->fDrawHits     == 0 ||
       rawOpt->fDrawRawOrReco <  1 ) return;

    int h = 0;
    for(auto const& which : recoOpt->fHitLabels) {
      
      std::vector<const rec::Hit*> hits;
      this->GetHits(evt, which, hits);
      
      this->DrawHit3D(hits, view, h%evd::kNCOLS);
      ++h;
    } // loop on imod folders
    
    return;
  }
  
  //......................................................................
  void RecoBaseDrawer::DrawHit3D(std::vector<const rec::Hit*> const& hits,
                                 evdb::View3D                      * view,
                                 int                                 color,
                                 int                                 /*marker*/,
                                 int                                 /*size*/)
  {
    //auto const* detp = gar::providerFrom<detinfo::DetectorPropertiesService>();

    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    
    // Make and fill a polymarker.
    TPolyMarker3D& pm = view->AddPolyMarker3D(hits.size(), color, 1, 3);
    
    // Display all hits on the 3D view
    size_t p = 0;
    for(auto itr : hits){
      
      // Try to get the "best" charge measurement, ie. the one last in
      // the calibration chain
      //double dQdX = itr->Signal() / detp->ElectronsToADC();
      
      // Try to get the "best" charge measurement, ie. the one last in
      // the calibration chain
      auto const* pos = itr->Position();
      
      pm.SetPoint(p, xcent+pos[0], ycent+pos[1], zcent+pos[2]);
      ++p;
    }

    return;
  }
  
  //----------------------------------------------------------------------------
  void RecoBaseDrawer::Track3D(const art::Event& evt,
                               evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;
    
    if(rawOpt->fDrawRawOrReco <  1 ) return;
    
    
    for(auto const& which : recoOpt->fTrackLabels){
      art::View<rec::Track> trackView;
      this->GetTracks(evt, which, trackView);

      if(!trackView.isValid()) continue;
      
      art::FindMany<rec::Hit> fmh(trackView, evt, which);
      
      if(!fmh.isValid()) continue;
      
      for(size_t t = 0; t < trackView.size(); ++t){
        
        int color  = evd::kColor[t%evd::kNCOLS];
        int marker = kFullDotMedium;
        int size   = 2;
        
        // Draw track using only embedded information.
        auto const& hits = fmh.at(t);
        
        DrawHit3D(hits, view, color, marker, size);
      }
    }

    return;
  }
  

  //......................................................................
  int RecoBaseDrawer::GetHits(art::Event                   const& evt,
                              std::string                  const& which,
                              std::vector<const rec::Hit*>      & hits)
  {
    hits.clear();

    std::vector<const rec::Hit*> temp;

    try{
      evt.getView(which, temp);
      for(size_t t = 0; t < temp.size(); ++t){
        hits.push_back(temp[t]);
      }
    }
    catch(cet::exception& e){
      writeErrMsg("GetHits", e);
    }
    
    return hits.size();
  }
  
  //......................................................................
  int RecoBaseDrawer::GetTracks(art::Event            const& evt,
                                std::string           const& which,
                                art::View<rec::Track>      & track)
  {
    try{
      evt.getView(which,track);
    }
    catch(cet::exception& e){
      writeErrMsg("GetTracks", e);
    }
    
    return track.vals().size();
  }

  //......................................................................
  int RecoBaseDrawer::GetShowers(art::Event             const& evt,
                                 std::string            const& which,
                                 art::View<rec::Shower>      & shower)
  {
    try{
      evt.getView(which,shower);
    }
    catch(cet::exception& e){
      writeErrMsg("GetShowers", e);
    }
    
    return shower.vals().size();
  }

}
}// namespace
////////////////////////////////////////////////////////////////////////
