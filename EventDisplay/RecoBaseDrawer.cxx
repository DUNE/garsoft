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
#include "ReconstructionDataProducts/Vertex.h"
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

  void RecoBaseDrawer::DrawVertex3D(const float *pos, 
				    evdb::View3D*       view,
				    int                 color,
				    int                 marker,
				    int                 size )
  {
    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;
    TPolyMarker3D& pm = view->AddPolyMarker3D(1, color, marker, size);
    pm.SetPoint(0,xcent + pos[0], ycent + pos[1], zcent + pos[2] );
  } 
  
  void RecoBaseDrawer::DrawHelix3D(const float *trackpar,
				   const float xpar,
				   const float xother,
                                   evdb::View3D                      * view,
				   int                                 color)
  {
    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;

    float r = 0;
    if (trackpar[2] != 0) r=1.0/trackpar[2];
    float si = 0;
    if (trackpar[4] != 0) si=1.0/trackpar[4];

    float ycc = trackpar[0] + r*TMath::Cos(trackpar[3]);
    float zcc = trackpar[1] - r*TMath::Sin(trackpar[3]);

    float dphimax = TMath::Pi();
    float phi1 = 0;
    float phi2 = (xother-xpar)*trackpar[2]*trackpar[4];
    if (phi2-phi1>dphimax) phi2 = phi1+dphimax;
    if (phi2-phi1<-dphimax) phi2 = phi1-dphimax;
    if (phi2-phi1==0) phi2=phi1+0.01;

    int nptshelix=100;
    float dphihelix=(phi2-phi1)/((float) nptshelix);

    TPolyLine3D& tpoly = view->AddPolyLine3D(nptshelix,color,1,1);

    for (int ipoint=0;ipoint<nptshelix;++ipoint)
      {
	float philoc = phi1 + ipoint*dphihelix;
	float xl = xpar + r*si*(philoc);
	float yl = ycc - r*TMath::Cos(philoc + trackpar[3]);
	float zl = zcc + r*TMath::Sin(philoc + trackpar[3]);
	tpoly.SetPoint(ipoint,xcent+xl,ycent+yl,zcent+zl);
      }
  }


  void RecoBaseDrawer::DrawArrow3D(const float *startpos,
	                     	   const float *arrowvec,
                                   evdb::View3D*       view,
                                   int                 color,
		                   float lengthscale)
  {
    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;
    TVector3 cent(xcent,ycent,zcent);

    TPolyLine3D& tpoly = view->AddPolyLine3D(5,color,1,1);
    TVector3 spv(startpos);
    TVector3 spcv = startpos + cent;
    TVector3 av(arrowvec);
    TVector3 endpos = spcv + lengthscale*av;

    tpoly.SetPoint(0,spcv.X(),spcv.Y(),spcv.Z());
    tpoly.SetPoint(1,endpos.X(),endpos.Y(),endpos.Z());

    TVector3 xhat(1,0,0);
    TVector3 yhat(0,1,0);
    TVector3 zhat(0,0,1);

    TVector3 perp1 = av.Cross(xhat);
    if (perp1.Mag() == 0)
      {
	perp1 = av.Cross(yhat);
      }
    if (perp1.Mag() == 0)
      {
	perp1 = av.Cross(zhat);
      }
    if (perp1.Mag() == 0)
      {
         tpoly.SetPoint(2,endpos.X(),endpos.Y(),endpos.Z());
         tpoly.SetPoint(3,endpos.X(),endpos.Y(),endpos.Z());
         tpoly.SetPoint(4,endpos.X(),endpos.Y(),endpos.Z());
	 return;
      }
    perp1 = perp1*(1.0/perp1.Mag());

    TVector3 arpd1 = endpos -lengthscale*0.1*av + perp1*lengthscale*0.1*av.Mag();
    TVector3 arpd2 = endpos -lengthscale*0.1*av - perp1*lengthscale*0.1*av.Mag();
    tpoly.SetPoint(2,arpd1.X(),arpd1.Y(),arpd1.Z());
    tpoly.SetPoint(3,endpos.X(),endpos.Y(),endpos.Z());
    tpoly.SetPoint(4,arpd2.X(),arpd2.Y(),arpd2.Z());
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
    xcent = 0;
    ycent = 0;
    zcent = 0;

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
        int marker = 20;
        int size   = 2;
        
        // Draw track using only embedded information.
        auto const& hits = fmh.at(t);
        
        DrawHit3D(hits, view, color, marker, size);

      }

      size_t icounter=0;
      for (auto tv = trackView.begin(); tv != trackView.end(); ++tv)
	{
          int color  = evd::kColor[icounter%evd::kNCOLS];
	  DrawHelix3D((*tv)->TrackParBeg(), (*tv)->Vertex()[0], (*tv)->End()[0], view, color);
	  DrawHelix3D((*tv)->TrackParEnd(), (*tv)->End()[0], (*tv)->Vertex()[0], view, color);

	  DrawArrow3D((*tv)->Vertex(),(*tv)->VtxDir(),view,0);
	  DrawArrow3D((*tv)->End(),(*tv)->EndDir(),view,0);

	  ++icounter;
	}
    }

    return;
  }


  //----------------------------------------------------------------------------
  void RecoBaseDrawer::Vertex3D(const art::Event& evt,
                                evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;
    
    if(rawOpt->fDrawRawOrReco <  1 ) return;
    
    for(auto const& which : recoOpt->fVertexLabels){
      art::View<rec::Vertex> vertexView;
      GetVertices(evt, which, vertexView);

      if(!vertexView.isValid()) continue;
            
      size_t icounter=0;
      for (auto tv = vertexView.begin(); tv != vertexView.end(); ++tv)
	{
          int color  = evd::kColor[icounter%evd::kNCOLS];
	  color = 5;  // hardwire yellow for now
	  DrawVertex3D( (*tv)->Position(), view, color, 20, 1 );
	  ++icounter;
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
  int RecoBaseDrawer::GetVertices(art::Event            const& evt,
                                  std::string           const& which,
                                  art::View<rec::Vertex>      & vertex)
  {
    try{
      evt.getView(which,vertex);
    }
    catch(cet::exception& e){
      writeErrMsg("GetVertices", e);
    }
    
    return vertex.vals().size();
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
