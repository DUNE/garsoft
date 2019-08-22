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
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Vertex.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/Shower.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/CaloHit.h"
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
      
      std::vector<const gar::rec::Hit*> hits;
      this->GetHits(evt, which, hits);
      
      this->DrawHit3D(hits, view, h%evd::kNCOLS);
      ++h;
    } // loop on imod folders
    
    return;
  }


  //......................................................................
  ///
  /// Render TPCCluster objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void RecoBaseDrawer::TPCCluster3D(const art::Event& evt,
                             evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;

    if(recoOpt->fDrawTPCClusters     == 0 ||
       rawOpt->fDrawRawOrReco <  1 ) return;

    int h = 0;
    for(auto const& which : recoOpt->fTPCClusterLabels) {
      
      std::vector<const gar::rec::TPCCluster*> TPCClusters;
      this->GetTPCClusters(evt, which, TPCClusters);
      
      this->DrawTPCCluster3D(TPCClusters, view, h%evd::kNCOLS);
      ++h;
    } 
    
    return;
  }


  //......................................................................
  ///
  /// Render Calorimeter Cluster objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  ///
  void RecoBaseDrawer::CaloCluster3D(const art::Event& evt,
                                     evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;

    if(recoOpt->fDrawCaloClusters     == 0 ||
       rawOpt->fDrawRawOrReco <  1 ) return;

    int h = 0;
    for(auto const& which : recoOpt->fCaloClusterLabels) {
      
      std::vector<const gar::rec::Cluster*> Clusters;
      this->GetCaloClusters(evt, which, Clusters);
      
      //this->DrawCaloCluster3D(Clusters, view, h%evd::kNCOLS);  
      this->DrawCaloCluster3D(Clusters, view, 3);  // think about colors. currently just green
      ++h;
    } 
    
    return;
  }


  //......................................................................
  ///
  /// Render Calorimeter Hit objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  ///
  void RecoBaseDrawer::CaloHit3D(const art::Event& evt,
                                     evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;

    if(recoOpt->fDrawCaloHits     == 0 ||
       rawOpt->fDrawRawOrReco <  1 ) return;

    int h = 0;
    for(auto const& which : recoOpt->fCaloHitLabels) {
      
      std::vector<const gar::rec::CaloHit*> CaloHits;
      this->GetCaloHits(evt, which, CaloHits);
      
      this->DrawCaloHit3D(CaloHits, view);  
      ++h;
    } 
    
    return;
  }

  void RecoBaseDrawer::DrawVertex3D(const float *pos, 
				    evdb::View3D* view,
				    int color,
				    int marker,
				    int size )
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
                                   evdb::View3D *view,
				   int color)
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
    float si = TMath::Tan(trackpar[4]);

    float ycc = trackpar[0] + r*TMath::Cos(trackpar[3]);
    float zcc = trackpar[1] - r*TMath::Sin(trackpar[3]);

    float dphimax = TMath::Pi();
    float phi1 = 0;
    float phi2 = (xother-xpar)*trackpar[2]/si;
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
  void RecoBaseDrawer::DrawHit3D(std::vector<const gar::rec::Hit*> const& hits,
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


  //......................................................................
  void RecoBaseDrawer::DrawTPCCluster3D(std::vector<const gar::rec::TPCCluster*> const& TPCClusters,
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
    TPolyMarker3D& pm = view->AddPolyMarker3D(TPCClusters.size(), color, 1, 3);
    
    // Display all TPCClusters on the 3D view
    size_t p = 0;
    for(auto itr : TPCClusters){
      
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


  //......................................................................
  void RecoBaseDrawer::DrawCaloCluster3D(std::vector<const gar::rec::Cluster*> const& Clusters,
					 evdb::View3D                      * view,
					 int                               color)
  {
    //auto const* detp = gar::providerFrom<detinfo::DetectorPropertiesService>();

    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;

    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

    // Make and fill polylines
    
    for(auto itr : Clusters){
      TPolyLine3D& pl = view->AddPolyLine3D(14, color, 1, 1);
      
      auto const* pos = itr->Position();
      auto const* dir = itr->EigenVectors();
      //auto const* shape = itr->Shape();
      float length = itr->Energy()*recoOpt->fCaloClusterScale;  

      // vertices of an octahedron
      float poslist[6][3];
      for (int i=0;i<3;++i)
	{
	  poslist[0][i] = pos[i];
	  poslist[5][i] = pos[i]+dir[i]*length;
	  poslist[1][i] = pos[i]+(dir[i] + dir[3+i]*0.2)*length/2.0; 
	  poslist[2][i] = pos[i]+(dir[i] + dir[6+i]*0.2)*length/2.0; 
	  poslist[3][i] = pos[i]+(dir[i] - dir[3+i]*0.2)*length/2.0; 
	  poslist[4][i] = pos[i]+(dir[i] - dir[6+i]*0.2)*length/2.0; 
	}
      for (int j=0;j<6; ++j)
	{
	  poslist[j][0] += xcent;
	  poslist[j][1] += ycent;
	  poslist[j][2] += zcent;
	}
      pl.SetPoint(0, poslist[0][0], poslist[0][1], poslist[0][2]);
      pl.SetPoint(1, poslist[1][0], poslist[1][1], poslist[1][2]);
      pl.SetPoint(2, poslist[5][0], poslist[5][1], poslist[5][2]);
      pl.SetPoint(3, poslist[2][0], poslist[2][1], poslist[2][2]);
      pl.SetPoint(4, poslist[0][0], poslist[0][1], poslist[0][2]);
      pl.SetPoint(5, poslist[3][0], poslist[3][1], poslist[3][2]);
      pl.SetPoint(6, poslist[5][0], poslist[5][1], poslist[5][2]);
      pl.SetPoint(7, poslist[4][0], poslist[4][1], poslist[4][2]);
      pl.SetPoint(8, poslist[0][0], poslist[0][1], poslist[0][2]);
      pl.SetPoint(9, poslist[1][0], poslist[1][1], poslist[1][2]);
      pl.SetPoint(10, poslist[2][0], poslist[2][1], poslist[2][2]);
      pl.SetPoint(11, poslist[3][0], poslist[3][1], poslist[3][2]);
      pl.SetPoint(12, poslist[4][0], poslist[4][1], poslist[4][2]);
      pl.SetPoint(13, poslist[1][0], poslist[1][1], poslist[1][2]);
    }

    return;
  }



  //......................................................................
  void RecoBaseDrawer::DrawCaloHit3D(std::vector<const gar::rec::CaloHit*> const& CaloHits,
                                     evdb::View3D                      * view)
  {
    //auto const* detp = gar::providerFrom<detinfo::DetectorPropertiesService>();

    const int colorheatindex[16] = { kBlue+2, kBlue+1, kBlue, 
				     kCyan+2, kCyan+1, kCyan,
				     kGreen+2, kGreen+1, kGreen, 
				     kRed-2, kRed+1, kRed, 
                                     kYellow+2, kYellow+1, kYellow, 
				     kWhite };

    art::ServiceHandle<geo::Geometry> geo;
    double xcent = geo->TPCXCent();
    double ycent = geo->TPCYCent();
    double zcent = geo->TPCZCent();
    xcent = 0;
    ycent = 0;
    zcent = 0;

    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

    // Make and fill polylines -- cubes in this case
    
    for(auto itr : CaloHits){
      
      auto const* pos = itr->Position();
      float energy=itr->Energy();
      if (energy<=0) energy=1E-3;
      energy = TMath::Log(energy);
      float normenergy = (energy + recoOpt->fCaloHitScale)/recoOpt->fCaloHitScale;
      if (normenergy < 0)  normenergy = 0;
      if (normenergy > 1) normenergy = 1;
      int icolor = normenergy*16.0;
      if (icolor>15) icolor = 15;
      float halflength = 0.5;  // in cm
      TPolyLine3D& pl = view->AddPolyLine3D(14, colorheatindex[icolor], 1, 1);

      //std::cout << "Drawing calo hit at : " << pos[0] << " " << pos[1] << " " << pos[2] << " " << 
      //	"with normenergy: " << normenergy << std::endl;

      // vertices of a cube

      float poslist[8][3];
      std::set<int> bot = {0, 1, 2, 3};
      std::set<int> left = {0, 3, 4, 7};
      std::set<int> front = {0, 1, 5, 4};

      for (int i=0; i<8; ++i)
	{
	  if (bot.find(i) == bot.end())  // we're on top
	    {
	      poslist[i][2] = pos[2] + halflength;
	    }
	  else
	    {
	      poslist[i][2] = pos[2] - halflength;
	    }
	  if (left.find(i) == left.end())  // we're on the right-hand side
	    {
	      poslist[i][0] = pos[0] + halflength;
	    }
	  else
	    {
	      poslist[i][0] = pos[0] - halflength;
	    }
	  if (front.find(i) == front.end())  // we're in back
	    {
	      poslist[i][1] = pos[1] + halflength;
	    }
	  else
	    {
	      poslist[i][1] = pos[1] - halflength;
	    }
	}
      // move the origin as needed
      for (int j=0;j<8; ++j)
	{
	  poslist[j][0] += xcent;
	  poslist[j][1] += ycent;
	  poslist[j][2] += zcent;
	}

      // a little double-coverage here 

      pl.SetPoint(0, poslist[0][0], poslist[0][1], poslist[0][2]);
      pl.SetPoint(1, poslist[1][0], poslist[1][1], poslist[1][2]);
      pl.SetPoint(2, poslist[5][0], poslist[5][1], poslist[5][2]);
      pl.SetPoint(3, poslist[6][0], poslist[6][1], poslist[6][2]);
      pl.SetPoint(4, poslist[2][0], poslist[2][1], poslist[2][2]);

      pl.SetPoint(5, poslist[3][0], poslist[3][1], poslist[3][2]);
      pl.SetPoint(6, poslist[0][0], poslist[0][1], poslist[0][2]);
      pl.SetPoint(7, poslist[4][0], poslist[4][1], poslist[4][2]);
      pl.SetPoint(8, poslist[7][0], poslist[7][1], poslist[7][2]);
      pl.SetPoint(9, poslist[6][0], poslist[6][1], poslist[6][2]);

      pl.SetPoint(10, poslist[7][0], poslist[7][1], poslist[7][2]);
      pl.SetPoint(11, poslist[3][0], poslist[3][1], poslist[3][2]);
      pl.SetPoint(12, poslist[2][0], poslist[2][1], poslist[2][2]);
      pl.SetPoint(13, poslist[1][0], poslist[1][1], poslist[1][2]);
    }

    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawVecHit3D(std::vector<const gar::rec::VecHit*> const& vechits,
				    evdb::View3D *view,
				    int color,
				    int /*marker*/,
				    int /*size*/)
  {
    size_t p =0;
    for(unsigned int idx = 0; idx < vechits.size(); ++idx){

      int color  = evd::kColor[idx%evd::kNCOLS];

      // Make and fill a polyline: the vechit extent in both directions
      TPolyLine3D& pl = view->AddPolyLine3D(2, color, 1, 1);
      // Make and fill a polymarker: the vechit position
      // TODO: Could add one but how best to display?
      //TPolyMarker3D& pm = view->AddPolyMarker3D(1, 6, 33, 1);

      const float *pos = vechits[idx]->Position();
      const float *dir = vechits[idx]->Direction();

      const double len = vechits[idx]->Length();

      //pm.SetPoint(0, pos[0], pos[1], pos[2]);

      pl.SetPoint(0, pos[0] - dir[0]*(len/2),
		     pos[1] - dir[1]*(len/2),
		     pos[2] - dir[2]*(len/2));

      // This would be the center, but don't need to draw that
      //pl.SetPoint(1, pos[0], pos[1], pos[2]);
      
      pl.SetPoint(1, pos[0] + dir[0]*(len/2),
		     pos[1] + dir[1]*(len/2),
		     pos[2] + dir[2]*(len/2));
      ++p;
    }
  }

  //----------------------------------------------------------------------------
  void RecoBaseDrawer::Track3D(const art::Event& evt,
                               evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions>  rawOpt;
    
    if(rawOpt->fDrawRawOrReco <  1 ) return;
    bool drawTracks = (recoOpt->fDrawTracks != 0);    

    if(drawTracks){
      for(auto const& which : recoOpt->fTrackLabels){
	art::View<gar::rec::Track> trackView;
	this->GetTracks(evt, which, trackView);
	
	if(!trackView.isValid()) continue;
	
	art::FindMany<gar::rec::TPCCluster> fmc(trackView, evt, which);
	
	if(!fmc.isValid()) continue;
	
	for(size_t t = 0; t < trackView.size(); ++t){
	  
	  int color  = evd::kColor[t%evd::kNCOLS];
	  int marker = 20;
	  int size   = 2;
        
	  // Draw track using only embedded information.
	  auto const& TPCClusters = fmc.at(t);
	  
	  DrawTPCCluster3D(TPCClusters, view, color, marker, size);
	  
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
    bool drawVertices = (recoOpt->fDrawVertices != 0);

    if(drawVertices){
      for(auto const& which : recoOpt->fVertexLabels){
	art::View<gar::rec::Vertex> vertexView;
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
    }
    
    return;
  }  

  //----------------------------------------------------------------------------
  void RecoBaseDrawer::VecHit3D(const art::Event& evt,
                                evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    art::ServiceHandle<evd::RawDrawingOptions> rawOpt;

    if(rawOpt->fDrawRawOrReco < 1) return;
    
    bool drawVecHit = (recoOpt->fDrawVecHits != 0);

    if(drawVecHit){
      for(auto const& which : recoOpt->fVecHitLabels){
	std::vector<const gar::rec::VecHit*> vhit;
	this->GetVecHits(evt, which, vhit);
	DrawVecHit3D(vhit, view);
      } // loop on imod folders
    }
  }
  
  //......................................................................
  int RecoBaseDrawer::GetHits(art::Event                   const& evt,
                              std::string                  const& which,
                              std::vector<const gar::rec::Hit*>      & hits)
  {
    hits.clear();

    std::vector<const gar::rec::Hit*> temp;

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
  int RecoBaseDrawer::GetTPCClusters(art::Event                   const& evt,
                              std::string                  const& which,
                              std::vector<const gar::rec::TPCCluster*>      & TPCClusters)
  {
    TPCClusters.clear();

    std::vector<const gar::rec::TPCCluster*> temp;

    try{
      evt.getView(which, temp);
      for(size_t t = 0; t < temp.size(); ++t){
        TPCClusters.push_back(temp[t]);
      }
    }
    catch(cet::exception& e){
      writeErrMsg("GetTPCClusters", e);
    }
    
    return TPCClusters.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetCaloClusters(art::Event           const& evt,
                              std::string                  const& which,
                              std::vector<const gar::rec::Cluster*>      & Clusters)
  {
    Clusters.clear();

    std::vector<const gar::rec::Cluster*> temp;

    try{
      evt.getView(which, temp);
      for(size_t t = 0; t < temp.size(); ++t){
        Clusters.push_back(temp[t]);
      }
    }
    catch(cet::exception& e){
      writeErrMsg("GetCaloClusters", e);
    }
    
    return Clusters.size();
  }


  //......................................................................
  int RecoBaseDrawer::GetCaloHits(art::Event           const& evt,
                              std::string                  const& which,
                              std::vector<const gar::rec::CaloHit*>      & CaloHits)
  {
    CaloHits.clear();

    std::vector<const gar::rec::CaloHit*> temp;

    try{
      evt.getView(which, temp);
      for(size_t t = 0; t < temp.size(); ++t){
        CaloHits.push_back(temp[t]);
      }
    }
    catch(cet::exception& e){
      writeErrMsg("GetCaloHits", e);
    }
    
    return CaloHits.size();
  }
  
  //......................................................................
  int RecoBaseDrawer::GetTracks(art::Event            const& evt,
                                std::string           const& which,
                                art::View<gar::rec::Track>      & track)
  {
    try{
      evt.getView(which, track);
    }
    catch(cet::exception& e){
      writeErrMsg("GetTracks", e);
    }
    
    return track.vals().size();
  }

  //......................................................................
  int RecoBaseDrawer::GetVecHits(art::Event                 const& evt,
				 std::string                const& which,
				 std::vector<const gar::rec::VecHit*> & vechit)
  {
    std::vector<const gar::rec::VecHit*> temp(vechit);
    try{
      evt.getView(which, temp);
      temp.swap(vechit);
    }
    catch(cet::exception& e){
      writeErrMsg("GetVecHits", e);
    }
    
    return vechit.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetVertices(art::Event            const& evt,
                                  std::string           const& which,
                                  art::View<gar::rec::Vertex>     & vertex)
  {
    try{
      evt.getView(which, vertex);
    }
    catch(cet::exception& e){
      writeErrMsg("GetVertices", e);
    }
    
    return vertex.vals().size();
  }

  //......................................................................
  int RecoBaseDrawer::GetShowers(art::Event             const& evt,
                                 std::string            const& which,
                                 art::View<gar::rec::Shower>      & shower)
  {
    try{
      evt.getView(which, shower);
    }
    catch(cet::exception& e){
      writeErrMsg("GetShowers", e);
    }
    
    return shower.vals().size();
  }
  
}
}// namespace
////////////////////////////////////////////////////////////////////////
