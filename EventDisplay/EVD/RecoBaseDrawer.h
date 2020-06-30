/// \file    RecoBaseDrawer.h
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  messier@indiana.edu
/// \version $Id: RecoBaseDrawer.h,v 1.3 2010/11/11 22:47:20 p-novaart Exp $
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOBASEDRAWER_H

#include <vector>

#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Principal/View.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Cluster.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"

#ifdef __ROOTCLING__

//namespace art { 
//    class Event;
//    class ServiceHandle;
//}

#else
#include "art/Framework/Services/Registry/ServiceHandle.h"
#endif

class TVector3;
class TH1F;
namespace evdb { 
    class View2D;
    class View3D;
}

namespace gar {

namespace rec {
    class Hit;
    class TPCCluster;
    class Track;
    class Shower;
    class VecHit;
    class Vertex;
    class Cluster;
}


namespace detinfo {
    class GArPropertiesService;
    class DetectorProperties;
}

namespace evd {

class ColorDrawingOptions;
class RawDrawingOptions;
class RecoDrawingOptions;

/// Aid in the rendering of RecoBase objects
class RecoBaseDrawer
{
public:
    RecoBaseDrawer();
    ~RecoBaseDrawer();

public:

  void Hit3D(art::Event  const& evt,
             evdb::View3D*      view);

  void TPCCluster3D(art::Event  const& evt,
             evdb::View3D*      view);

  void CaloCluster3D(art::Event  const& evt,
                     evdb::View3D*      view);

  void CaloHit3D(art::Event  const& evt,
                     evdb::View3D*      view);

  void Track3D(art::Event const& evt,
               evdb::View3D*     view);

  void VecHit3D(art::Event const& evt,
		evdb::View3D*     view);
  
  void Vertex3D(art::Event const& evt,
                evdb::View3D*     view);

  void DrawHit3D(std::vector<const rec::Hit*> const& hit,
                 evdb::View3D*                       view,
                 int                                 color,
                 int                                 marker = 2,
                 int                                 size = 2);

  void DrawTPCCluster3D(std::vector<const gar::rec::TPCCluster*> const& TPCClusters,
                                    evdb::View3D                      * view,
                                    int                                 color,
                                    int                                 marker = 2,
                                    int                                 size = 2);

  void DrawTrackPolyLine3D(std::vector<const gar::rec::TrackTrajectory*> const& trajectories,
                                    evdb::View3D                      * view,
                                    int                                 color,
                                    int                                 width = 2);

  void DrawCaloCluster3D(std::vector<const gar::rec::Cluster*> const& Clusters,
                         evdb::View3D                      * view,
			 int                                 color);

  void DrawCaloHit3D(std::vector<const gar::rec::CaloHit*> const& CaloHits,
                         evdb::View3D                      * view);

  void DrawTrack3D(rec::Track   const& track,
                   evdb::View3D*       view,
                   int                 color,
                   int                 marker = 1,
                   int                 size = 2);

  void DrawHelix3D(const float *trackpar,
		   const float xpar,
		   const float xother,
                   evdb::View3D* view,
                   int color,
		   int width);

  void DrawVecHit3D(std::vector<const rec::VecHit*> const& vechits,
		    evdb::View3D* view,
		    int color = 6,
		    int marker = 28,
		    int size = 2);

  void DrawVertex3D(const float *pos, 
		    evdb::View3D*       view,
		    int                 color = 5,
		    int                 marker = 20,
		    int                 size = 1);


  void DrawArrow3D(const float *startpos,
		   const float *arrowvec,
                   evdb::View3D*       view,
                   int                 color,
		   float lengthscale = 1.0);
  
  private:

  int GetHits(art::Event                   const& evt,
              std::string                  const& which,
              std::vector<const rec::Hit*>      & hits);

  int GetTPCClusters(art::Event                   const& evt,
              std::string                  const& which,
              std::vector<const rec::TPCCluster*>      & TPCClusters);
  
  int GetTracks(art::Event            const& evt,
                std::string           const& which,
                art::View<rec::Track>      & track);

  int GetVecHits(art::Event             const& evt,
		 std::string            const& which,
		 std::vector<const rec::VecHit*> & vechits);

  int GetVertices(art::Event             const& evt,
                  std::string            const& which,
                  art::View<rec::Vertex>      & vertex);
  
  int GetShowers(art::Event             const& evt,
                 std::string            const& which,
                 art::View<rec::Shower>      & shower);

  int GetCaloClusters(art::Event                   const& evt,
              std::string                  const& which,
              std::vector<const rec::Cluster*>      & Clusters);

  int GetCaloHits(art::Event                   const& evt,
              std::string                  const& which,
              std::vector<const rec::CaloHit*>      & CaloHits);
  
  private:

    
  };
}
}
#endif
////////////////////////////////////////////////////////////////////////
