/// \file    RecoBaseDrawer.h
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  messier@indiana.edu
/// \version $Id: RecoBaseDrawer.h,v 1.3 2010/11/11 22:47:20 p-novaart Exp $
#ifndef EVD_RECOBASEDRAWER_H
#define EVD_RECOBASEDRAWER_H

#include <vector>

#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Principal/View.h"

#ifdef __ROOTCLING__

namespace art { 
    class Event;
    class ServiceHandle;
}

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
    class Track;
    class Shower;
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

  void Track3D(art::Event const& evt,
               evdb::View3D*     view);

  void DrawHit3D(std::vector<const rec::Hit*> const& hit,
                 evdb::View3D*                       view,
                 int                                 color,
                 int                                 marker = 2,
                 int                                 size = 2);

  void DrawTrack3D(rec::Track   const& track,
                   evdb::View3D*       view,
                   int                 color,
                   int                 marker = 1,
                   int                 size = 2);

  void DrawHelix3D(const float  *trackpar,
		   const float xpar,
		   const float xother,
                   evdb::View3D*       view,
                   int                 color);

  void DrawArrow3D(const float *startpos,
		   const float *arrowvec,
                   evdb::View3D*       view,
                   int                 color,
		   float lengthscale = 1.0);
  
  private:

  int GetHits(art::Event                   const& evt,
              std::string                  const& which,
              std::vector<const rec::Hit*>      & hits);
  
  int GetTracks(art::Event            const& evt,
                std::string           const& which,
                art::View<rec::Track>      & track);
  
  int GetShowers(art::Event             const& evt,
                 std::string            const& which,
                 art::View<rec::Shower>      & shower);

  
  private:

    
  };
}
}
#endif
////////////////////////////////////////////////////////////////////////
