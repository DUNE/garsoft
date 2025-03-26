#ifndef TRACKPROPAGATOR_H
#define TRACKPROPAGATOR_H

////////////////////////////////////////////////////////////////////////
/// \file  TrackPropagator.h
/// \brief Interface to propagate a Track to the specific point
///
/// \version $Id:  $
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace gar {
  namespace rec {
    class Track;
  }
}

namespace util {
  class TrackPropagator {

  public:
    //No constructor
    TrackPropagator() = delete;
    //No copy constructor
    TrackPropagator(const TrackPropagator&) = delete;
    //No destructor
    ~TrackPropagator() = delete;

    /** IMPORTANT: Make sure that you have a consistent set of coordinate
        systems in the arguments passed to all these methods.  As a rule,
        track parameters are expressed in the coordinate system implicitly
        defined by the geometry gdml file specified in the ---.fcl file.  But
        in writing algorithms that check if some track extrapolates out to some
        point or whatever, you might be thinking in another coordinate system,
        e.g. one where the center of the MPD is at (0, 0, 0).  That would be
        an error. Consider trying instead something like

        int result = util::TrackPropagator::PropagateToCylinder(trk->TrackParEnd(), trk->End(),
            fGeo->GetECALInnerBarrelRadius(),fGeo->TPCYCent(),fGeo->TPCZCent(), xyz_intersection,
            fGeo->GetECALEndcapStartX() );

        where fGeo is an

        const gar::geo::GeometryCore* fGeo = gar::providerFrom<geo::GeometryGAr>();
        */

    /** Finds two intersections of a helix with an infinite cylinder of radius
        rCyl, centered at yCyl, zCyl and parallel to the x axis, storing the
        "closest" intersection in retXYZ1, meaning that the value of phi at the
        intersection is closest to what it is at Xpoint.  The other intersection is
        returned in retXYZ2.  Status codes: 0 is success; 1,2, 3 and 4 are
        rno-intersection conditions.  If the helix cylinder and given cylinder
        are disjoint, returns 1; if one is entirely inside the other,
        status is 2; if the cylinders have a common axes, status is 3.  Status
        code 4 results when the track has zero curvature and the line does not
        intersect the cylinder.  If the intersections are found, and Xmax > 0, this
        function checks if the selected intersection is outside the range [-Xmax,+Xmax]
        and returns status code 5 if either of them are.  In this case the return
        value retXYZ is still valid.  retXYZ must be float[3] in the calling code;
        trackpar and Xpoint are as in PropagateToX and DistXYZ. */

    static int PropagateToCylinder(const float* trackpar,
                                   const float* Xpoint,
                                   const float rCyl,
                                   const float yCyl,
                                   const float zCyl,
                                   float* retXYZ1,
                                   float* retXYZ2,
                                   const float Xmax = 0.0,
                                   const float epsilon = 2.0e-5);

    /** Finds the point on the track in 3d, given the x value.  You might call it like this:

        TrackPropagator::PropagateToX(Track::TrackParBeg(),Track::Vertex(),float,float[3])

        where the first float is the x input and the float[3] will be the computed (x,y,z)
        point.  Or, you may call it like this:

        TrackPropagator::PropagateToX(Track::TrackParEnd(),Track::End(),float,float[3],250.0)

        but don't mix the ends!  Track::Vertex() is used to give the Xo for the parameters
        of Track::TrackParBeg() and Track::End() specifies the Xo where the track
        parameters of Track::TrackParEnd() are valid.

        Return values are 0 for good finish, 1 for straight track; if Rmax > 0 , returns
        2 if y^2 + z^2 > Rmax^2 -- means you might want to use PropagateToCylinder
        */
    static int PropagateToX(const float* trackpar,
                            const float* Xpoint,
                            const float x,
                            float* retXYZ,
                            const float Rmax = 0.0);

    /** Finds the distance from the point xyz to this track.  You may call it so:

        util::TrackPropagator::DistXYZ( Track::TrackParBeg(), Track::Vertex(), float[3], &float)

        where the float[3] is the input (x,y,z) point and last float will be the computed
        distance.  Or, you may call it like this:

        util::TrackPropagator::DistXYZ( Track::TrackParEnd(), Track::End(), float[3], &float)

        but don't mix the ends!  The 2nd argument is an initial guess for the iteration to
        the point of closest approach.

        Return values are 0 for good finish, 1 for straight track, or 2 for track
        parallel to the x axis; but the computed distance is OK for all return values.
        */
    static int DistXYZ(const float* trackpar,
                       const float* Xpoint,
                       const float* xyz,
                       float& retDist);

    /** Two methods for direction of track at an arbitrary x or phi.
        Differs from Track::FindDirectionFromTrackParameters in that that method
        only finds directions at the two ends of the track and has a different
        sign convention.  One evaluates at a value of x  and the
        other evaluates at a specific value of phi.  You may call them like so:

        util::TrackPropagator::DirectionX  (Track::TrackParBeg(),Track::Vertex(), float&,float[3])
        util::TrackPropagator::DirectionPar(Track::TrackParBeg(),Track::Vertex(), float&,float[3])

        where Track::Vertex() would specify the Xo where the track parameters of
        the first argument are valid.  The float& is either the x or phi value at
        which the direction is wanted and the float[3] will be the returned unit
        direction vector.

        The return value is 0 for good return; DirectionX can return 1 for zero
        curvature.
        */
    static int DirectionX(const float* trackpar, const float* Xpoint, float& xEval, float* retXYZ);
    static int DirectionPhi(const float* trackpar,
                            const float* Xpoint,
                            float& phiEval,
                            float* retXYZ);

  private:
    // Used by DistXYZ:
    static float d2(float xt,
                    float yt,
                    float zt,
                    float x0,
                    float yc,
                    float zc,
                    float r,
                    float s,
                    float phi,
                    float phi0);
  };

} //namespace util

#endif
