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

#include <cmath>
#include <algorithm>
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

        static const float TWO_PI;
        static const float PI;



        /** Finds an intersections of a helix with an infinite cylinder of radius 
            rCyl, centered at yCyl, zCyl and parallel to the x axis, storing the 
            "closest" intersection in retXYZ, meaning that the value of phi at the
            intersection is closest to at Xpoint.  Status codes: 0 is success; 1,2 
            and 3 are no-intersection conditions.  If the helix cylinder and given 
            cylinder are disjoint, returns 1; if one is entirely inside the other, 
            status is 2; if the cylinders have a common axes, status is 3.  Status 
            code 4 results when the phi calculation has to take the acos of a number 
            outside -1 to +1.  If the intersections are found, and Xmax > 0, this 
            function checks if the selected intersection is outside the rangec -Xmax to +Xmax, 
            and returns status code 5 if either of them are.  In this case the return 
            values retXYZ is still valid.  retXYZ must be float[3] in the calling code:
            trackpar and Xpoint are as in PropagateToX and DistXYZ.
        */
        static int PropagateToCylinder(const float* trackpar, const float* Xpoint,
                                       const float rCyl, const float yCyl, const float zCyl,
                                       float* retXYZ,    const float Xmax=0.0, const float epsilon = 2.0e-5);



         /** Finds the point on the track in 3d, given the x value.  You might call it like this:
                util::TrackPropagator::PropagateToX( Track::TrackParBeg(), Track::Vertex(),
                                                  float, float[3])
            where the first float is the x input and the float[3] will be the computed (X,y,z)
            point.  Or, you may call it like this:
                util::TrackPropagator::PropagateToX( Track::TrackParEnd(), Track::End(),
                                                  float, float[3], 250.0)
            but don't mix the ends!

            Return values are 0 for good finish, 1 for straight track; if Rmax > 0 , returns
            -1 if y^2 + z^2 > Rmax^2 -- means you might want to use PropagateToCylinder
        */
        static int PropagateToX(const float* trackpar, const float* Xpoint, const float x,
                                float* retXYZ,         const float Rmax=0.0 );



        /** Finds the distance from the point setXYZ to this track.  You may call it so:
               util::TrackPropagator::DistXYZ( Track::TrackParBeg(), Track::Vertex(), float[3],
                                                &float)
            where the float[3] is the input (x,y,z) point and last float will be the computed
            distance.  Or, you may call it like this:
               util::TrackPropagator::DistXYZ( Track::TrackParEnd(), Track::End(), float[3],
                                                &float)
            but don't mix the ends!

            Return values are 0 for good finish, 1 for straight track, or 2 for track
            parallel to the x axis; but the computed distance is OK for all return values.
        */
        static int DistXYZ(const float* trackpar, const float* Xpoint, 
                           const float* xyz,      float& retDist);
    private:
        // Used by DistXYZ:
        static float d2(float xt,float yt,float zt, float x0,float yc,float zc,
                        float r, float s, float phi,float phi0);
   };

} //namespace util

#endif
