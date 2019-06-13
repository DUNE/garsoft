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



        /** STATUS: Compiles, untested
        */
        static int PropagateToCylinder(const float* trackpar, float xplane,
            float r0, float y0, float z0, float xpar, float* retXYZ, double epsilon = 1.0e-8);



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
                           const float *xyz,      float* retDist);
    private:
        // Used by DistXYZ:
        static float d2(float xt,float yt,float zt, float x0,float yc,float zc,
                        float r, float s, float phi,float phi0);
   };

} //namespace util

#endif
