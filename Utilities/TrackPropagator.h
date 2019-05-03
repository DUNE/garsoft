#ifndef TRACKPROPAGATOR_H
#define TRACKPROPAGATOR_H

////////////////////////////////////////////////////////////////////////
/// \file  TrackPropagator.h
/// \brief Interface to propagate a Track to the specific point
///
/// \version $Id:  $
/// \author  eldwan.brianne@desy.de
////////////////////////////////////////////////////////////////////////

namespace gar{
  namespace rec{
    class Track;
  }
}

namespace util
{
  class TrackPropagator {

  public:
    //No constructor
    TrackPropagator() = delete;
    //No copy constructor
    TrackPropagator(const TrackPropagator&) = delete;
    //No destructor
    ~TrackPropagator() = delete;

    /** Returns the status of the intersection (0 is success) of a track to an infinite cylinder of radius r0, centered in y0, z0 and parallel to the x axis, store the intersection coordinates in xyz
    */
    static int PropagateToCylinder(const float* trackpar, float xplane, float r0, float y0, float z0, float xpar, float* xyz, double epsilon = 1.0e-8);

    /** Returns the status of the intersection (0 is success) of a track to an infinite plane located at the xplane position perpendicular to the x axis, store the intersection coordinates in xyz
    */
    static int PropagateToXPlane(const float* trackpar, float xplane, float rmax, float xpar, float* xyz);

  private:
    static const float TWO_PI;
  };

} //namespace util

#endif
