#ifndef ClusterShapes_h
#define ClusterShapes_h

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstdlib>
#include <math.h>
#include <vector>

namespace util {

  /**
  *    Utility class to derive properties of clusters, such as centre of gravity,
  *    axes of inertia, fits of the cluster shape and so on. All the details are
  *    explained in the documentation of the methods. Several classes of the GSL
  *    (GNU Scientific Library) are needed in this class.
  *
  *    @authors V. Morgunov (ITEP/DESY), A. Raspereza (DESY), O. Wendt (DESY)
  *    @version $Id: ClusterShapes.h,v 1.14 2007-04-27 13:56:53 owendt Exp $
  *
  */
  class ClusterShapes {

  public:

    /**
    *    Constructor
    *    @param nhits : number of hits in the cluster
    *    @param a     : amplitudes of elements ('cells') of the cluster. Stored in
    *                   an array, with one entry for each element ('cell'). Each entry
    *                   is depending on coordinates x,y,z (Cartesian), which are stored
    *                   in the arrays x,y,z.
    *    @param x,y,z : array of coordinates corresponding to the array of amplitudes a.
    *
    *
    */
    ClusterShapes(int nhits, float* a, float* x, float* y, float* z);

    ClusterShapes(int nhits, float* a, float* t, float* x, float* y, float* z);

    ClusterShapes(int nhits, std::vector<float> a, std::vector<float> x, std::vector<float> y, std::vector<float> z);

    ClusterShapes(int nhits, std::vector<float> a, std::vector<float> t, std::vector<float> x, std::vector<float> y, std::vector<float> z);

    /**
    *    Destructor
    */
    ~ClusterShapes();

    /**
    * returns the number of elements of the cluster
    */
    int getNumberOfHits();

    /**
    * returns the accumulated amplitude for the whole cluster (just the sum of the
    * energy in all the entries of the cluster)
    */
    float getTotalAmplitude();

    /**
    * returns the average time over the whole cluster
    */
    float getAverageTime();

    /**
    * returns an array, which represents a vector from the origin of the
    * coordiante system, i.\ e.\ IP, to the centre of gravity of the cluster. The centre
    * of gravity is calculated with the energy of the entries of the cluster.
    */
    float* getCenterOfGravity();

    /**
    * array of the inertias of mass (i.\ e.\ energy) corresponding to the three main axes
    * of inertia. The array is sorted in ascending order.
    */
    float* getEigenValInertia();

    /**
    * array of the three main axes of inertia (9 entries) starting
    * with the axis corresponding to the smallest inertia of mass
    * (main principal axis). All axes are normalised to a length
    * of 1.
    */
    float* getEigenVecInertia();

    /**
    * 'mean' width of the cluster perpendicular to the main
    * principal axis, defined as:
    * width := sqrt( 1/TotalAmplitude * Sum(a[i]*d[i]*d[i]) ),
    * where d[i] is the distance of the i-th point to the main
    * principal axis.
    */
    float getWidth();

    /**
    * distance to the centre of gravity measured from IP
    * (absolut value of the vector to the centre of gravity)
    */
    inline float radius() { return _radius; }

    /**
    * largest spatial axis length of the ellipsoid derived
    * by the inertia tensor (by their eigenvalues and eigen-
    * vectors)
    */
    inline float getElipsoid_r1() { return _r1; }

    /**
    * medium spatial axis length of the ellipsoid derived
    * by the inertia tensor (by their eigenvalues and eigen-
    * vectors)
    */
    inline float getElipsoid_r2() { return _r2; }

    /**
    * smallest spatial axis length of the ellipsoid derived
    * by the inertia tensor (by their eigenvalues and eigen-
    * vectors)
    */
    inline float getElipsoid_r3() { return _r3; }

    /**
    * volume of the ellipsoid
    */
    inline float getElipsoid_vol() { return _vol; }

    /**
    * average radius of the ellipsoid (qubic root of volume)
    */
    inline float getElipsoid_r_ave() { return _r_ave; }

    /**
    * density of the ellipsoid defined by: totAmpl/vol
    */
    inline float getElipsoid_density() { return _density; }

    /**
    * eccentricity of the ellipsoid defined by:
    * Width/r1
    */
    inline float getElipsoid_eccentricity() { return _eccentricity; }

    /**
    * distance from centre of gravity to the point most far
    * away from IP projected on the main principal axis
    */
    inline float getElipsoid_r_forw() { return _r1_forw; }

    /**
    * distance from centre of gravity to the point nearest
    * to IP projected on the main principal axis
    */
    inline float getElipsoid_r_back() { return _r1_back; }

  private:

    int _nHits;

    std::vector<float> _aHit;
    std::vector<float> _tHit;
    std::vector<float> _xHit;
    std::vector<float> _yHit;
    std::vector<float> _zHit;

    int   _ifNotGravity = 1;
    float _totAmpl = 0.0;
    float _totTime = 0.0;
    float _radius = 0.0;
    float _xgr = 0.0;
    float _ygr = 0.0;
    float _zgr = 0.0;
    float _analogGravity[3] = {0.0, 0.0,0.0};

    int   _ifNotWidth = 1;
    float _analogWidth = 0.0;

    int   _ifNotInertia = 1;
    float _ValAnalogInertia[3];
    float _VecAnalogInertia[9];

    int   _ifNotElipsoid = 1;
    float _r1            = 0.0;  // Cluster spatial axis length -- the largest
    float _r2            = 0.0;  // Cluster spatial axis length -- less
    float _r3            = 0.0;  // Cluster spatial axis length -- less
    float _vol           = 0.0;  // Cluster ellipsoid volume
    float _r_ave         = 0.0;  // Cluster average radius  (qubic root)
    float _density       = 0.0;  // Cluster density
    float _eccentricity  = 0.0;  // Cluster Eccentricity
    float _r1_forw       = 0.0;
    float _r1_back       = 0.0;

    void  findElipsoid();
    void  findGravity();
    void  findInertia();
    void  findWidth();
    float findDistance(int i);

  };
}

#endif
