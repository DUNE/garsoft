#include "Utilities/TrackPropagator.h"

#include "ReconstructionDataProducts/Track.h"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace util {

    const float TrackPropagator::TWO_PI = static_cast<float>(2. * std::acos(-1.0));

    //----------------------------------------------------------------------
    int TrackPropagator::PropagateToCylinder(const float* trackpar, float xplane, float r0, float y0, float z0, float xpar, float* xyz, double epsilon)
    {
        //track fitted parameters (y, z, curvature, phi, lambda)
        const float yref = trackpar[0];
        const float zref = trackpar[1];
        const float omega = trackpar[2];
        const float phi0 = trackpar[3];
        const float tanl = std::tan(trackpar[4]);

        // Radius of curvature of track
        float radius = 0;
        if (omega != 0) radius = 1.0 / std::fabs(omega);

        //Coordinate of the center of the circle of radius r
        const float ycc = yref + radius * std::cos(phi0);
        const float zcc = zref - radius * std::sin(phi0);

        /* dx and dy are the vertical and horizontal distances between
        * the circle centers.
        */
        const float dy = ycc - y0;
        const float dz = zcc - z0;

        /* Determine the straight-line distance between the centers. */
        const float d = std::hypot(dy, dz);

        /* Check for solvability. */
        if (d > (r0 + radius))
        {
            /* no solution. circles do not intersect. */
            return 1;
        }
        if (d < fabs(r0 - radius))
        {
            /* no solution. one circle is contained in the other */
            return 2;
        }
        if (d < epsilon)
        {
            /* no solution. circles have common centre */
            return 3;
        }

        //There is solution
        //Calculate the intersection point between the cylinder and the track
        float phiStar = (d*d + r0*r0 - radius*radius);
        phiStar /= 2 * std::max(1.e-20f, r0 * d);

        if (phiStar > 1.f)
        phiStar = 0.9999999f;

        if (phiStar < -1.f)
        phiStar = -0.9999999f;

        phiStar = std::acos(phiStar);
        const float phiCentre(std::atan2(dy, dz));

        // two solutions -- pick the one closest to a linear extrapolation from the point where we are
        //First intersection
        const float yy1 = y0 + r0 * std::sin(phiCentre + phiStar);
        const float zz1 = z0 + r0 * std::cos(phiCentre + phiStar);

        //Second intersection
        const float yy2 = y0 + r0 * std::sin(phiCentre-phiStar);
        const float zz2 = z0 + r0 * std::cos(phiCentre-phiStar);

        float dir[3];
        gar::rec::FindDirectionFromTrackParameters(trackpar, dir);
        const float t1 = (yy1-yref) * (-dir[1]) + (zz1-zref) * (-dir[2]);
        const float t2 = (yy2-yref) * (-dir[1]) + (zz2-zref) * (-dir[2]);

        if (t1 > t2)
        {
            xyz[1] = yy1;
            xyz[2] = zz1;
        }
        else
        {
            xyz[1] = yy2;
            xyz[2] = zz2;
        }

        //Calculate the xpos
        const float acy = (ycc - xyz[1]) / radius;

        if (acy < -1 || acy > 1)  // bail out if we can't do this
        {
            xyz[0] = 0.;
            xyz[1] = 0.;
            xyz[2] = 0.;
            return 4;
        }

        float dphi = std::acos(acy) - phi0;

        if (dphi > TWO_PI)
        {
            int nturns = dphi/TWO_PI;
            dphi = dphi - TWO_PI*nturns;
        }
        if (dphi < -TWO_PI)
        {
            int nturns = -dphi/TWO_PI;
            dphi = dphi + TWO_PI*nturns;
        }

        xyz[0] = xpar + radius * tanl * dphi;

        //Check if xpar is bigger than endcap - if yes track is likely endcap not barrel
        if(std::fabs(xyz[0]) > std::fabs(xplane))
        return 5;

        return 0;
    }

    //----------------------------------------------------------------------
    int TrackPropagator::PropagateToXPlane(const float* trackpar, float xplane, float rmax, float xpar, float* xyz)
    {
        const float yref = trackpar[0];
        const float zref = trackpar[1];
        const float omega = trackpar[2];
        const float phi0 = trackpar[3];
        const float tanl = std::tan(trackpar[4]);

        // Radius of curvature of track
        float radius = 0;
        if (omega != 0) radius = 1.0 / std::fabs(omega);

        //Coordinate of the center of the circle of radius r
        const float ycc = yref + radius * std::cos(phi0);
        const float zcc = zref - radius * std::sin(phi0);

        // get path length to crossing point
        const float xref = xpar;
        const float s = ( xplane - xref ) / tanl;
        const float delta_phi_half = (omega * s)/2.0 ;

        xyz[0] = xplane;
        xyz[1] = ycc + radius * std::sin(phi0 - delta_phi_half);
        xyz[2] = zcc + radius * std::cos(phi0 - delta_phi_half);

        if(std::sqrt( xyz[1]*xyz[1] + xyz[2]*xyz[2] ) > rmax)
        return 1;

        return 0;
    }

} //namespace util
