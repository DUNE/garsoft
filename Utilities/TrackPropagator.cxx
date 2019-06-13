#include "Utilities/TrackPropagator.h"

#include "ReconstructionDataProducts/Track.h"

#include <cmath>
#include <algorithm>
#include <iostream>



namespace util {

    const float TrackPropagator::TWO_PI = static_cast<float>(2. * std::acos(-1.0));



    //----------------------------------------------------------------------
    int TrackPropagator::PropagateToCylinder(const float* trackpar, float xplane, 
        float r0, float y0, float z0, float xpar, float* retXYZ, double epsilon) {

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
            retXYZ[1] = yy1;
            retXYZ[2] = zz1;
        }
        else
        {
            retXYZ[1] = yy2;
            retXYZ[2] = zz2;
        }

        //Calculate the xpos
        const float acy = (ycc - retXYZ[1]) / radius;

        if (acy < -1 || acy > 1)  // bail out if we can't do this
        {
            retXYZ[0] = 0.0;
            retXYZ[1] = 0.0;
            retXYZ[2] = 0.0;
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

        retXYZ[0] = xpar + radius * tanl * dphi;

        //Check if xpar is bigger than endcap - if yes track is likely endcap not barrel
        if(std::fabs(retXYZ[0]) > std::fabs(xplane))
        return 5;

        return 0;
    }






    //----------------------------------------------------------------------
    int TrackPropagator::PropagateToX(const float* trackpar, const float* Xpoint,  const float x,
                                      float* retXYZ,         const float Rmax) {

        int retval = 0;
        float y, z, s;

        if (trackpar[2] == 0) {
            y = 0;
            z = 0;
            retval = 1;
        } else {

            float phi0 = trackpar[3];
            float r = 1.0/trackpar[2];
            float ZCent = trackpar[1] - r*TMath::Sin(phi0);
            float YCent = trackpar[0] + r*TMath::Cos(phi0);

            s = TMath::Tan( trackpar[4] );
            if (s != 0) { s = 1.0/s;
            }    else   { s = 1E9;  retval = 1;}

            float phi = (x -Xpoint[0]) * s / r + phi0;
            y = YCent - r * TMath::Cos(phi);
            z = ZCent + r * TMath::Sin(phi);
 
            if ( Rmax > 0) {
                if ( (retXYZ[1]*retXYZ[1] +retXYZ[2]*retXYZ[2]) > Rmax*Rmax ) retval = -1;
            }
        }

        retXYZ[0] = x;    retXYZ[1] = y;   retXYZ[2] = z;

        return retval;
    }





    //----------------------------------------------------------------------
    int TrackPropagator::DistXYZ(const float* trackpar, const float* Xpoint, const float *xyz, 
                                 float* retDist) {

        *retDist = 0;
    
        // just to make formulas more readable.  (xt,yt,zt) = test point.
        float xt = xyz[0];
        float yt = xyz[1];
        float zt = xyz[2];

        float x0 = Xpoint[0];
        float y0 = trackpar[0];
        float z0 = trackpar[1];
        float curv = trackpar[2];
        float phi0 = trackpar[3];

        float s  = TMath::Tan(trackpar[4]);
        if (s != 0) { s = 1.0/s;
        } else      { s = 1E9;  }

        float sinphi0 = TMath::Sin(phi0);
        float cosphi0 = TMath::Cos(phi0);
        float zc = trackpar[1] - sinphi0 / curv;
        float yc = trackpar[0] + cosphi0 / curv;

        if (curv == 0) {
            // distance from a point to a line -- use the norm of the cross product of the
            // unit vector along the line and the displacement between the test point and
            // a point on the line

            float h   = TMath::Sqrt(1.0 + s*s);
            float xhc = s*(yt-y0)*(cosphi0/h) - s*(sinphi0/h)*(zt-z0);
            float yhc = (zt-z0)/h             - s*(xt-x0)*(cosphi0/h);
            float zhc = s*(xt-x0)*(sinphi0/h) - (yt-y0)/h;

            *retDist = TMath::Sqrt( xhc*xhc + yhc*yhc + zhc*zhc );
            return 1;
        }

        if (s == 0) {
            // zero slope.  The track is another line, this time along the x axis

            *retDist = TMath::Sqrt( TMath::Sq(yt-y0) + TMath::Sq(zt-z0) );
            return 2;
        }

        // general case -- need to compute distance from a point to a helix
        float r = 1.0/curv;
        float span = TMath::Pi();
        float gold = 1.61803398875;

        // First guess is phi for this xt; but try also +/- a half turn as
        // that could easily be closer
        float phicent = (s*curv)*(xt-x0) +phi0;
        float d2cent  = d2(xt,yt,zt, x0,yc,zc, r,s, phicent, phi0);
        float philow  = phicent -span;
        float d2low   = d2(xt,yt,zt, x0,yc,zc, r,s, philow,  phi0);
        float phihi   = phicent +span;
        float d2hi    = d2(xt,yt,zt, x0,yc,zc, r,s, phihi,   phi0);
        if ( d2low<d2cent ) {
            // In fact, d2low==d2high at this point.  Pick one solution,
            // somewhat arbitrarily; go for the one closest to phi0.
            if ( fabs(philow -phi0)<fabs(phihi -phi0) ) {
                phicent = philow;
            } else {
                phicent = phihi;
            }
            philow = phicent -span;
            phihi  = phicent +span;
        }


        // solve for phi of the closest point using quadratic fit; 18 iters gives
        // phi to 1mrad accuracy
        for (size_t iter=0; iter<17; ++iter) {

            // 2 coefficients of the quadratic; requires (phihi -phicent) == (phicent -philow)
            float B  = (d2hi -d2low) / (2*span);
            float A  = (d2hi +d2low -2*d2cent) / ( 2*span*span );
            if ( A == 0 ) break;    // d2hi, d2low & d2cent all pretty close- just quit
            phicent += -B / (2*A);
            span    /= gold;
            philow   = phicent -span;
            phihi    = phicent +span;
            d2cent = d2(xt,yt,zt, x0,yc,zc, r,s, phicent, phi0);
            d2low  = d2(xt,yt,zt, x0,yc,zc, r,s, philow,  phi0);
            d2hi   = d2(xt,yt,zt, x0,yc,zc, r,s, phihi,   phi0);
        }

        float xp = x0 + r*(phicent -phi0)/s;
        float yp = yc - r*TMath::Cos(phicent);
        float zp = zc + r*TMath::Sin(phicent);
        *retDist = TMath::Sqrt( TMath::Sq(xt-xp) + TMath::Sq(yt-yp) + TMath::Sq(zt-zp) );

        return 0;
    }

    float TrackPropagator::d2(float xt, float yt, float zt, float x0, float yc, float zc,
                              float r, float s, float phi, float phi0) {
        float dx = (xt -x0) - (r/s)*(phi -phi0);
        float dy = (yt -yc) + r*TMath::Cos(phi);
        float dz = (zt -zc) - r*TMath::Sin(phi);
        return dx*dx + dy*dy + dz*dz;
    }



} //namespace util
