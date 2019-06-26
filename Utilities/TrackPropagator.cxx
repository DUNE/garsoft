#include "Utilities/TrackPropagator.h"
#include "ReconstructionDataProducts/Track.h"



namespace util {

    const float TrackPropagator::TWO_PI = static_cast<float>(2. * std::acos(-1.0));



    //----------------------------------------------------------------------
    int TrackPropagator::PropagateToCylinder(const float* trackpar, const float* Xpoint,
                                             const float rCyl, const float yCyl, const float zCyl,
                                             float* retXYZ,    const float Xmax, const float epsilon) {

        //track fitted parameters (y, z, curvature, phi, lambda)
        const float y0   = trackpar[0];
        const float z0   = trackpar[1];
        const float curv = trackpar[2];
        const float phi0 = trackpar[3];
        const float tanl = std::tan(trackpar[4]);

        // Radius of curvature of track
        float radius = 0;
        if (curv != 0) radius = 1.0 / curv;

        //Coordinate of the center of the circle of radius r
        const float zcc = z0 - radius * std::sin(phi0);
        const float ycc = y0 + radius * std::cos(phi0);
 
        // dz and dy are the 'horizontal' and 'vertical' distances between
        // the circle centers.
        const float dz = zcc - zCyl;
        const float dy = ycc - yCyl;

        /* Determine the straight-line distance between the centers. */
        const float d = std::hypot(dy, dz);



        /* Check for solvability. */
        if ( d > (rCyl + radius) ) {
            /* no solution. circles do not intersect. */
            return 1;
        }
        if (d < std::fabs(rCyl - radius)) {
            /* no solution. one circle is contained in the other */
            return 2;
        }
        if (d < epsilon) {
            /* no solution. circles have common centre */
            return 3;
        }



        //There is solution
        //Calculate the intersection point between the cylinder and the track
        float phiStar = (d*d + rCyl*rCyl - radius*radius);
        phiStar /= 2 * std::max(1.e-20f, rCyl * d);
        if (phiStar > +1.f) phiStar = 0.9999999f;
        if (phiStar < -1.f) phiStar = -0.9999999f;
        phiStar = std::acos(phiStar);

        const float phiCentre(std::atan2(dy, dz));

        // two solutions -- pick the one closest to a linear extrapolation from the point where we are
        //First intersection (compute using eqn of cylinder, not helix!
        const float zz1 = zCyl + rCyl * std::cos(phiCentre + phiStar);
        const float yy1 = yCyl + rCyl * std::sin(phiCentre + phiStar);
 
        //Second intersection
        const float zz2 = zCyl + rCyl * std::cos(phiCentre-phiStar);
        const float yy2 = yCyl + rCyl * std::sin(phiCentre-phiStar);

        // Select which intersection differently for track originating inside 
        // cylinder (usual case) vs outside cylinder (could happen!)
        float dir[3];    float t1,t2;
        gar::rec::FindDirectionFromTrackParameters(trackpar, dir);
        if ( std::hypot( z0-zCyl, y0-yCyl) <= rCyl ) {
            // originates inside
            t1 = (yy1-y0) * dir[1] + (zz1-z0) * dir[2];
            t2 = (yy2-y0) * dir[1] + (zz2-z0) * dir[2];
        } else {
            // Originates outside
            t1 = (yy2-yy1) * dir[1] + (zz2-zz1) * dir[2];
            t2 = -t1;
        }
        if (t1 > t2) {
            retXYZ[1] = yy1;
            retXYZ[2] = zz1;
        } else {
            retXYZ[1] = yy2;
            retXYZ[2] = zz2;
        }

        //Calculate the xpos
        const float acz = retXYZ[2] -zcc;
        const float acy = retXYZ[1] -ycc;
        const float phiInt = atan2(acy, acz) +TWO_PI/4.0;    // Need the +/- range of atan2

        // Select the intersection within 2pi of phi0
        float dphi = phiInt - phi0;
        if (dphi > TWO_PI) {
            int nturns = dphi/TWO_PI;
            dphi = dphi - TWO_PI*nturns;
        }
        if (dphi < -TWO_PI) {
            int nturns = -dphi/TWO_PI;
            dphi = dphi + TWO_PI*nturns;
        }

        retXYZ[0] = Xpoint[0] + radius * tanl * dphi;

        // If Xmax > 0, check if found intersection has x beyond it - i.e.
        // it might be in the endcap
        if (fabs(Xmax)>0) {
            if ( std::fabs(retXYZ[0]) > std::fabs(Xmax) ) return 5;
        }
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
            float ZCent = trackpar[1] - r*std::sin(phi0);
            float YCent = trackpar[0] + r*std::cos(phi0);

            s = TMath::Tan( trackpar[4] );
            if (s != 0) { s = 1.0/s;
            }    else   { s = 1E9;  retval = 1;}

            float phi = (x -Xpoint[0]) * s / r + phi0;
            y = YCent - r * std::cos(phi);
            z = ZCent + r * std::sin(phi);
 
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
        float r = 1.0/curv;
        float phi0 = trackpar[3];

        float s  = std::tan(trackpar[4]);
        if (s != 0) { s = 1.0/s;
        } else      { s = 1E9;  }

        float sinphi0 = std::sin(phi0);
        float cosphi0 = std::cos(phi0);
        float zc = trackpar[1] - r * sinphi0;
        float yc = trackpar[0] + r * cosphi0;

        if (curv == 0) {
            // distance from a point to a line -- use the norm of the cross product of the
            // unit vector along the line and the displacement between the test point and
            // a point on the line

            float h   = std::sqrt(1.0 + s*s);
            float xhc = s*(yt-y0)*(cosphi0/h) - s*(sinphi0/h)*(zt-z0);
            float yhc = (zt-z0)/h             - s*(xt-x0)*(cosphi0/h);
            float zhc = s*(xt-x0)*(sinphi0/h) - (yt-y0)/h;

            *retDist = std::sqrt( xhc*xhc + yhc*yhc + zhc*zhc );
            return 1;
        }

        if (s == 0) {
            // zero slope.  The track is another line, this time along the x axis

            *retDist = std::sqrt( TMath::Sq(yt-y0) + TMath::Sq(zt-z0) );
            return 2;
        }

        // general case -- need to compute distance from a point to a helix
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
            if ( std::fabs(philow -phi0) < std::fabs(phihi -phi0) ) {
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
        *retDist = std::sqrt( TMath::Sq(xt-xp) + TMath::Sq(yt-yp) + TMath::Sq(zt-zp) );

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
