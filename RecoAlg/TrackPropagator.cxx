#include "RecoAlg/TrackPropagator.h"
#include "ReconstructionDataProducts/Track.h"

namespace util {

  //----------------------------------------------------------------------
  int TrackPropagator::PropagateToCylinder(const float* trackpar,
                                           const float* Xpoint,
                                           const float rCyl,
                                           const float yCyl,
                                           const float zCyl,
                                           float* retXYZ1,
                                           float* retXYZ2,
                                           const float Xmax,
                                           const float epsilon)
  {
    //track fitted parameters (y, z, curvature, phi, lambda)
    // Don't forget r can be negative!  But rCyl can't :)
    const float y0 = trackpar[0];
    const float z0 = trackpar[1];
    const float curv = trackpar[2];
    const float phi0 = trackpar[3];
    const float tanl = std::tan(trackpar[4]);

    // Radius of curvature of track
    float radius;
    if (curv != 0) { radius = 1.0 / curv; }
    else {
      // Straight-line algorithm not tested.  Probably not needed, really.
      // From Wolfram MathWorld
      float zL[2];
      float yL[2];
      zL[0] = zCyl - rCyl - 1.0;
      zL[1] = zCyl + rCyl + 1.0;
      for (int i = 0; i < 2; ++i)
        yL[i] = y0 + (zL[i] - z0) * std::tan(phi0);

      float Dx = zL[1] - zL[0];
      float Dy = yL[1] - yL[0];
      float Dr = std::hypot(Dx, Dy);
      float D = zL[0] * yL[1] - zL[1] * yL[0];
      float det = (rCyl * Dr) * (rCyl * Dr) - D * D;
      if (det < 0) return 4;

      float zz1 = (D * Dy + Dx * det) / (Dr * Dr);
      float zz2 = (D * Dy - Dx * det) / (Dr * Dr);
      float yy1 = y0 + (zz1 - z0) * std::tan(phi0);
      float yy2 = y0 + (zz2 - z0) * std::tan(phi0);
      if (std::hypot(zz1 - z0, yy1 - y0) < std::hypot(zz2 - z0, yy2 - y0)) {
        retXYZ1[1] = yy1;
        retXYZ1[2] = zz1;
        retXYZ2[1] = yy2;
        retXYZ2[2] = zz2;
      }
      else {
        retXYZ1[1] = yy2;
        retXYZ1[2] = zz2;
        retXYZ2[1] = yy1;
        retXYZ2[2] = zz1;
      }
      retXYZ1[0] = (retXYZ1[2] - z0) / std::cos(phi0);
      retXYZ2[0] = (retXYZ2[2] - z0) / std::cos(phi0);
      return 0;
    }

    //Coordinate of the center of the circle of radius r
    float zcc = z0 - radius * std::sin(phi0);
    float ycc = y0 + radius * std::cos(phi0);

    // dz and dy are the 'horizontal' and 'vertical' distances between
    // the circle centers.
    float dz = zcc - zCyl;
    float dy = ycc - yCyl;

    /* Determine the straight-line distance between the centers. */
    const float d = std::hypot(dy, dz);

    /* Check for solvability. */
    if (d > (rCyl + abs(radius))) {
      /* no solution. circles do not intersect. */
      return 1;
    }
    if (d < std::fabs(rCyl - abs(radius))) {
      /* no solution. one circle is contained in the other */
      return 2;
    }
    if (d < epsilon) {
      /* no solution. circles have common centre */
      return 3;
    }

    // Calculate the intersection point between the cylinder and the track
    float phiStar = (d * d + rCyl * rCyl - radius * radius);
    phiStar /= 2 * std::max(1.e-20f, rCyl * d);
    if (phiStar > +1.f) phiStar = +0.9999999f;
    if (phiStar < -1.f) phiStar = -0.9999999f;
    phiStar = std::acos(phiStar);

    float phiCentre(std::atan2(dy, dz));

    // There are two solutions which can not differ in phi from
    // phi0 by more than 2pi.
    float zz1 = zCyl + rCyl * std::cos(phiCentre + phiStar);
    float yy1 = yCyl + rCyl * std::sin(phiCentre + phiStar);
    // Add pi/2 or 3pi/2 because of definition of where phi = 0
    float phi1 = atan2(yy1 - ycc, zz1 - zcc);
    phi1 += (curv >= 0) ? M_PI / 2.0 : 3.0 * M_PI / 2.0;
    float zz2 = zCyl + rCyl * std::cos(phiCentre - phiStar);
    float yy2 = yCyl + rCyl * std::sin(phiCentre - phiStar);
    float phi2 = atan2(yy2 - ycc, zz2 - zcc);
    phi2 += (curv >= 0) ? M_PI / 2.0 : 3.0 * M_PI / 2.0;

    float dphi1 = phi1 - phi0;
    float dphi2 = phi2 - phi0;
    int nturns = dphi1 / (2.0 * M_PI);
    // bring dphi1 in the range 0-2pi, then into range -pi to +pi
    if (dphi1 > +(2.0 * M_PI)) dphi1 -= (2.0 * M_PI) * nturns;
    if (dphi1 < -(2.0 * M_PI)) dphi1 += (2.0 * M_PI) * nturns;
    if (dphi1 > +M_PI) dphi1 -= 2.0 * M_PI;
    if (dphi1 < -M_PI) dphi1 += 2.0 * M_PI;
    // likewise dphi2
    nturns = dphi2 / (2.0 * M_PI);
    if (dphi2 > +(2.0 * M_PI)) dphi2 -= (2.0 * M_PI) * nturns;
    if (dphi2 < -(2.0 * M_PI)) dphi2 += (2.0 * M_PI) * nturns;
    if (dphi2 > +M_PI) dphi2 -= 2.0 * M_PI;
    if (dphi2 < -M_PI) dphi2 += 2.0 * M_PI;
    if (fabs(dphi1) < fabs(dphi2)) {
      retXYZ1[1] = yy1;
      retXYZ1[2] = zz1;
      retXYZ1[0] = Xpoint[0] + radius * tanl * dphi1;
      retXYZ2[1] = yy2;
      retXYZ2[2] = zz2;
      retXYZ2[0] = Xpoint[0] + radius * tanl * dphi2;
    }
    else {
      retXYZ1[1] = yy2;
      retXYZ1[2] = zz2;
      retXYZ1[0] = Xpoint[0] + radius * tanl * dphi2;
      retXYZ2[1] = yy1;
      retXYZ2[2] = zz1;
      retXYZ2[0] = Xpoint[0] + radius * tanl * dphi1;
    }

    // If Xmax > 0, check if found intersection has x beyond it - i.e.
    // it might be in the endcap
    if (std::fabs(Xmax) > 0) {
      if (std::fabs(retXYZ1[0]) > std::fabs(Xmax)) {
        // std::cout << "5 - x > Xmax, might be in the endcap" << std::endl;
        return 5;
      }
      if (std::fabs(retXYZ2[0]) > std::fabs(Xmax)) { return 5; }
    }

    return 0;
  }

  //----------------------------------------------------------------------
  int TrackPropagator::PropagateToX(const float* trackpar,
                                    const float* Xpoint,
                                    const float x,
                                    float* retXYZ,
                                    const float Rmax)
  {
    int retval = 0;
    float y, z;

    const float y0 = trackpar[0];
    const float z0 = trackpar[1];
    const float phi0 = trackpar[3];
    const float curv = trackpar[2];
    float radius = 0;
    if (curv != 0) radius = 1.0 / curv;

    float ZCent = trackpar[1] - radius * std::sin(phi0);
    float YCent = trackpar[0] + radius * std::cos(phi0);
    float tanl = std::tan(trackpar[4]);

    if (radius == 0) {
      // Straight-line case not tested yet
      float dT = (x - Xpoint[0]) / tanl;
      y = y0 + dT * std::sin(phi0);
      z = z0 + dT * std::cos(phi0);
      retval = 1;
    }
    else {
      float s;
      if (tanl != 0) { s = 1.0 / tanl; }
      else {
        s = 1E9;
        retval = 1;
      }

      float phi = (x - Xpoint[0]) * s / radius + phi0;
      y = YCent - radius * std::cos(phi);
      z = ZCent + radius * std::sin(phi);
    }

    if (Rmax > 0) {
      if (std::hypot(z, y) > Rmax) retval = 2;
    }

    retXYZ[0] = x;
    retXYZ[1] = y;
    retXYZ[2] = z;

    return retval;
  }

  //----------------------------------------------------------------------
  int TrackPropagator::DistXYZ(const float* trackpar,
                               const float* Xpoint,
                               const float* xyz,
                               float& retDist)
  {

    retDist = 0;

    // just to make formulas more readable.  (xt,yt,zt) = test point.
    float xt = xyz[0];
    float yt = xyz[1];
    float zt = xyz[2];

    float x0 = Xpoint[0];
    float y0 = trackpar[0];
    float z0 = trackpar[1];
    float curv = trackpar[2];
    float r = 1.0 / curv;
    float phi0 = trackpar[3];

    float s = std::tan(trackpar[4]);
    if (s != 0) { s = 1.0 / s; }
    else {
      s = 1E9;
    }

    float sinphi0 = std::sin(phi0);
    float cosphi0 = std::cos(phi0);
    float zc = trackpar[1] - r * sinphi0;
    float yc = trackpar[0] + r * cosphi0;

    if (curv == 0) {
      // distance from a point to a line -- use the norm of the cross product of the
      // unit vector along the line and the displacement between the test point and
      // a point on the line

      float h = std::sqrt(1.0 + s * s);
      float xhc = s * (yt - y0) * (cosphi0 / h) - s * (sinphi0 / h) * (zt - z0);
      float yhc = (zt - z0) / h - s * (xt - x0) * (cosphi0 / h);
      float zhc = s * (xt - x0) * (sinphi0 / h) - (yt - y0) / h;

      retDist = std::sqrt(xhc * xhc + yhc * yhc + zhc * zhc);
      return 1;
    }

    if (s == 0) {
      // zero slope.  The track is another line, this time along the x axis
      retDist = std::sqrt(TMath::Sq(yt - y0) + TMath::Sq(zt - z0));
      return 2;
    }

    // general case -- need to compute distance from a point to a helix
    float span = M_PI;
    float gold = 1.61803398875;

    // First guess is phi for this xt; but try also +/- a half turn as
    // that could easily be closer
    float phicent = (s * curv) * (xt - x0) + phi0;
    float d2cent = d2(xt, yt, zt, x0, yc, zc, r, s, phicent, phi0);
    float philow = phicent - span;
    float d2low = d2(xt, yt, zt, x0, yc, zc, r, s, philow, phi0);
    float phihi = phicent + span;
    float d2hi = d2(xt, yt, zt, x0, yc, zc, r, s, phihi, phi0);
    if (d2low < d2cent) {
      // In fact, d2low==d2high at this point.  Pick one solution,
      // somewhat arbitrarily; go for the one closest to phi0.
      if (std::fabs(philow - phi0) < std::fabs(phihi - phi0)) { phicent = philow; }
      else {
        phicent = phihi;
      }
      philow = phicent - span;
      phihi = phicent + span;
    }

    // solve for phi of the closest point using quadratic fit; 18 iters gives
    // phi to 1mrad accuracy
    for (size_t iter = 0; iter < 18; ++iter) {
      // 2 coefficients of the quadratic; requires (phihi -phicent) == (phicent -philow)
      float B = (d2hi - d2low) / (2 * span);
      float A = (d2hi + d2low - 2 * d2cent) / (2 * span * span);
      if (A == 0) break; // d2hi, d2low & d2cent all pretty close- just quit
      phicent += -B / (2 * A);
      span /= gold;
      philow = phicent - span;
      phihi = phicent + span;
      d2cent = d2(xt, yt, zt, x0, yc, zc, r, s, phicent, phi0);
      d2low = d2(xt, yt, zt, x0, yc, zc, r, s, philow, phi0);
      d2hi = d2(xt, yt, zt, x0, yc, zc, r, s, phihi, phi0);
    }

    float xp = x0 + r * (phicent - phi0) / s;
    float yp = yc - r * TMath::Cos(phicent);
    float zp = zc + r * TMath::Sin(phicent);
    retDist = std::sqrt(TMath::Sq(xt - xp) + TMath::Sq(yt - yp) + TMath::Sq(zt - zp));

    return 0;
  }

  //----------------------------------------------------------------------
  int TrackPropagator::DirectionX(const float* trackpar,
                                  const float* Xpoint,
                                  float& xEval,
                                  float* retXYZ)
  {

    if (trackpar[2] == 0) return 1;
    float r = 1.0 / trackpar[2];
    float s = std::tan(trackpar[4]);
    float x0 = Xpoint[0];

    float phi = (xEval - x0) / (r * s);
    DirectionPhi(trackpar, Xpoint, phi, retXYZ);
    return 0;
  }

  int TrackPropagator::DirectionPhi(const float* trackpar,
                                    const float*, // (unused) Xpoint,
                                    float& phiEval,
                                    float* retXYZ)
  {

    retXYZ[0] = std::tan(trackpar[4]);
    retXYZ[1] = std::sin(phiEval);
    retXYZ[2] = std::cos(phiEval);
    float norm = std::hypot(retXYZ[0], retXYZ[1], retXYZ[2]);
    for (int i = 0; i < 3; ++i)
      retXYZ[i] /= norm;

    return 0;
  }

  //----------------------------------------------------------------------
  float TrackPropagator::d2(float xt,
                            float yt,
                            float zt,
                            float x0,
                            float yc,
                            float zc,
                            float r,
                            float s,
                            float phi,
                            float phi0)
  {
    float dx = (xt - x0) - (r / s) * (phi - phi0);
    float dy = (yt - yc) + r * TMath::Cos(phi);
    float dz = (zt - zc) - r * TMath::Sin(phi);
    return dx * dx + dy * dy + dz * dz;
  }

} //namespace util
