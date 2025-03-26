#include "GArObjects/Helix.h"

#include <cmath>
#include <iostream>
#include <limits>

namespace gar_content {

  const float Helix::FCT = 2.99792458E-4f;
  const float Helix::TWO_PI = static_cast<float>(2. * std::acos(-1.0));
  const float Helix::HALF_PI = static_cast<float>(0.5 * std::acos(-1.0));

  //------------------------------------------------------------------------------------------------------------------------------------------

  Helix::Helix(const float phi0,
               const float d0,
               const float x0,
               const float omega,
               const float tanLambda,
               const float bField)
    : m_referencePoint(0.f, 0.f, 0.f)
    , m_momentum(0.f, 0.f, 0.f)
    , m_phi0(phi0)
    , m_d0(d0)
    , m_x0(x0)
    , m_omega(omega)
    , m_tanLambda(tanLambda)
  {
    if ((bField < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(omega) < std::numeric_limits<float>::epsilon())) {
      std::cout << "Helix, invalid parameter: bField " << bField << ", omega " << omega
                << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    m_charge = omega / std::fabs(omega);
    m_radius = 1.f / std::fabs(omega);
    m_zAtPCA = -m_d0 * std::sin(m_phi0);
    m_yAtPCA = m_d0 * std::cos(m_phi0);
    m_referencePoint.SetValues(m_zAtPCA, m_yAtPCA, m_x0);
    m_pzy = FCT * bField * m_radius;
    m_momentum.SetValues(m_pzy * std::cos(m_phi0), m_pzy * std::sin(m_phi0), m_tanLambda * m_pzy);
    m_pzAtPCA = m_momentum.GetZ();
    m_pyAtPCA = m_momentum.GetY();
    m_phiMomRefPoint = std::atan2(m_pyAtPCA, m_pzAtPCA);
    m_zCentre = m_referencePoint.GetZ() + m_radius * std::cos(m_phi0 - HALF_PI * m_charge);
    m_yCentre = m_referencePoint.GetY() + m_radius * std::sin(m_phi0 - HALF_PI * m_charge);
    m_phiAtPCA = std::atan2(-m_yCentre, -m_zCentre);
    m_phiRefPoint = m_phiAtPCA;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  Helix::Helix(const pandora::CartesianVector& position,
               const pandora::CartesianVector& momentum,
               const float charge,
               const float bField)
    : m_referencePoint(position), m_momentum(momentum), m_charge(charge)
  {
    const double pz(momentum.GetZ()), py(momentum.GetY());
    const double pzy(std::sqrt(pz * pz + py * py));

    if ((bField < std::numeric_limits<float>::epsilon()) ||
        (std::fabs(pzy) < std::numeric_limits<float>::epsilon())) {
      std::cout << "Helix, invalid parameter: bField " << bField << ", pzy " << pzy << std::endl;
      throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    const double radius(pzy / (FCT * bField));
    m_pzy = static_cast<float>(pzy);
    m_radius = static_cast<float>(radius);

    m_omega = charge / m_radius;
    m_tanLambda = momentum.GetX() / m_pzy;
    m_phiMomRefPoint = static_cast<float>(std::atan2(py, pz));

    const double z(position.GetZ()), y(position.GetY());
    const double zCentre(z + radius *
                               static_cast<double>(std::cos(m_phiMomRefPoint - HALF_PI * charge)));
    const double yCentre(y + radius *
                               static_cast<double>(std::sin(m_phiMomRefPoint - HALF_PI * charge)));
    m_zCentre = static_cast<float>(zCentre);
    m_yCentre = static_cast<float>(yCentre);

    double d0;
    if (charge > 0) {
      d0 = static_cast<double>(charge) * radius -
           static_cast<double>(std::sqrt(zCentre * zCentre + yCentre * yCentre));
    }
    else {
      d0 = static_cast<double>(charge) * radius +
           static_cast<double>(std::sqrt(zCentre * zCentre + yCentre * yCentre));
    }
    m_d0 = static_cast<float>(d0);

    m_phiRefPoint = static_cast<float>(std::atan2(y - yCentre, z - zCentre));
    m_phiAtPCA = static_cast<float>(std::atan2(-yCentre, -zCentre));
    m_phi0 = -HALF_PI * charge + m_phiAtPCA;

    while (m_phi0 < 0.)
      m_phi0 += TWO_PI;

    while (m_phi0 >= TWO_PI)
      m_phi0 -= TWO_PI;

    m_zAtPCA = m_zCentre + m_radius * std::cos(m_phiAtPCA);
    m_yAtPCA = m_yCentre + m_radius * std::sin(m_phiAtPCA);
    m_pzAtPCA = m_pzy * std::cos(m_phi0);
    m_pyAtPCA = m_pzy * std::sin(m_phi0);

    float deltaPhi = m_phiRefPoint - m_phiAtPCA;
    float zCircles = (-position.GetX() * charge) / (TWO_PI * ((m_radius * m_tanLambda) - deltaPhi));

    int n1, n2;
    if (zCircles >= std::numeric_limits<float>::epsilon()) {
      n1 = static_cast<int>(zCircles);
      n2 = n1 + 1;
    }
    else {
      n1 = static_cast<int>(zCircles) - 1;
      n2 = n1 + 1;
    }

    const int nCircles((std::fabs(n1 - zCircles) < std::fabs(n2 - zCircles) ? n1 : n2));
    m_x0 = position.GetX() + m_radius * m_tanLambda * charge * (deltaPhi + TWO_PI * nCircles);
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode Helix::GetPointInZY(const float z0,
                                          const float y0,
                                          const float az,
                                          const float ay,
                                          const pandora::CartesianVector& referencePoint,
                                          pandora::CartesianVector& intersectionPoint,
                                          float& genericTime) const
  {
    const float AA(std::sqrt(az * az + ay * ay));

    if (AA <= 0) return pandora::STATUS_CODE_FAILURE;

    const float BB((az * (z0 - m_zCentre) + ay * (y0 - m_yCentre)) / AA);
    const float CC(((z0 - m_zCentre) * (z0 - m_zCentre) + (y0 - m_yCentre) * (y0 - m_yCentre) -
                    m_radius * m_radius) /
                   AA);

    const float DET(BB * BB - CC);

    if (DET < 0) return pandora::STATUS_CODE_NOT_FOUND;

    float tt1(-BB + std::sqrt(DET));
    float tt2(-BB - std::sqrt(DET));

    const float zz1(z0 + tt1 * az);
    const float yy1(y0 + tt1 * ay);
    const float zz2(z0 + tt2 * az);
    const float yy2(y0 + tt2 * ay);

    const float phi1(std::atan2(yy1 - m_yCentre, zz1 - m_zCentre));
    const float phi2(std::atan2(yy2 - m_yCentre, zz2 - m_zCentre));
    const float phi0(
      std::atan2(referencePoint.GetY() - m_yCentre, referencePoint.GetZ() - m_zCentre));

    float dphi1(phi1 - phi0);
    float dphi2(phi2 - phi0);

    if (dphi1 < 0 && m_charge < 0) { dphi1 = dphi1 + TWO_PI; }
    else if (dphi1 > 0 && m_charge > 0) {
      dphi1 = dphi1 - TWO_PI;
    }

    if (dphi2 < 0 && m_charge < 0) { dphi2 = dphi2 + TWO_PI; }
    else if (dphi2 > 0 && m_charge > 0) {
      dphi2 = dphi2 - TWO_PI;
    }

    // Calculate generic time
    tt1 = -m_charge * dphi1 * m_radius / m_pzy;
    tt2 = -m_charge * dphi2 * m_radius / m_pzy;

    if ((tt1 < 0.) || (tt2 < 0.))
      std::cout << "Helix:: GetPointInXY, warning - negative generic time, tt1 " << tt1 << ", tt2 "
                << tt2 << std::endl;

    if (tt1 < tt2) {
      genericTime = tt1;
      intersectionPoint.SetValues(
        referencePoint.GetX() + genericTime * m_momentum.GetX(), yy1, zz1);
    }
    else {
      genericTime = tt2;
      intersectionPoint.SetValues(
        referencePoint.GetX() + genericTime * m_momentum.GetX(), yy2, zz2);
    }

    return pandora::STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode Helix::GetPointInX(const float xPlane,
                                         const pandora::CartesianVector& referencePoint,
                                         pandora::CartesianVector& intersectionPoint,
                                         float& genericTime) const
  {
    if (std::fabs(m_momentum.GetX()) < std::numeric_limits<float>::epsilon())
      return pandora::STATUS_CODE_NOT_FOUND;

    genericTime = (xPlane - referencePoint.GetX()) / m_momentum.GetX();

    const float phi0(
      std::atan2(referencePoint.GetY() - m_yCentre, referencePoint.GetZ() - m_zCentre));
    const float phi(phi0 - m_charge * m_pzy * genericTime / m_radius);

    intersectionPoint.SetValues(
      xPlane, m_yCentre + m_radius * std::sin(phi), m_zCentre + m_radius * std::cos(phi));

    return pandora::STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode Helix::GetPointOnCircle(const float radius,
                                              const pandora::CartesianVector& referencePoint,
                                              pandora::CartesianVector& intersectionPoint,
                                              float& genericTime) const
  {
    const float distCenterToIP(std::sqrt(m_zCentre * m_zCentre + m_yCentre * m_yCentre));

    if (((distCenterToIP + m_radius) < radius) || ((m_radius + radius) < distCenterToIP))
      return pandora::STATUS_CODE_NOT_FOUND;

    const float phiCentre(std::atan2(m_yCentre, m_zCentre));

    float phiStar(radius * radius + distCenterToIP * distCenterToIP - m_radius * m_radius);
    phiStar = 0.5f * phiStar / std::max(1.e-20f, radius * distCenterToIP);

    if (phiStar > 1.f) phiStar = 0.9999999f;

    if (phiStar < -1.f) phiStar = -0.9999999f;

    phiStar = std::acos(phiStar);

    const float zz1(radius * std::cos(phiCentre + phiStar));
    const float yy1(radius * std::sin(phiCentre + phiStar));

    const float zz2(radius * std::cos(phiCentre - phiStar));
    const float yy2(radius * std::sin(phiCentre - phiStar));

    const float phi1(std::atan2(yy1 - m_yCentre, zz1 - m_zCentre));
    const float phi2(std::atan2(yy2 - m_yCentre, zz2 - m_zCentre));
    const float phi0(
      std::atan2(referencePoint.GetY() - m_yCentre, referencePoint.GetZ() - m_zCentre));

    float dphi1(phi1 - phi0);
    float dphi2(phi2 - phi0);

    if (dphi1 < 0 && m_charge < 0) { dphi1 = dphi1 + TWO_PI; }
    else if (dphi1 > 0 && m_charge > 0) {
      dphi1 = dphi1 - TWO_PI;
    }

    if (dphi2 < 0 && m_charge < 0) { dphi2 = dphi2 + TWO_PI; }
    else if (dphi2 > 0 && m_charge > 0) {
      dphi2 = dphi2 - TWO_PI;
    }

    // Calculate generic time
    const float tt1(-m_charge * dphi1 * m_radius / m_pzy);
    const float tt2(-m_charge * dphi2 * m_radius / m_pzy);

    if ((tt1 < 0.) || (tt2 < 0.))
      std::cout << "Helix:: GetPointOnCircle, warning - negative generic time, tt1 " << tt1
                << ", tt2 " << tt2 << std::endl;

    // Previously returned both xx1, xx2, etc.
    if (tt1 < tt2) {
      genericTime = tt1;
      intersectionPoint.SetValues(
        referencePoint.GetX() + genericTime * m_momentum.GetX(), yy1, zz1);
    }
    else {
      genericTime = tt2;
      intersectionPoint.SetValues(
        referencePoint.GetX() + genericTime * m_momentum.GetX(), yy2, zz2);
    }

    return pandora::STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode Helix::GetDistanceToPoint(const pandora::CartesianVector& point,
                                                pandora::CartesianVector& distance,
                                                float& genericTime) const
  {
    const float phi(std::atan2(point.GetY() - m_yCentre, point.GetZ() - m_zCentre));
    const float phi0(
      std::atan2(m_referencePoint.GetY() - m_yCentre, m_referencePoint.GetZ() - m_zCentre));

    int nCircles = 0;
    if (std::fabs(m_tanLambda * m_radius) > 1.e-20) {
      const float zCircles(
        (phi0 - phi -
         m_charge * (point.GetX() - m_referencePoint.GetX()) / (m_tanLambda * m_radius)) /
        TWO_PI);

      int n1, n2;
      if (zCircles >= std::numeric_limits<float>::epsilon()) {
        n1 = static_cast<int>(zCircles);
        n2 = n1 + 1;
      }
      else {
        n1 = static_cast<int>(zCircles) - 1;
        n2 = n1 + 1;
      }

      nCircles = ((std::fabs(n1 - zCircles) < std::fabs(n2 - zCircles) ? n1 : n2));
    }

    const float dPhi(TWO_PI * (static_cast<float>(nCircles)) + phi - phi0);
    const float xOnHelix(m_referencePoint.GetX() - m_charge * m_radius * m_tanLambda * dPhi);

    const float distZ(std::fabs(m_zCentre - point.GetZ()));
    const float distY(std::fabs(m_yCentre - point.GetY()));
    const float distX(std::fabs(xOnHelix - point.GetX()));

    float distZY(std::sqrt(distZ * distZ + distY * distY));
    distZY = std::fabs(distZY - m_radius);

    distance.SetValues(distX, distZY, std::sqrt(distZY * distZY + distZ * distZ));

    if (std::fabs(m_momentum.GetX()) > 0) {
      genericTime = (xOnHelix - m_referencePoint.GetX()) / m_momentum.GetX();
    }
    else {
      genericTime = m_charge * m_radius * dPhi / m_pzy;
    }

    return pandora::STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::StatusCode Helix::GetDistanceToHelix(const gar_content::Helix* const pHelix,
                                                pandora::CartesianVector& positionOfClosestApproach,
                                                pandora::CartesianVector& v0momentum,
                                                float& helixDistance) const
  {
    if (this == pHelix) return pandora::STATUS_CODE_INVALID_PARAMETER;

    const float z01(this->GetZCentre());
    const float y01(this->GetYCentre());

    const float z02(pHelix->GetZCentre());
    const float y02(pHelix->GetYCentre());

    const float rad1(this->GetRadius());
    const float rad2(pHelix->GetRadius());

    const float distance(std::sqrt((z01 - z02) * (z01 - z02) + (y01 - y02) * (y01 - y02)));

    bool singlePoint(true);
    float phi1(0.), phi2(0.);

    if (rad1 + rad2 < distance) {
      phi1 = std::atan2(y02 - y01, z02 - z01);
      phi2 = std::atan2(y01 - y02, z01 - z02);
    }
    else if (distance + rad2 < rad1) {
      phi1 = std::atan2(y02 - y01, z02 - z01);
      phi2 = phi1;
    }
    else if (distance + rad1 < rad2) {
      phi1 = std::atan2(y01 - y02, z01 - z02);
      phi2 = phi1;
    }
    else {
      if ((std::fabs(distance) < std::numeric_limits<float>::epsilon()) ||
          (std::fabs(rad2) < std::numeric_limits<float>::epsilon()))
        return pandora::STATUS_CODE_FAILURE;

      singlePoint = false;
      float cosAlpha = 0.5f * (distance * distance + rad2 * rad2 - rad1 * rad1) / (distance * rad2);
      float alpha = std::acos(cosAlpha);
      float phi0 = std::atan2(y01 - y02, z01 - z02);
      phi1 = phi0 + alpha;
      phi2 = phi0 - alpha;
    }

    const pandora::CartesianVector& referencePoint1(m_referencePoint);
    const pandora::CartesianVector& referencePoint2(pHelix->GetReferencePoint());

    pandora::CartesianVector position1(0.f, 0.f, 0.f), position2(0.f, 0.f, 0.f);

    if (singlePoint) {
      const float zSect1(z01 + rad1 * std::cos(phi1));
      const float ySect1(y01 + rad1 * std::sin(phi1));
      const float zSect2(z02 + rad2 * std::cos(phi2));
      const float ySect2(y02 + rad2 * std::sin(phi2));
      const float r1(std::sqrt(zSect1 * zSect1 + ySect1 * ySect1));
      const float r2(std::sqrt(zSect2 * zSect2 + ySect2 * ySect2));

      PANDORA_RETURN_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=, this->GetPointOnCircle(r1, referencePoint1, position1));
      PANDORA_RETURN_RESULT_IF(
        pandora::STATUS_CODE_SUCCESS, !=, pHelix->GetPointOnCircle(r2, referencePoint2, position2));
    }
    else {
      const float zSect1(z02 + rad2 * std::cos(phi1));
      const float ySect1(y02 + rad2 * std::sin(phi1));
      const float zSect2(z02 + rad2 * std::cos(phi2));
      const float ySect2(y02 + rad2 * std::sin(phi2));

      const float phiI2(std::atan2(referencePoint2.GetY() - y02, referencePoint2.GetZ() - z02));
      const float phiF21(std::atan2(ySect1 - y02, zSect1 - z02));
      const float phiF22(std::atan2(ySect2 - y02, zSect2 - z02));
      const float charge2(pHelix->GetCharge());

      float deltaPhi21(phiF21 - phiI2);
      float deltaPhi22(phiF22 - phiI2);

      if (deltaPhi21 < 0 && charge2 < 0) { deltaPhi21 += TWO_PI; }
      else if (deltaPhi21 > 0 && charge2 > 0) {
        deltaPhi21 -= TWO_PI;
      }

      if (deltaPhi22 < 0 && charge2 < 0) { deltaPhi22 += TWO_PI; }
      else if (deltaPhi22 > 0 && charge2 > 0) {
        deltaPhi22 -= TWO_PI;
      }

      const float pzy2(pHelix->GetPzy());
      const float genericTime21(-charge2 * deltaPhi21 * rad2 / pzy2);
      const float genericTime22(-charge2 * deltaPhi22 * rad2 / pzy2);

      const float px2(pHelix->GetMomentum().GetX());
      const float X21(referencePoint2.GetX() + genericTime21 * px2);
      const float X22(referencePoint2.GetX() + genericTime22 * px2);

      const pandora::CartesianVector temp21(X21, ySect1, zSect1);
      const pandora::CartesianVector temp22(X22, ySect2, zSect2);

      const float phiI1(std::atan2(referencePoint1.GetY() - y01, referencePoint1.GetZ() - z01));
      const float phiF11(std::atan2(ySect1 - y01, zSect1 - z01));
      const float phiF12(std::atan2(ySect2 - y01, zSect2 - z01));
      const float charge1(m_charge);

      float deltaPhi11(phiF11 - phiI1);
      float deltaPhi12(phiF12 - phiI1);

      if (deltaPhi11 < 0 && charge1 < 0) { deltaPhi11 += TWO_PI; }
      else if (deltaPhi11 > 0 && charge1 > 0) {
        deltaPhi11 -= TWO_PI;
      }

      if (deltaPhi12 < 0 && charge1 < 0) { deltaPhi12 += TWO_PI; }
      else if (deltaPhi12 > 0 && charge1 > 0) {
        deltaPhi12 -= TWO_PI;
      }

      const float pzy1(m_pzy);
      const float genericTime11(-charge1 * deltaPhi11 * rad1 / pzy1);
      const float genericTime12(-charge1 * deltaPhi12 * rad1 / pzy1);

      const float px1(m_momentum.GetX());
      const float X11(referencePoint1.GetX() + genericTime11 * px1);
      const float X12(referencePoint1.GetX() + genericTime12 * px1);

      const pandora::CartesianVector temp11(X11, ySect1, zSect1);
      const pandora::CartesianVector temp12(X12, ySect2, zSect2);

      const float dist1((temp11 - temp21).GetMagnitudeSquared());
      const float dist2((temp12 - temp22).GetMagnitudeSquared());

      if (dist1 < dist2) {
        position1 = temp11;
        position2 = temp21;
      }
      else {
        position1 = temp12;
        position2 = temp22;
      }
    }

    const pandora::CartesianVector momentum1(this->GetExtrapolatedMomentum(position1));
    const pandora::CartesianVector momentum2(pHelix->GetExtrapolatedMomentum(position2));

    helixDistance = (position1 - position2).GetMagnitude();
    positionOfClosestApproach = (position1 + position2) * 0.5;
    v0momentum = momentum1 + momentum2;

    return pandora::STATUS_CODE_SUCCESS;
  }

  //------------------------------------------------------------------------------------------------------------------------------------------

  pandora::CartesianVector Helix::GetExtrapolatedMomentum(
    const pandora::CartesianVector& position) const
  {
    float phi = std::atan2(position.GetY() - m_yCentre, position.GetZ() - m_zCentre);

    if (phi < 0.) phi += TWO_PI;

    phi = phi - m_phiAtPCA + m_phi0;

    return pandora::CartesianVector(
      m_momentum.GetX(), m_pzy * std::sin(phi), m_pzy * std::cos(phi));
  }

} // namespace gar_content
