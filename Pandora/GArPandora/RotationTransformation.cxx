#include "RotationTransformation.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace gar_pandora {

    const int RotationTransformation::kAxisX = 0;
    const int RotationTransformation::kAxisY = 1;
    const int RotationTransformation::kAxisZ = 2;

    //------------------------------------------------------------------------------------------------------------------------------------------

    RotationTransformation::RotationTransformation(int axis, float angle)
      : m_RotationSet(false), m_RotationAngle(angle * CLHEP::pi / 180.)
    {
      switch (axis) {
      case kAxisX:
        this->SetRotationX();
        MF_LOG_INFO("RotationTransformation")
          << "Rotation set around the x axis with angle " << m_RotationAngle << " rad";
        m_RotationSet = true;
        break;
      case kAxisY:
        this->SetRotationY();
        MF_LOG_INFO("RotationTransformation")
          << "Rotation set around the y axis with angle " << m_RotationAngle << " rad";
        m_RotationSet = true;
        break;
      case kAxisZ:
        this->SetRotationZ();
        MF_LOG_INFO("RotationTransformation")
          << "Rotation set around the z axis with angle " << m_RotationAngle << " rad";
        m_RotationSet = true;
        break;
      default: throw; break;
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    RotationTransformation::~RotationTransformation() {}

    //------------------------------------------------------------------------------------------------------------------------------------------

    const pandora::CartesianVector RotationTransformation::MakeRotation(
      const pandora::CartesianVector& initialVec) const
    {
      if (m_RotationSet) {
        return pandora::CartesianVector(
          initialVec.GetX() * fRotMatrix[0] + initialVec.GetY() * fRotMatrix[3] +
            initialVec.GetZ() * fRotMatrix[6],
          initialVec.GetX() * fRotMatrix[1] + initialVec.GetY() * fRotMatrix[4] +
            initialVec.GetZ() * fRotMatrix[7],
          initialVec.GetX() * fRotMatrix[2] + initialVec.GetY() * fRotMatrix[5] +
            initialVec.GetZ() * fRotMatrix[8]);
      }
      else {
        MF_LOG_INFO("RotationTransformation")
          << "No rotation axis was defined - operation aborted!";
        return pandora::CartesianVector(initialVec.GetX(), initialVec.GetY(), initialVec.GetZ());
      }
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void RotationTransformation::SetRotationX()
    {
      //col x
      fRotMatrix[0] = 1.;
      fRotMatrix[1] = 0.;
      fRotMatrix[2] = 0.;

      //col y
      fRotMatrix[3] = 0.;
      fRotMatrix[4] = std::cos(m_RotationAngle);
      fRotMatrix[5] = std::sin(m_RotationAngle);

      //col z
      fRotMatrix[6] = 0.;
      fRotMatrix[7] = -std::sin(m_RotationAngle);
      fRotMatrix[8] = std::cos(m_RotationAngle);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void RotationTransformation::SetRotationY()
    {
      //col x
      fRotMatrix[0] = std::cos(m_RotationAngle);
      fRotMatrix[1] = 0.;
      fRotMatrix[2] = -std::sin(m_RotationAngle);

      //col y
      fRotMatrix[3] = 0.;
      fRotMatrix[4] = 1.;
      fRotMatrix[5] = 0.;

      //col z
      fRotMatrix[6] = std::sin(m_RotationAngle);
      fRotMatrix[7] = 0.;
      fRotMatrix[8] = std::cos(m_RotationAngle);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    void RotationTransformation::SetRotationZ()
    {
      //col x
      fRotMatrix[0] = std::cos(m_RotationAngle);
      fRotMatrix[1] = std::sin(m_RotationAngle);
      fRotMatrix[2] = 0.;

      //col y
      fRotMatrix[3] = -std::sin(m_RotationAngle);
      fRotMatrix[4] = std::cos(m_RotationAngle);
      fRotMatrix[5] = 0.;

      //col z
      fRotMatrix[6] = 0.;
      fRotMatrix[7] = 0.;
      fRotMatrix[8] = 1.;
    }
  }
}
