#ifndef ROTATIONTRANSFORMATION_H
#define ROTATIONTRANSFORMATION_H 1

#include "Api/PandoraApi.h"

#include "CLHEP/Units/PhysicalConstants.h"

namespace gar {
  namespace gar_pandora {

    class RotationTransformation {
    public:
      RotationTransformation(int kAxis, float angle);

      virtual ~RotationTransformation();

      const pandora::CartesianVector MakeRotation(const pandora::CartesianVector& initialVec) const;

      static const int kAxisX;
      static const int kAxisY;
      static const int kAxisZ;

    private:
      void SetRotationX();
      void SetRotationY();
      void SetRotationZ();

      bool m_RotationSet;
      float m_RotationAngle;
      float fRotMatrix[9]; //the rotation matrix
    };
  }
}

#endif
