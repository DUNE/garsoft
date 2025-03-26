////////////////////////////////////////////////////////////////////////
/// \file  AuxDetSensitiveGeo.h
/// \brief Encapsulate the geometry of the sensitive portion of an auxiliary detector
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GEO_AuxDetSensitiveGeo_H
#define GEO_AuxDetSensitiveGeo_H
#include "TGeoVolume.h"
#include <vector>

class TGeoNode;
class TGeoHMatrix;
class TGeoMatrix;

namespace gar {
  namespace geo {

    class AuxDetSensitiveGeo {
    public:
      AuxDetSensitiveGeo(std::vector<const TGeoNode*>& path, int depth);
      AuxDetSensitiveGeo(const TGeoVolume* volume, TGeoHMatrix* rotation);
      ~AuxDetSensitiveGeo();

      void GetCenter(double* xyz, double localz = 0.0) const;
      void GetNormalVector(double* xyzDir) const;

      //box geometry
      double Length() const { return fLength; }
      double HalfWidth1() const { return fHalfWidth1; }
      double HalfWidth2() const { return fHalfWidth2; }
      double HalfHeight() const { return fHalfHeight; }
      const TGeoVolume* TotalVolume() const { return fTotalVolume; }

      double DistanceToPoint(double* xyz) const;

      void LocalToWorld(const double* local, double* world) const;
      void LocalToWorldVect(const double* local, double* world) const;
      void WorldToLocal(const double* world, double* local) const;
      void WorldToLocalVect(const double* world, double* local) const;

    private:
      TGeoHMatrix* fGeoMatrix;        ///< Transformation matrix to world frame
      const TGeoVolume* fTotalVolume; ///< Total volume of AuxDet, called vol*
      double fLength;                 ///< length of volume, along z direction in local
      double fHalfWidth1;             ///< 1st half width of volume, at -z/2 in local coordinates
      double fHalfWidth2;             ///< 2nd half width (width1==width2 for boxes), at +z/2
      double fHalfHeight;             ///< half height of volume
    };
  }

} // gar
#endif
