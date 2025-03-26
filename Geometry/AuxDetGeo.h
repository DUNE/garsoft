////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeo.h
/// \brief Encapsulate the geometry of an auxiliary detector
///
/// \author  miceli@fnal.gov
////////////////////////////////////////////////////////////////////////

#ifndef GEO_AuxDetGeo_H
#define GEO_AuxDetGeo_H
#include "TGeoVolume.h"
#include <vector>

#include "Geometry/AuxDetSensitiveGeo.h"

class TGeoNode;
class TGeoHMatrix;
class TGeoMatrix;

namespace gar {
  namespace geo {

    class AuxDetGeoObjectSorter;

    class AuxDetGeo {
    public:
      AuxDetGeo(std::vector<const TGeoNode*>& path, int depth);
      ~AuxDetGeo();

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

      std::string Name() const { return fTotalVolume->GetName(); }

      // methods for the sensitive volumes in the aux det

      size_t FindSensitiveVolume(double const worldLoc[3]) const;
      AuxDetSensitiveGeo const& PositionToSensitiveVolume(double const worldLoc[3],
                                                          size_t& sv) const;
      AuxDetSensitiveGeo const& SensitiveVolume(size_t sv) const { return *fSensitive[sv]; }
      size_t NSensitiveVolume() const { return fSensitive.size(); }

      void SortSubVolumes(geo::AuxDetGeoObjectSorter const& sorter);

    private:
      void FindAuxDetSensitive(std::vector<const TGeoNode*>& path, unsigned int depth);
      void MakeAuxDetSensitive(std::vector<const TGeoNode*>& path, int depth);

      const TGeoVolume* fTotalVolume; ///< Total volume of AuxDet, called vol*
      TGeoHMatrix* fGeoMatrix;        ///< Transformation matrix to world frame
      double fLength;                 ///< length of volume, along z direction in local
      double fHalfWidth1;             ///< 1st half width of volume, at -z/2 in local coordinates
      double fHalfWidth2;             ///< 2nd half width (width1==width2 for boxes), at +z/2
      double fHalfHeight;             ///< half height of volume
      std::vector<AuxDetSensitiveGeo*> fSensitive; ///< sensitive volumes in the detector
    };
  }

} // gar

#endif
