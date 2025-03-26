////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeo.cxx
/// \brief Encapsulate the geometry of an auxilary detector
///
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <limits>

#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetGeoObjectSorter.h"

#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoNode.h"
#include "TGeoTrd1.h"
#include "TGeoTrd2.h"
#include "TGeoTube.h"
#include "TMath.h"
#include "TVector3.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace geo {

    //-----------------------------------------
    AuxDetGeo::AuxDetGeo(std::vector<const TGeoNode*>& path, int depth) : fTotalVolume(0)
    {

      TGeoVolume* vc = path[depth]->GetVolume();
      if (vc) {
        fTotalVolume = vc;
        if (!fTotalVolume) throw cet::exception("AuxDetGeo") << "cannot find AuxDet volume\n";

      } // end if found volume

      MF_LOG_DEBUG("Geometry") << "detector total  volume is " << fTotalVolume->GetName();

      // Build the matrix that takes us to the top world frame
      // build a matrix to take us from the local to the world coordinates
      // in one step
      fGeoMatrix = new TGeoHMatrix(*path[0]->GetMatrix());
      for (int i = 1; i <= depth; ++i) {
        fGeoMatrix->Multiply(path[i]->GetMatrix());
      }

      // look for sensitive volumes - if there are none then this aux det
      // could be from an older gdml file than the introduction of AuxDetSensitiveGeo
      // in that case assume the full AuxDetGeo is sensitive and copy its information
      // into a single AuxDetSensitive
      this->FindAuxDetSensitive(path, depth);
      if (fSensitive.size() < 1)
        fSensitive.push_back(new AuxDetSensitiveGeo(fTotalVolume, fGeoMatrix));

      // set the ends depending on whether the shape is a box or trapezoid
      std::string volName(fTotalVolume->GetName());
      if (volName.find("Trap") != std::string::npos) {

        //       Small Width
        //          ____          Height is the thickness
        //         /    \     T     of the trapezoid
        //        /      \    |
        //       /        \   | Length
        //      /__________\  _
        //         Width
        fHalfHeight = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDy1(); // same as Dy2()
        fLength = 2.0 * ((TGeoTrd2*)fTotalVolume->GetShape())->GetDz();
        fHalfWidth1 = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx1(); // at -Dz
        fHalfWidth2 = ((TGeoTrd2*)fTotalVolume->GetShape())->GetDx2(); // at +Dz
      }
      else {
        fHalfWidth1 = ((TGeoBBox*)fTotalVolume->GetShape())->GetDX();
        fHalfHeight = ((TGeoBBox*)fTotalVolume->GetShape())->GetDY();
        fLength = 2.0 * ((TGeoBBox*)fTotalVolume->GetShape())->GetDZ();
        fHalfWidth2 = fHalfWidth1;
      }
    }

    //......................................................................
    AuxDetGeo::~AuxDetGeo()
    {
      if (fGeoMatrix) delete fGeoMatrix;
      return;
    }

    //......................................................................
    void AuxDetGeo::FindAuxDetSensitive(std::vector<const TGeoNode*>& path, unsigned int depth)
    {
      // Check if the current node is a senstive volume - we already know
      // we are in an Auxiliary detector
      std::string pathName(path[depth]->GetName());
      if (pathName.find("Sensitive") != std::string::npos) {
        this->MakeAuxDetSensitive(path, depth);
        return;
      }

      // Explore the next layer down
      unsigned int deeper = depth + 1;
      if (deeper >= path.size()) {
        throw cet::exception("ExceededMaxDepth") << "Exceeded maximum depth\n";
      }
      const TGeoVolume* v = path[depth]->GetVolume();
      int nd = v->GetNdaughters();
      for (int i = 0; i < nd; ++i) {
        path[deeper] = v->GetNode(i);
        this->FindAuxDetSensitive(path, deeper);
      }
    }

    //......................................................................
    void AuxDetGeo::MakeAuxDetSensitive(std::vector<const TGeoNode*>& path, int depth)
    {
      fSensitive.push_back(new AuxDetSensitiveGeo(path, depth));
    }

    //......................................................................
    /// Transform a position from local frame to world frame
    /// \param local : 3D array. Position in the local frame  Input.
    /// \param world : 3D array. Position in the world frame. Returned.
    void AuxDetGeo::LocalToWorld(const double* local, double* world) const
    {
      fGeoMatrix->LocalToMaster(local, world);
    }

    //......................................................................
    /// Transform a 3-vector from local frame to world frame
    /// \param local : 3D array. Position in the local frame  Input.
    /// \param world : 3D array. Position in the world frame. Returned.
    void AuxDetGeo::LocalToWorldVect(const double* local, double* world) const
    {
      fGeoMatrix->LocalToMasterVect(local, world);
    }

    //......................................................................
    /// Transform a position from world frame to local frame
    /// \param world : 3D array. Position in the world frame. Input.
    /// \param local : 3D array. Position in the local frame  Returned.
    void AuxDetGeo::WorldToLocal(const double* world, double* local) const
    {
      fGeoMatrix->MasterToLocal(world, local);
    }

    //......................................................................
    /// Transform a 3-vector from world frame to local frame
    /// \param world : 3D array. Position in the world frame. Input.
    /// \param local : 3D array. Position in the local frame  Returned.
    void AuxDetGeo::WorldToLocalVect(const double* world, double* local) const
    {
      fGeoMatrix->MasterToLocalVect(world, local);
    }

    //......................................................................
    /// Return the center position of an AuxDet
    /// \param xyz : 3-D array. The returned location.
    /// \param localz : Distance along the length of the volume (z)
    /// (cm). Default is 0.
    void AuxDetGeo::GetCenter(double* xyz, double localz) const
    {
      double xyzLocal[3] = {0., 0., localz};
      this->LocalToWorld(xyzLocal, xyz);
    }

    //......................................................................
    // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
    void AuxDetGeo::GetNormalVector(double* xyzDir) const
    {
      double normal[3] = {0., 0., 1.};
      this->LocalToWorldVect(normal, xyzDir);
    }

    //......................................................................
    // Get the distance from some point to this detector, xyz in global coordinates
    double AuxDetGeo::DistanceToPoint(double* xyz) const
    {
      double Center[3];
      GetCenter(Center);
      return std::sqrt((Center[0] - xyz[0]) * (Center[0] - xyz[0]) +
                       (Center[1] - xyz[1]) * (Center[1] - xyz[1]) +
                       (Center[2] - xyz[2]) * (Center[2] - xyz[2]));
    }

    //......................................................................
    size_t AuxDetGeo::FindSensitiveVolume(double const worldPos[3]) const
    {
      double local[3] = {0.};
      for (unsigned int a = 0; a < fSensitive.size(); ++a) {

        this->SensitiveVolume(a).WorldToLocal(worldPos, local);
        double HalfCenterWidth =
          (this->SensitiveVolume(a).HalfWidth1() + this->SensitiveVolume(a).HalfWidth2()) / 2;

        if (local[2] >= -this->SensitiveVolume(a).Length() / 2 &&
            local[2] <= this->SensitiveVolume(a).Length() / 2 &&
            local[1] >= -this->SensitiveVolume(a).HalfHeight() &&
            local[1] <= this->SensitiveVolume(a).HalfHeight() &&
            // if SensitiveVolume a is a box, then HalfSmallWidth = HalfWidth
            local[0] >= -HalfCenterWidth +
                          local[2] * (HalfCenterWidth - this->SensitiveVolume(a).HalfWidth2()) /
                            (this->SensitiveVolume(a).Length() / 2) &&
            local[0] <= HalfCenterWidth -
                          local[2] * (HalfCenterWidth - this->SensitiveVolume(a).HalfWidth2()) /
                            (this->SensitiveVolume(a).Length() / 2))
          return a;

      } // for loop over AuxDetSensitive a

      throw cet::exception("Geometry") << "Can't find AuxDetSensitive for position (" << worldPos[0]
                                       << "," << worldPos[1] << "," << worldPos[2] << ")\n";

      return UINT_MAX;
    }

    //......................................................................
    AuxDetSensitiveGeo const& AuxDetGeo::PositionToSensitiveVolume(double const worldLoc[3],
                                                                   size_t& sv) const
    {
      sv = this->FindSensitiveVolume(worldLoc);
      return this->SensitiveVolume(sv);
    }

    //......................................................................
    void AuxDetGeo::SortSubVolumes(AuxDetGeoObjectSorter const& sorter)
    {
      sorter.SortAuxDetSensitive(fSensitive);

      return;
    }

  }
} // gar
////////////////////////////////////////////////////////////////////////
