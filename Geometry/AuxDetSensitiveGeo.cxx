////////////////////////////////////////////////////////////////////////
/// \file  AuxDetSensitiveGeo.cxx
/// \brief Encapsulate the geometry of the sensitive portion of an auxilary detector
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "Geometry/AuxDetSensitiveGeo.h"

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
    AuxDetSensitiveGeo::AuxDetSensitiveGeo(std::vector<const TGeoNode*>& path, int depth)
      : fTotalVolume(0)
    {

      TGeoVolume* vc = path[depth]->GetVolume();
      if (vc) {
        fTotalVolume = vc;
        if (!fTotalVolume)
          throw cet::exception("AuxDetSensitiveGeo") << "cannot find AuxDetSensitive volume\n";

      } // end if found volume

      MF_LOG_DEBUG("Geometry") << "detector sensitive total  volume is " << fTotalVolume->GetName();

      // Build the matrix that takes us to the top world frame
      // build a matrix to take us from the local to the world coordinates
      // in one step
      fGeoMatrix = new TGeoHMatrix(*path[0]->GetMatrix());
      for (int i = 1; i <= depth; ++i) {
        fGeoMatrix->Multiply(path[i]->GetMatrix());
      }

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

    //-----------------------------------------
    AuxDetSensitiveGeo::AuxDetSensitiveGeo(const TGeoVolume* volume, TGeoHMatrix* rotation)
      : fGeoMatrix(rotation), fTotalVolume(volume)
    {

      MF_LOG_DEBUG("Geometry") << "detector sensitive total  volume is " << fTotalVolume->GetName();

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
    AuxDetSensitiveGeo::~AuxDetSensitiveGeo()
    {
      if (fGeoMatrix) delete fGeoMatrix;
      return;
    }

    //......................................................................
    /// Transform a position from local frame to world frame
    /// \param local : 3D array. Position in the local frame  Input.
    /// \param world : 3D array. Position in the world frame. Returned.
    void AuxDetSensitiveGeo::LocalToWorld(const double* local, double* world) const
    {
      fGeoMatrix->LocalToMaster(local, world);
    }

    //......................................................................
    /// Transform a 3-vector from local frame to world frame
    /// \param local : 3D array. Position in the local frame  Input.
    /// \param world : 3D array. Position in the world frame. Returned.
    void AuxDetSensitiveGeo::LocalToWorldVect(const double* local, double* world) const
    {
      fGeoMatrix->LocalToMasterVect(local, world);
    }

    //......................................................................
    /// Transform a position from world frame to local frame
    /// \param world : 3D array. Position in the world frame. Input.
    /// \param local : 3D array. Position in the local frame  Returned.
    void AuxDetSensitiveGeo::WorldToLocal(const double* world, double* local) const
    {
      fGeoMatrix->MasterToLocal(world, local);
    }

    //......................................................................
    /// Transform a 3-vector from world frame to local frame
    /// \param world : 3D array. Position in the world frame. Input.
    /// \param local : 3D array. Position in the local frame  Returned.
    void AuxDetSensitiveGeo::WorldToLocalVect(const double* world, double* local) const
    {
      fGeoMatrix->MasterToLocalVect(world, local);
    }

    //......................................................................
    /// Return the center position of an AuxDet
    /// \param xyz : 3-D array. The returned location.
    /// \param localz : Distance along the length of the volume (z)
    /// (cm). Default is 0.
    void AuxDetSensitiveGeo::GetCenter(double* xyz, double localz) const
    {
      double xyzLocal[3] = {0., 0., localz};
      this->LocalToWorld(xyzLocal, xyz);
    }

    //......................................................................
    // Return the unit normal vector (0,0,1) in local coordinates to global coordinates
    void AuxDetSensitiveGeo::GetNormalVector(double* xyzDir) const
    {
      double normal[3] = {0., 0., 1.};
      this->LocalToWorldVect(normal, xyzDir);
    }

    //......................................................................
    // Get the distance from some point to this detector, xyz in global coordinates
    double AuxDetSensitiveGeo::DistanceToPoint(double* xyz) const
    {
      double Center[3];
      GetCenter(Center);
      return std::sqrt((Center[0] - xyz[0]) * (Center[0] - xyz[0]) +
                       (Center[1] - xyz[1]) * (Center[1] - xyz[1]) +
                       (Center[2] - xyz[2]) * (Center[2] - xyz[2]));
    }
  }
} // gar
////////////////////////////////////////////////////////////////////////
