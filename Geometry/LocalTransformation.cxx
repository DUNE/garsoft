/**
* @file   larcorealg/Geometry/LocalTransformation.cxx
* @brief  Class containing local-to-world transformations (implementation file)
* @author Gianluca Petrillo (petrillo@fnal.gov)
* @date   November 30, 2016
* @see    LocalTransformation.h
*
*/
// our header
#include "LocalTransformation.h"
// ROOT
#include "TGeoMatrix.h" // TGeoHMatrix
// CLHEP
#include "CLHEP/Geometry/Transform3D.h" // HepGeom::Transform3D
#include "CLHEP/Vector/Rotation.h" // CLHEP::HepRotation
#include "CLHEP/Vector/RotationInterfaces.h" // CLHEP::HepRep3x3
#include "CLHEP/Vector/ThreeVector.h" // CLHEP::Hep3Vector
//------------------------------------------------------------------------------
// specialisations
//
namespace gar {
    namespace geo {

        template <>
        TGeoHMatrix LocalTransformation<TGeoHMatrix>::transformationFromPath
        (std::vector<TGeoNode const*> const& path, size_t depth)
        {

            TGeoHMatrix matrix = *(path[0]->GetMatrix());
            for(size_t i = 1; i <= depth; ++i)
            matrix.Multiply(path[i]->GetMatrix());
            return matrix;

        } // gar::geo::LocalTransformation::transformationFromPath()

        //----------------------------------------------------------------------------
        template <>
        HepGeom::Transform3D
        LocalTransformation<HepGeom::Transform3D>::transformationFromPath
        (std::vector<TGeoNode const*> const& path, size_t depth)
        {

            auto const mat =
            LocalTransformation<TGeoHMatrix>::transformationFromPath
            (path, depth);
            const Double_t* translation = mat.GetTranslation();
            return HepGeom::Transform3D(
            CLHEP::HepRotation(CLHEP::HepRep3x3(mat.GetRotationMatrix())),
            CLHEP::Hep3Vector(translation[0], translation[1], translation[2])
            );

        } // gar::geo::LocalTransformation::transformationFromPath()
    } // namespace geo
} // namespace gar

//------------------------------------------------------------------------------
