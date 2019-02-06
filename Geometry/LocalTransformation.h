/**
* @file   larcorealg/Geometry/LocalTransformation.h
* @brief  Class containing local-to-world transformations
* @author Gianluca Petrillo (petrillo@fnal.gov)
* @date   November 30, 2016
* @ingroup Geometry
*
*/

#ifndef GARSOFT_GEOMETRY_LOCALTRANSFORMATION_H
#define GARSOFT_GEOMETRY_LOCALTRANSFORMATION_H

// ROOT libraries
// (none)
// C/C++ standard libraries
#include <vector>
#include <utility> // std::move()
#include <cstdlib> // std::size_t

// forward declarations
class TGeoNode;

namespace gar {
    namespace geo {

        /**
        * @brief Class to transform between world and local coordinates
        * @tparam StoredMatrix type of transformation matrix internally stored
        * @ingroup Geometry
        *
        * This class provides two directions of transformations (world to local and
        * the other way around), for points and for vectors.
        * The vector version of the transformation does not apply translation.
        *
        * @note In the class method examples, the following definition is assumed:
        * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
        * using LocalTransformation_t = geo::LocalTransformation<TGeoHMatrix>;
        * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        template <typename StoredMatrix>
        class LocalTransformation {
        public:

            /// Type of transformation matrix
            using TransformationMatrix_t = StoredMatrix;

            //@{
            /**
            * @brief Constructor: uses the specified transformation matrix
            * @param matrix the transformation matrix to be used
            *
            * The specified matrix is copied into a local copy.
            */
            LocalTransformation(TransformationMatrix_t const& matrix)
            : fGeoMatrix(matrix) {}
            LocalTransformation(TransformationMatrix_t&& matrix)
            : fGeoMatrix(std::move(matrix)) {}
            //@}

            /**
            * @brief Constructor: uses the specified transformation matrix
            * @param path the path of ROOT geometry nodes
            * @param depth the index in the path of the last node to be considered
            *
            * The specified matrix is copied into a local copy.
            */
            LocalTransformation(std::vector<TGeoNode const*> const& path, size_t depth)
            : fGeoMatrix(transformationFromPath(path, depth)) {}

            /**
            * @brief Transforms a point from local frame to world frame
            * @param local local coordinates: [0] x, [1] y, [2] z [cm]
            * @param world (output) corresponding world coordinates [cm]
            *
            * The full transformation is applied. Fox example:
            * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            * LocalTransformation_t trans( ... ); // with proper initialisation
            *
            * std::array<double, 3U> origin, center;
            * origin.fill(0.);
            * trans.LocalToWorld(origin.data(), center.data());
            * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            * `center` will contain the world coordinates of the center of the volume,
            * which is usually represented by the origin in the local coordinates.
            *
            * In-place replacement is *not* supported: `world` and `local` buffers are
            * assumed not to, and must not, overlap.
            */
            void LocalToWorld(double const* local, double* world) const;

            /**
            * @brief Transforms a vector from local frame to world frame
            * @param local local coordinates: [0] x, [1] y, [2] z [cm]
            * @param world (output) corresponding world coordinates [cm]
            *
            * The translation is not applied, since the argument is supposed to be a
            * vector, relative difference between two points.
            *
            * In-place replacement is *not* supported: `world` and `local` buffers are
            * assumed not to, and must not, overlap.
            */
            void LocalToWorldVect(double const* local, double* world) const;

            /**
            * @brief Transforms a point from world frame to local frame
            * @param world world coordinates: [0] x, [1] y, [2] z [cm]
            * @param local (output) corresponding local coordinates [cm]
            *
            * The full transformation is applied. Fox example:
            * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            * LocalTransformation_t trans( ... ); // with proper initialisation
            *
            * std::array<double, 3U> world{ 4.0, 5.0, -2.5 }, local;
            * trans.WorldToLocal(world.data(), local.data());
            * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            * `local` will contain the local coordinates of the specified point.
            *
            * In-place replacement is *not* supported: `world` and `local` buffers are
            * assumed not to, and must not, overlap.
            */
            void WorldToLocal(double const* world, double* local) const;

            /**
            * @brief Transforms a vector from world frame to local frame
            * @param world world coordinates: [0] x, [1] y, [2] z [cm]
            * @param local (output) corresponding local coordinates [cm]
            *
            * The translation is not applied, since the argument is supposed to be a
            * vector, relative difference between two points.
            *
            * In-place replacement is *not* supported: `world` and `local` buffers are
            * assumed not to, and must not, overlap.
            */
            void WorldToLocalVect(const double* world, double* local) const;


            /// Direct access to the transformation matrix
            TransformationMatrix_t const& Matrix() const { return fGeoMatrix; }


            /// Builds a matrix to go from local to world coordinates in one step
            static TransformationMatrix_t transformationFromPath
            (std::vector<TGeoNode const*> const& path, size_t depth);


        protected:

            TransformationMatrix_t fGeoMatrix; ///< local to world transform

        }; // class LocalTransformation<>


    } // namespace geo
} //namespace gar

//------------------------------------------------------------------------------
// template implementation
#include "LocalTransformation.tcc"
//------------------------------------------------------------------------------
#endif // GARSOFT_GEOMETRY_LOCALTRANSFORMATION_H
