/**
* @file   larcorealg/Geometry/LocalTransformation.tcc
* @brief  Class containing local-to-world transformations
*         (template implementation)
* @author Gianluca Petrillo (petrillo@fnal.gov)
* @date   November 30, 2016
* @see    LocalTransformation.h
*
* This file is expected to be included directly in the header.
*
*/

#ifndef GARSOFT_GEOMETRY_LOCALTRANSFORMATION_TCC
#define GARSOFT_GEOMETRY_LOCALTRANSFORMATION_TCC

// framework libraries
#include "cetlib_except/exception.h"
// ROOT
#include "TGeoNode.h"
#include "TGeoMatrix.h"
// CLHEP
#include "CLHEP/Geometry/Transform3D.h"
// C standard library
#include <cassert>

//------------------------------------------------------------------------------

namespace gar {
    namespace geo {
        namespace details {

            template <typename T, std::size_t SrcN = 3, std::size_t DestN = SrcN>
            bool doBuffersOverlap(T const* src, T const* dest)
            { return (dest < (src + SrcN)) && (src < (dest + DestN)); }

            template <typename T, std::size_t SrcN = 3, std::size_t DestN = SrcN>
            void checkVectorBufferOverlap(T const* src, T const* dest) {
                if (doBuffersOverlap<T, SrcN, DestN>(src, dest)) {
                    throw cet::exception("LocalTransformation")
                    << "source " << SrcN << "@[" << ((void*) src) << "]"
                    << " and destination " << DestN << "@[" << ((void*) dest) << "]"
                    << " buffers overlap!\n";
                }
                assert(!doBuffersOverlap(src, dest));
            } // checkVectorBufferOverlap()

        } // namespace details
    } // namespace geo
} //namespace gar

//------------------------------------------------------------------------------

template <typename Matrix>
void gar::geo::LocalTransformation<Matrix>::LocalToWorld
(double const* local, double* world) const
{
    details::checkVectorBufferOverlap(local, world);
    fGeoMatrix.LocalToMaster(local, world);
} // geo::LocalTransformation::LocalToWorld()

//------------------------------------------------------------------------------

template <typename Matrix>
void gar::geo::LocalTransformation<Matrix>::LocalToWorldVect
(double const* local, double* world) const
{
    details::checkVectorBufferOverlap(local, world);
    fGeoMatrix.LocalToMasterVect(local, world);
} // geo::LocalTransformation::LocalToWorldVect()

//------------------------------------------------------------------------------

template <typename Matrix>
void gar::geo::LocalTransformation<Matrix>::WorldToLocal
(double const* world, double* local) const
{
    details::checkVectorBufferOverlap(local, world);
    fGeoMatrix.MasterToLocal(world, local);
} // geo::LocalTransformation::WorldToLocal()

//------------------------------------------------------------------------------

template <typename Matrix>
void gar::geo::LocalTransformation<Matrix>::WorldToLocalVect
(const double* world, double* local) const
{
    details::checkVectorBufferOverlap(local, world);
    fGeoMatrix.MasterToLocalVect(world, local);
} // geo::LocalTransformation::WorldToLocalVect()

//------------------------------------------------------------------------------
// specialisations (forward declarations)
//
namespace gar {
    namespace geo {

        template <>
        TGeoHMatrix LocalTransformation<TGeoHMatrix>::transformationFromPath
        (std::vector<TGeoNode const*> const& path, size_t depth);

        template <>
        HepGeom::Transform3D
        LocalTransformation<HepGeom::Transform3D>::transformationFromPath
        (std::vector<TGeoNode const*> const& path, size_t depth);

    } // namespace geo
} // namespace gar

//------------------------------------------------------------------------------

#endif // GARSOFT_GEOMETRY_LOCALTRANSFORMATION_TCC
