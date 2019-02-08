/**
*  @file   LCContent/src/LCUtility/KDTreeLinkerToolsT.cc
*
*  @brief  Implementation of the kd tree linker tools template class
*
*  $Log: $
*/

#include "Utilities/KDTreeLinkerToolsT.h"

namespace util
{

    std::pair<float, float> minmax(const float a, const float b)
    {
        return ((b < a) ? std::pair<float, float>(b, a) : std::pair<float, float>(a, b));
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    KDTreeTesseract fill_and_bound_4d_kd_tree(const std::list<const gar::rec::InternalCaloHit *> &points, std::vector<KDTreeNodeInfoT<const gar::rec::InternalCaloHit*, 4> > &nodes)
    {
        std::array<double, 4> minpos{ {0.f, 0.f, 0.f, 0.f} }, maxpos{ {0.f, 0.f, 0.f, 0.f} };

        unsigned i = 0;

        for (const gar::rec::InternalCaloHit *const point : points)
        {
            const CLHEP::Hep3Vector &pos = kdtree_type_adaptor<const gar::rec::InternalCaloHit>::position(point);
            const double layer = static_cast<double>(point->GetLayer());
            nodes.emplace_back(point, pos.x(), pos.y(), pos.z(), layer);

            if (0 == i)
            {
                minpos[0] = pos.x(); minpos[1] = pos.y(); minpos[2] = pos.z(); minpos[3] = layer;
                maxpos[0] = pos.x(); maxpos[1] = pos.y(); maxpos[2] = pos.z(); maxpos[3] = layer;
            }
            else
            {
                minpos[0] = std::min(pos.x(), minpos[0]);
                minpos[1] = std::min(pos.y(), minpos[1]);
                minpos[2] = std::min(pos.z(), minpos[2]);
                minpos[3] = std::min(layer, minpos[3]);
                maxpos[0] = std::max(pos.x(), maxpos[0]);
                maxpos[1] = std::max(pos.y(), maxpos[1]);
                maxpos[2] = std::max(pos.z(), maxpos[2]);
                maxpos[3] = std::max(layer, maxpos[3]);
            }

            ++i;
        }

        return KDTreeTesseract(minpos[0], maxpos[0], minpos[1], maxpos[1], minpos[2], maxpos[2], minpos[3], maxpos[3]);
    }


    //------------------------------------------------------------------------------------------------------------------------------------------

    KDTreeCube build_3d_kd_search_region(const gar::rec::InternalCaloHit *const point, const float x_span, const float y_span, const float z_span)
    {
        const CLHEP::Hep3Vector &pos = point->GetPositionVector();

        const auto x_side = minmax(pos.x() + x_span, pos.x() - x_span);
        const auto y_side = minmax(pos.y() + y_span, pos.y() - y_span);
        const auto z_side = minmax(pos.z() + z_span, pos.z() - z_span);

        return KDTreeCube(x_side.first, x_side.second, y_side.first, y_side.second, z_side.first, z_side.second);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    KDTreeTesseract build_4d_kd_search_region(const gar::rec::InternalCaloHit *const point, const float x_span, const float y_span, const float z_span,
    const float search_layer)
    {
        return build_4d_kd_search_region(point->GetPositionVector(), x_span, y_span, z_span, search_layer);
    }

    //------------------------------------------------------------------------------------------------------------------------------------------

    KDTreeTesseract build_4d_kd_search_region(const CLHEP::Hep3Vector &pos, const float x_span, const float y_span, const float z_span,
    const float search_layer)
    {
        const auto x_side = minmax(pos.x() + x_span, pos.x() - x_span);
        const auto y_side = minmax(pos.y() + y_span, pos.y() - y_span);
        const auto z_side = minmax(pos.z() + z_span, pos.z() - z_span);
        const auto layer_side = minmax(search_layer + 0.5f, search_layer - 0.5f);

        return KDTreeTesseract(x_side.first, x_side.second, y_side.first, y_side.second, z_side.first, z_side.second, layer_side.first, layer_side.second);
    }

} // namespace lc_content
