#ifndef ECALSEGMENTATIONMULTIGRIDSTRIPXYALG_H
#define ECALSEGMENTATIONMULTIGRIDSTRIPXYALG_H

#include "Geometry/ECALSegmentationAlg.h"

#include <string>
#include <vector>

namespace fhicl{
    class ParameterSet;
}

#define CM_2_MM 10

namespace gar {
    namespace geo {

        class ECALSegmentationMultiGridStripXYAlg: public ECALSegmentationAlg {

        public:
            ECALSegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset);

            ECALSegmentationMultiGridStripXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

            virtual ~ECALSegmentationMultiGridStripXYAlg();

            virtual void reconfigure(fhicl::ParameterSet const& pset);

            virtual ROOT::Math::XYZVector position(const gar::geo::GeometryCore& geo, const long64& cID) const;

            virtual long64 cellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const ROOT::Math::XYZVector& localPosition) const;

            virtual int getIDbyCellID(const long64& cID, const char* id) const;

            virtual void Initialize(const gar::geo::GeometryCore & geo);

            double gridSizeX() const { return _gridSizeX; }

            double gridSizeY() const { return _gridSizeY; }

            double stripSizeX() const { return _stripSizeX; }

            double stripSizeY() const { return _stripSizeY; }

            double offsetX() const { return _offsetX; }

            double offsetY() const { return _offsetY; }

            double layerDimX() const { return _layer_dim_X; }

            double layerDimY() const { return _layer_dim_Y; }

            const std::string& fieldNameX() const { return _xId; }

            const std::string& fieldNameY() const { return _yId; }

            const std::string& fieldNameLayer() const { return _layerId; }

            const std::string& fieldNameSlice() const { return _sliceId; }

            void setGridSizeX(double cellSize) { _gridSizeX = cellSize; }

            void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }

            void setStripSizeX(double stripSize) { _stripSizeX = stripSize; }

            void setStripSizeY(double stripSize) { _stripSizeY = stripSize; }

            void setOffsetX(double offset) { _offsetX = offset; }

            void setOffsetY(double offset) { _offsetY = offset; }

            void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

            void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

            void setFieldNameLayer(const std::string& fieldName) { _layerId = fieldName; }

            void setFieldNameSlice(const std::string& fieldName) { _sliceId = fieldName; }

            void setLayerDimXY(const double& dimX, const double& dimY) const { _layer_dim_X = dimX; _layer_dim_Y = dimY; }

        protected:

            virtual void PrintParameters() const;

            void TokenizeLayerVectors(std::vector<std::string> grid, std::vector<std::string> strip) const;

            void CheckLayerConfiguration(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer) const;

            bool GetSegmentationConfig(const unsigned int& det_id, const unsigned int& layer) const;

            void UpdateGeometrySegmentation(const gar::geo::GeometryCore& geo) const;

            /// the grid size in X
            double _gridSizeX;
            /// the coordinate offset in X
            double _offsetX;
            /// the grid size in Y
            double _gridSizeY;
            /// the coordinate offset in Y
            double _offsetY;
            /// the field name used for X
            std::string _xId;
            /// the field name used for Y
            std::string _yId;
            /// the field name used for layer
            std::string _layerId;
            /// the field name used for slice
            std::string _sliceId;
            /// the strip size in X
            double _stripSizeX;
            /// the strip size in Y
            double _stripSizeY;

            /// the layers (start::end) for grid in Barrel
            std::vector<std::string> _gridBarrelLayers;
            /// the layers (start::end) for strips in Barrel
            std::vector<std::string> _stripBarrelLayers;
            /// the layers (start::end) for grid in Endcap
            std::vector<std::string> _gridEndcapLayers;
            /// the layers (start::end) for strips in Endcap
            std::vector<std::string> _stripEndcapLayers;

            /// vectors containing the layers that are tiled or stripped
            mutable std::vector<unsigned int> _gridFirst;
            mutable std::vector<unsigned int> _gridLast;
            mutable std::vector<unsigned int> _stripFirst;
            mutable std::vector<unsigned int> _stripLast;

            /// mutable boolean to check if the layer is composed of tiles or strips
            mutable bool _isTile;
            /// layer dimension in X
            mutable double _layer_dim_X;
            /// layer dimension in Y
            mutable double _layer_dim_Y;
        };
    }
}

#endif
