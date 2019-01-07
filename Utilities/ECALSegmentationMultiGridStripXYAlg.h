#ifndef ECALSEGMENTATIONMULTIGRIDSTRIPXYALG_H
#define ECALSEGMENTATIONMULTIGRIDSTRIPXYALG_H

#include "Utilities/ECALSegmentationAlg.h"

#include <string>
#include <vector>

#include "TGeoShape.h"

namespace fhicl{
    class ParameterSet;
}

#define CM_2_MM 10

namespace util {

    class ECALSegmentationMultiGridStripXYAlg: public ECALSegmentationAlg {

    public:
        ECALSegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset);

        ECALSegmentationMultiGridStripXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

        virtual ~ECALSegmentationMultiGridStripXYAlg();

        virtual void reconfigure(fhicl::ParameterSet const& pset);

        virtual ROOT::Math::XYZVector position(const long64& cID) const;

        virtual long64 cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const;

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

        void setLayerDimXY(const double dimX, const double dimY) { _layer_dim_X = dimX; _layer_dim_Y = dimY; }

    protected:

        virtual void PrintParameters() const;

        void TokenizeLayerVectors(std::vector<std::string> grid, std::vector<std::string> strip) const;

        void CheckLayerConfiguration(const unsigned int layer) const;

        void SetLayerParameters(TGeoShape *shape);

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
        /// the layers (start::end) for grid
        std::vector<std::string> _gridLayers;
        /// the layers (start::end) for strips
        std::vector<std::string> _stripLayers;

        mutable std::vector<unsigned int> _gridFirst;
        mutable std::vector<unsigned int> _gridLast;

        mutable std::vector<unsigned int> _stripFirst;
        mutable std::vector<unsigned int> _stripLast;

        mutable bool _isTile;

        double _layer_dim_X;

        double _layer_dim_Y;
    };
}

#endif
