#ifndef SEGMENTATIONSTRIPXALG_H
#define SEGMENTATIONSTRIPXALG_H

#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"

#include <string>
#include <vector>
#include <array>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {
        namespace seg {

            class SegmentationStripXAlg: public SegmentationAlg {

            public:
                SegmentationStripXAlg(fhicl::ParameterSet const& pset);

                SegmentationStripXAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

                ~SegmentationStripXAlg();

                void reconfigure(fhicl::ParameterSet const& pset) override;

                void Initialize(const gar::geo::GeometryCore& geo) override;

                std::array<double, 3> GetPosition(const gar::geo::GeometryCore& geo, const gar::raw::CellID_t& cID) const override;

                gar::raw::CellID_t GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const override;

                bool isTile(const gar::raw::CellID_t& cID) const override;

                bool isBarrel(const gar::raw::CellID_t& cID) const override;

                double getStripLength(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::pair<TVector3, TVector3> getStripEnds(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::pair<float, float> CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const gar::raw::CellID_t& cID) const override;

                std::array<double, 3> ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t& cID) const override;

                const double& stripSizeX() const override { return _stripSizeX; }

                const double& layerDimX() const { return _layer_dim_X; }

                const double& layerDimY() const { return _layer_dim_Y; }

                const std::string& fieldNameX() const { return _xId; }

                const std::string& fieldNameY() const { return _yId; }

                const std::string& fieldNameLayer() const { return _layerId; }

                const std::string& fieldNameSlice() const { return _sliceId; }

                const unsigned int& nLayers() const override { return _nLayers; }

                void setStripSizeX(double stripSize) { _stripSizeX = stripSize; }

                void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

                void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

                void setFieldNameLayer(const std::string& fieldName) { _layerId = fieldName; }

                void setFieldNameSlice(const std::string& fieldName) { _sliceId = fieldName; }

                void setLayerDimXY(const double& dimX, const double& dimY) const override { _layer_dim_X = dimX; _layer_dim_Y = dimY; }
		//unused variables: innerangle, endcapsidelength
                void setVariables(const double& , const double & ) const override { /* no op */ }

            protected:

                void PrintParameters() const override;

                /// the field name used for X
                std::string _xId;
                /// the field name used for Y
                std::string _yId;
                /// the field name used for layer
                std::string _layerId;
                /// the field name used for slice
                std::string _sliceId;
                /// the encoding string
                std::string _encoding;
                /// the strip size in X
                double _stripSizeX;
                /// fraction of tiles to remove at the edge
                double _frac;
                /// number of layers
                unsigned int _nLayers;
                /// layer dimension in X
                mutable double _layer_dim_X;
                /// layer dimension in Y
                mutable double _layer_dim_Y;
            };
        }
    }
}

#endif
