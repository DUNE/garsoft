#ifndef ECALSEGMENTATIONGRIDXYALG_H
#define ECALSEGMENTATIONGRIDXYALG_H

#include "Geometry/ECALSegmentationAlg.h"
#include "Geometry/GeometryCore.h"

#include <string>
#include <vector>

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {

        class ECALSegmentationGridXYAlg: public ECALSegmentationAlg {

        public:
            ECALSegmentationGridXYAlg(fhicl::ParameterSet const& pset);

            ECALSegmentationGridXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

            virtual ~ECALSegmentationGridXYAlg();

            virtual void reconfigure(fhicl::ParameterSet const& pset);

            virtual void Initialize(const gar::geo::GeometryCore& geo);

            virtual G4ThreeVector position(const gar::geo::GeometryCore& geo, const raw::CellID_t& cID) const;

            virtual raw::CellID_t cellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const G4ThreeVector& localPosition) const;

            virtual int getIDbyCellID(const raw::CellID_t& cID, const char* id) const;

            virtual bool isTile(const raw::CellID_t& cID) const { return true; }

            virtual double getStripLength(const gar::geo::GeometryCore& geo, const raw::CellID_t& cID) const { return 0.; }

            virtual std::pair<float, float> CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3U> &local, const raw::CellID_t& cID) const { return std::make_pair(0., 0.); }

            virtual std::array<double, 3U> ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3U> &local, const float &xlocal, const raw::CellID_t& cID) const { return std::array<double, 3U>{ {0., 0., 0.} }; }

            virtual void PrintParameters() const;

            const double& gridSizeX() const { return _gridSizeX; }

            const double& gridSizeY() const { return _gridSizeY; }

            const double& offsetX() const { return _offsetX; }

            const double& offsetY() const { return _offsetY; }

            const std::string& fieldNameX() const { return _xId; }

            const std::string& fieldNameY() const { return _yId; }

            void setGridSizeX(double cellSize) { _gridSizeX = cellSize; }

            void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }

            void setOffsetX(double offset) { _offsetX = offset; }

            void setOffsetY(double offset) { _offsetY = offset; }

            void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

            void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

            void setLayerDimXY(const double& dimX, const double& dimY) const { return; }

        protected:
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
        };
    }
}

#endif
