#include "Geometry/ECALSegmentationGridXYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

namespace gar {
    namespace geo {

        //----------------------------------------------------------------------------
        /// Default constructor used by derived classes passing the encoding string
        ECALSegmentationGridXYAlg::ECALSegmentationGridXYAlg(fhicl::ParameterSet const& pset)
        : ECALSegmentationAlg(pset)
        {
            _type = "GridXY";
            _description = "Cartesian segmentation in the local XY-plane";

            std::cout << " ######### util::ECALSegmentationGridXYAlg() " << std::endl ;

            this->reconfigure(pset);
        }

        //----------------------------------------------------------------------------
        /// Default constructor used by derived classes passing an existing decoder
        ECALSegmentationGridXYAlg::ECALSegmentationGridXYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
        : ECALSegmentationAlg(decode, pset)
        {
            _type = "GridXY";
            _description = "Cartesian segmentation in the local XY-plane";

            std::cout << " ######### util::ECALSegmentationGridXYAlg() " << std::endl ;

            this->reconfigure(pset);
        }

        //----------------------------------------------------------------------------
        ECALSegmentationGridXYAlg::~ECALSegmentationGridXYAlg()
        {
        }

        //----------------------------------------------------------------------------
        void ECALSegmentationGridXYAlg::reconfigure(fhicl::ParameterSet const& pset)
        {
            _xId = pset.get<std::string>("identifier_x");
            _yId = pset.get<std::string>("identifier_y");

            _gridSizeX = pset.get<double>("grid_size_x");
            _gridSizeY = pset.get<double>("grid_size_y");

            _offsetX = pset.get<double>("offset_x");
            _offsetY = pset.get<double>("offset_y");

            PrintParameters();

            return;
        }

        //----------------------------------------------------------------------------
        void ECALSegmentationGridXYAlg::Initialize(const gar::geo::GeometryCore& geo)
        {

        }

        //----------------------------------------------------------------------------
        std::array<double, 3> ECALSegmentationGridXYAlg::GetPosition(const gar::geo::GeometryCore& geo, const long64& cID) const
        {
            std::array<double, 3> cellPosition;
            cellPosition[0] = binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX);
            cellPosition[1] = binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY);
            cellPosition[2] = 0.;

            return cellPosition;
        }

        //----------------------------------------------------------------------------
        /// determine the cell ID based on the position
        long64 ECALSegmentationGridXYAlg::GetCellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
        {
            long64 cID = 0;

            _decoder->set(cID, "system", det_id);
            _decoder->set(cID, "module", module);
            _decoder->set(cID, "stave", stave);
            _decoder->set(cID, "layer", layer);
            _decoder->set(cID, "slice", slice);

            _decoder->set(cID, _xId, positionToBin(localPosition[0], _gridSizeX, _offsetX));
            _decoder->set(cID, _yId, positionToBin(localPosition[1], _gridSizeY, _offsetY));

            return cID;
        }

        //----------------------------------------------------------------------------
        int ECALSegmentationGridXYAlg::getIDbyCellID(const long64& cID, const char* id) const
        {
            return _decoder->get(cID, id);
        }

        //----------------------------------------------------------------------------
        void ECALSegmentationGridXYAlg::PrintParameters() const
        {
            std::cout << "identifier_x: " << _xId << std::endl;
            std::cout << "identifier_y: " << _yId << std::endl;
            std::cout << "grid_size_x: " << _gridSizeX << " cm" << std::endl;
            std::cout << "grid_size_y: " << _gridSizeY << " cm" << std::endl;
            std::cout << "offset_x: " << _offsetX << " cm" << std::endl;
            std::cout << "offset_y: " << _offsetY << " cm" << std::endl;
        }
    } //geo
} //gar
