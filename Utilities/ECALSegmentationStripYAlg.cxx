#include "Utilities/ECALSegmentationStripYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

namespace util {

    /// Default constructor used by derived classes passing the encoding string
    ECALSegmentationStripYAlg::ECALSegmentationStripYAlg(fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(pset)
    {
        _type = "StripY";
        _description = "Cartesian segmentation in the local Y-plane";

        std::cout << " ######### util::ECALSegmentationStripYAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    /// Default constructor used by derived classes passing an existing decoder
    ECALSegmentationStripYAlg::ECALSegmentationStripYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(decode, pset)
    {
        _type = "StripY";
        _description = "Cartesian segmentation in the local Y-plane";

        std::cout << " ######### util::ECALSegmentationStripYAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    ECALSegmentationStripYAlg::~ECALSegmentationStripYAlg()
    {
    }

    void ECALSegmentationStripYAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
        _yId = pset.get<std::string>("identifier_y");
        _stripSizeY = pset.get<double>("strip_size_y");
        _offsetY = pset.get<double>("offset_y");

        PrintParameters();

        return;
    }

    ROOT::Math::XYZVector ECALSegmentationStripYAlg::position(const long64& cID) const
    {
        ROOT::Math::XYZVector cellPosition;
        cellPosition.SetXYZ(0., binToPosition(_decoder->get(cID, _yId), _stripSizeY, _offsetY), 0.);

        return cellPosition;
    }

    /// determine the cell ID based on the position
    long64 ECALSegmentationStripYAlg::cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const
    {
        long64 cID = 0;

        _decoder->set(cID, "system", det_id);
        _decoder->set(cID, "module", module);
        _decoder->set(cID, "stave", stave);
        _decoder->set(cID, "layer", layer);
        _decoder->set(cID, "slice", slice);

        _decoder->set(cID, _yId, positionToBin(localPosition.Y(), _stripSizeY, _offsetY));

        return cID;
    }

    void ECALSegmentationStripYAlg::PrintParameters() const
    {
        std::cout << "identifier_y: " << _yId << std::endl;
        std::cout << "strip_size_y: " << _stripSizeY << " mm" << std::endl;
        std::cout << "offset_y: " << _offsetY << " mm" << std::endl;
    }

} // util
