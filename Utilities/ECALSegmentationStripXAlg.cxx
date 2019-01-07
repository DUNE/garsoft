#include "Utilities/ECALSegmentationStripXAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

namespace util {

    /// Default constructor used by derived classes passing the encoding string
    ECALSegmentationStripXAlg::ECALSegmentationStripXAlg(fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(pset)
    {
        _type = "StripX";
        _description = "Cartesian segmentation in the local X-plane";

        std::cout << " ######### util::ECALSegmentationStripXAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    /// Default constructor used by derived classes passing an existing decoder
    ECALSegmentationStripXAlg::ECALSegmentationStripXAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(decode, pset)
    {
        _type = "StripX";
        _description = "Cartesian segmentation in the local X-plane";

        std::cout << " ######### util::ECALSegmentationStripXAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    ECALSegmentationStripXAlg::~ECALSegmentationStripXAlg()
    {
    }

    void ECALSegmentationStripXAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
        _xId = pset.get<std::string>("identifier_x");
        _stripSizeX = pset.get<double>("strip_size_x");
        _offsetX = pset.get<double>("offset_x");

        PrintParameters();

        return;
    }

    ROOT::Math::XYZVector ECALSegmentationStripXAlg::position(const long64& cID) const
    {
        ROOT::Math::XYZVector cellPosition;
        cellPosition.SetXYZ(binToPosition(_decoder->get(cID, _xId), _stripSizeX, _offsetX), 0., 0.);

        return cellPosition;
    }

    /// determine the cell ID based on the position
    long64 ECALSegmentationStripXAlg::cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const
    {
        long64 cID = 0;

        _decoder->set(cID, "system", det_id);
        _decoder->set(cID, "module", module);
        _decoder->set(cID, "stave", stave);
        _decoder->set(cID, "layer", layer);
        _decoder->set(cID, "slice", slice);

        _decoder->set(cID, _xId, positionToBin(localPosition.X(), _stripSizeX, _offsetX));

        return cID;
    }

    void ECALSegmentationStripXAlg::PrintParameters() const
    {
        std::cout << "identifier_x: " << _xId << std::endl;
        std::cout << "strip_size_x: " << _stripSizeX << " mm" << std::endl;
        std::cout << "offset_x: " << _offsetX << " mm" << std::endl;
    }

} // util
