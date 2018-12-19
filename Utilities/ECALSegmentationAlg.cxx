#include "Utilities/ECALSegmentationAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

namespace util {

    ECALSegmentationAlg::ECALSegmentationAlg(fhicl::ParameterSet const& pset) :
      _name("Segmentation"), _type("Segmentation"), _decoder(new BitFieldCoder(pset.get<std::string>("cellEncoding"))) {
    }

    ECALSegmentationAlg::ECALSegmentationAlg(const BitFieldCoder* newDecoder, fhicl::ParameterSet const& pset) :
      _name("Segmentation"), _type("Segmentation"), _decoder(newDecoder) {
    }

    ECALSegmentationAlg::~ECALSegmentationAlg() {
        if (_decoder != 0) {
            delete _decoder;
        }
    }

    void ECALSegmentationAlg::setDecoder(const BitFieldCoder* newDecoder) {
        if ( _decoder == newDecoder )
        return; //self assignment
        else {
            delete _decoder;
        }
        _decoder = newDecoder;
    }

    double ECALSegmentationAlg::binToPosition(long64 bin, double cellSize, double offset) {
        return bin * cellSize + offset;
    }

    int ECALSegmentationAlg::positionToBin(double position, double cellSize, double offset) {
        return int(floor((position + 0.5 * cellSize - offset) / cellSize));
    }

}// util
