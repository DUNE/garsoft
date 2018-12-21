#include "Utilities/ECALSegmentationMultiGridStripXYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

#include "TGeoBBox.h"

namespace util {

    /// Default constructor used by derived classes passing the encoding string
    ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(pset)
    {
        _type = "MultiGridStripXY";
        _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

        std::cout << " ######### util::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    /// Default constructor used by derived classes passing an existing decoder
    ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
    : ECALSegmentationAlg(decode, pset)
    {
        _type = "MultiGridStripXY";
        _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

        std::cout << " ######### util::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

        this->reconfigure(pset);
    }

    ECALSegmentationMultiGridStripXYAlg::~ECALSegmentationMultiGridStripXYAlg()
    {
    }

    void ECALSegmentationMultiGridStripXYAlg::reconfigure(fhicl::ParameterSet const& pset)
    {
        _gridSizeX = pset.get<double>("grid_size_x");
        _gridSizeY = pset.get<double>("grid_size_y");

        _stripSizeX = pset.get<double>("strip_size_x");
        _stripSizeY = pset.get<double>("strip_size_y");

        _offsetX = pset.get<double>("offset_x");
        _offsetY = pset.get<double>("offset_y");

        _xId = pset.get<std::string>("identifier_x");
        _yId = pset.get<std::string>("identifier_y");
        _layerId = pset.get<std::string>("identifier_layer");
        _sliceId = pset.get<std::string>("identifier_slice");

        _gridLayers = pset.get< std::vector<std::string> >("grid_layers");
        _stripLayers = pset.get< std::vector<std::string> >("strip_layers");

        TokenizeLayerVectors(_gridLayers, _stripLayers);

        PrintParameters();

        return;
    }

    ROOT::Math::XYZVector ECALSegmentationMultiGridStripXYAlg::position(const long64& cID) const
    {
        ROOT::Math::XYZVector cellPosition;

        //Need to differentiate case tile and strips based on layer and slice
        CheckLayerConfiguration(_decoder->get(cID, _layerId));

        if(_isTile){
            cellPosition.SetXYZ(binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX), binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY), 0.);
        }
        else{
            int cellIndexX = _decoder->get(cID,_xId);
            int cellIndexY = _decoder->get(cID,_yId);

            if(_decoder->get(cID, _sliceId) == 2)
            cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / 1. ), ( cellIndexY + 0.5 ) * (_layer_dim_Y / (_layer_dim_Y / _stripSizeY ) ), 0.);
            if(_decoder->get(cID, _sliceId) == 3)
            cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / (_layer_dim_X / _stripSizeX ) ), ( cellIndexY + 0.5 ) * (_layer_dim_Y / 1. ), 0.);
        }

        return cellPosition;
    }

    /// determine the cell ID based on the position
    long64 ECALSegmentationMultiGridStripXYAlg::cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const
    {
        long64 cID = 0;

        _decoder->set(cID, "system", det_id);
        _decoder->set(cID, "module", module);
        _decoder->set(cID, "stave", stave);
        _decoder->set(cID, "layer", layer);
        _decoder->set(cID, "slice", slice);

        //Need to differentiate case tile and strips based on layer and slice
        CheckLayerConfiguration(layer);

        if(_isTile)
        {
            _decoder->set(cID, _xId, positionToBin(localPosition.X(), _gridSizeX, _offsetX));
            _decoder->set(cID, _yId, positionToBin(localPosition.Y(), _gridSizeY, _offsetY));
        }
        else{
            if(slice == 2)
            {
                //Segmentation in Y
                int _cellIndexX = int ( localPosition.X() / ( _layer_dim_X / 1. ) );
                int _cellIndexY = int ( localPosition.Y() / ( _layer_dim_Y / (_layer_dim_Y / _stripSizeY) ) );

                _decoder->set(cID, _xId, _cellIndexX);
                _decoder->set(cID, _yId, _cellIndexY);
            }

            if(slice == 3)
            {
                //Segmentation in X
                int _cellIndexX = int ( localPosition.X() / ( _layer_dim_X / (_layer_dim_X / _stripSizeX) ) );
                int _cellIndexY = int ( localPosition.Y() / ( _layer_dim_Y / 1. ) );

                _decoder->set(cID, _xId, _cellIndexX);
                _decoder->set(cID, _yId, _cellIndexY);
            }
        }

        return cID;
    }

    void ECALSegmentationMultiGridStripXYAlg::PrintParameters() const
    {
        std::cout << "identifier_x: " << _xId << std::endl;
        std::cout << "identifier_y: " << _yId << std::endl;
        std::cout << "grid_size_x: " << _gridSizeX << " mm" << std::endl;
        std::cout << "grid_size_y: " << _gridSizeY << " mm" << std::endl;
        std::cout << "offset_x: " << _offsetX << " mm" << std::endl;
        std::cout << "offset_y: " << _offsetY << " mm" << std::endl;

        std::cout << "identifier_layer: " << _layerId << std::endl;
        std::cout << "identifier_slice: " << _sliceId << std::endl;

        std::cout << "grid start layers: ";
        for(unsigned int i = 0; i < _gridFirst.size(); i++)
        std::cout << _gridFirst.at(i) << " ";
        std::cout << std::endl;

        std::cout << "grid last layers: ";
        for(unsigned int i = 0; i < _gridLast.size(); i++)
        std::cout << _gridLast.at(i) << " ";
        std::cout << std::endl;

        std::cout << "strip start layers: ";
        for(unsigned int i = 0; i < _stripFirst.size(); i++)
        std::cout << _stripFirst.at(i) << " ";
        std::cout << std::endl;

        std::cout << "strip last layers: ";
        for(unsigned int i = 0; i < _stripLast.size(); i++)
        std::cout << _stripLast.at(i) << " ";
        std::cout << std::endl;
    }

    void ECALSegmentationMultiGridStripXYAlg::TokenizeLayerVectors(std::vector<std::string> grid, std::vector<std::string> strip) const
    {
        for(unsigned int i = 0; i < grid.size(); i++)
        {
            std::vector<std::string> descriptors ;
            Tokenizer t( descriptors ,':') ;

            std::string value = grid.at(i);
            std::for_each( value.begin(), value.end(), t ) ;

            _gridFirst.push_back(std::atoi(descriptors.at(0).c_str()));
            _gridLast.push_back(std::atoi(descriptors.at(1).c_str()));
        }

        for(unsigned int i = 0; i < strip.size(); i++)
        {
            std::vector<std::string> descriptors ;
            Tokenizer t( descriptors ,':') ;

            std::string value = strip.at(i);
            std::for_each( value.begin(), value.end(), t ) ;

            _stripFirst.push_back(std::atoi(descriptors.at(0).c_str()));
            _stripLast.push_back(std::atoi(descriptors.at(1).c_str()));
        }
    }

    void ECALSegmentationMultiGridStripXYAlg::CheckLayerConfiguration(const unsigned int layer) const
    {
        _isTile = false;
        //Check if it is tile configuration
        for(unsigned int i = 0; i < _gridFirst.size(); i++)
        {
            if( layer >= _gridFirst.at(i) && layer <= _gridLast.at(i) )
            {
                _isTile = true;
                break;
            }
            else{
                continue;
            }
        }
    }

    void ECALSegmentationMultiGridStripXYAlg::SetLayerParameters(TGeoShape *shape)
    {
        //Get layer length in X and Y
        const TGeoBBox *box = static_cast<TGeoBBox*>(shape);
        setLayerDimXY(box->GetDX()*2 * CM_2_MM, box->GetDY()*2 * CM_2_MM); // in mm
    }

} // util
