#include "Geometry/ECALSegmentationMultiGridStripXYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include "Geometry/GeometryCore.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

#include "TGeoBBox.h"
#include "TString.h"

namespace gar {
    namespace geo {

        /// Default constructor used by derived classes passing the encoding string
        ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(fhicl::ParameterSet const& pset)
        : ECALSegmentationAlg(pset)
        {
            _type = "MultiGridStripXY";
            _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

            std::cout << " ######### gar::geo::seg::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

            this->reconfigure(pset);
        }

        /// Default constructor used by derived classes passing an existing decoder
        ECALSegmentationMultiGridStripXYAlg::ECALSegmentationMultiGridStripXYAlg(const BitFieldCoder* decode, fhicl::ParameterSet const& pset)
        : ECALSegmentationAlg(decode, pset)
        {
            _type = "MultiGridStripXY";
            _description = "Cartesian segmentation in the local XY-plane, containing integer number of tiles/strips/cells";

            std::cout << " ######### gar::geo::seg::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;

            this->reconfigure(pset);
        }

        ECALSegmentationMultiGridStripXYAlg::~ECALSegmentationMultiGridStripXYAlg()
        {
        }

        void ECALSegmentationMultiGridStripXYAlg::Initialize(const gar::geo::GeometryCore& geo) {
            UpdateGeometrySegmentation(geo);
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

            _gridBarrelLayers = pset.get< std::vector<std::string> >("grid_barrel_layers");
            _stripBarrelLayers = pset.get< std::vector<std::string> >("strip_barrel_layers");

            _gridEndcapLayers = pset.get< std::vector<std::string> >("grid_endcap_layers");
            _stripEndcapLayers = pset.get< std::vector<std::string> >("strip_endcap_layers");

            PrintParameters();

            return;
        }

        ROOT::Math::XYZVector ECALSegmentationMultiGridStripXYAlg::position(const gar::geo::GeometryCore& geo, const long64& cID) const
        {
            ROOT::Math::XYZVector cellPosition;

            //Need to differentiate case tile and strips based on layer and slice
            CheckLayerConfiguration(geo, _decoder->get(cID, "system"), _decoder->get(cID, "stave"), _decoder->get(cID, "module"), _decoder->get(cID, _layerId));

            if(_isTile){
                cellPosition.SetXYZ(binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX), binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY), 0.);
            }
            else{
                int cellIndexX = _decoder->get(cID,_xId);
                int cellIndexY = _decoder->get(cID,_yId);

                if(_decoder->get(cID, _sliceId) == 2)
                {
                    int nCellsX = 1;
                    int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ) - (_layer_dim_X/2), ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ), 0.);
                    // cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ), ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ), 0.);
                }
                if(_decoder->get(cID, _sliceId) == 3)
                {
                    int nCellsX = int(_layer_dim_X / _stripSizeX);
                    int nCellsY = 1;

                    cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ), ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ) - (_layer_dim_Y/2), 0.);
                    // cellPosition.SetXYZ(( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ), ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ), 0.);
                }
            }

            return cellPosition;
        }

        /// determine the cell ID based on the position
        long64 ECALSegmentationMultiGridStripXYAlg::cellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const ROOT::Math::XYZVector& localPosition) const
        {
            long64 cID = 0;

            _decoder->set(cID, "system", det_id);
            _decoder->set(cID, "module", module);
            _decoder->set(cID, "stave", stave);
            _decoder->set(cID, "layer", layer);
            _decoder->set(cID, "slice", slice);

            //Need to differentiate case tile and strips based on layer and slice
            CheckLayerConfiguration(geo, det_id, stave, module, layer);

            double localX = localPosition.X();
            double localY = localPosition.Y();

            if(_isTile)
            {
                _decoder->set(cID, _xId, positionToBin(localX, _gridSizeX, _offsetX));
                _decoder->set(cID, _yId, positionToBin(localY, _gridSizeY, _offsetY));
            }
            else{
                if(slice == 2)
                {
                    //Segmentation in Y
                    int nCellsX = 1;
                    int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    _decoder->set(cID, _xId, _cellIndexX);
                    _decoder->set(cID, _yId, _cellIndexY);
                }

                if(slice == 3)
                {
                    //Segmentation in X
                    int nCellsX = int(_layer_dim_X / _stripSizeX);
                    int nCellsY = 1;

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    _decoder->set(cID, _xId, _cellIndexX);
                    _decoder->set(cID, _yId, _cellIndexY);
                }
            }

            return cID;
        }

        int ECALSegmentationMultiGridStripXYAlg::getIDbyCellID(const long64& cID, const char* id) const
        {
            return _decoder->get(cID, id);
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

            std::cout << "grid_barrel_layers: ";
            for(unsigned int i = 0; i < _gridBarrelLayers.size(); i++)
            std::cout << _gridBarrelLayers.at(i) << " ";
            std::cout << std::endl;

            std::cout << "strip_barrel_layers: ";
            for(unsigned int i = 0; i < _stripBarrelLayers.size(); i++)
            std::cout << _stripBarrelLayers.at(i) << " ";
            std::cout << std::endl;

            std::cout << "grid_endcap_layers: ";
            for(unsigned int i = 0; i < _gridEndcapLayers.size(); i++)
            std::cout << _gridEndcapLayers.at(i) << " ";
            std::cout << std::endl;

            std::cout << "strip_endcap_layers: ";
            for(unsigned int i = 0; i < _stripEndcapLayers.size(); i++)
            std::cout << _stripEndcapLayers.at(i) << " ";
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

        void ECALSegmentationMultiGridStripXYAlg::CheckLayerConfiguration(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer) const
        {
            if(det_id == 1)
            {
                TokenizeLayerVectors(_gridBarrelLayers, _stripBarrelLayers);
                std::string layer_name = std::string(TString::Format("BarrelECal_stave%02i_module%02i_layer_%02i", stave, module, layer));

                const std::shared_ptr<gar::geo::ECALLayerParamsClass> p = geo.GetECALLayerDimensions( layer_name );

                double dim_x = p->getLayerDimensionX();
                double dim_y = p->getLayerDimensionY();

                setLayerDimXY(dim_x * 2 * CM_2_MM, dim_y * 2 * CM_2_MM); // in mm
            }

            if(det_id == 2)
            {
                TokenizeLayerVectors(_gridEndcapLayers, _stripEndcapLayers);
                std::string layer_name = std::string(TString::Format("EndcapECal_stave%02i_module%02i_layer_%02i", stave, module, layer));
                const std::shared_ptr<gar::geo::ECALLayerParamsClass> p = geo.GetECALLayerDimensions( layer_name );

                double dim_x = p->getLayerDimensionX();
                double dim_y = p->getLayerDimensionY();

                setLayerDimXY(dim_x * 2 * CM_2_MM, dim_y * 2 * CM_2_MM); // in mm
            }

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

        void ECALSegmentationMultiGridStripXYAlg::UpdateGeometrySegmentation(const gar::geo::GeometryCore& geo) const
        {
            //Barrel
            for(int istave = 0; istave < 8; istave++)
            {
                for(int imodule = 0; imodule < 5; imodule++)
                {
                    for(unsigned int ilayer = 0; ilayer < 55; ilayer++)
                    {
                        std::string layer_name = std::string(TString::Format("BarrelECal_stave%02i_module%02i_layer_%02i", istave+1, imodule+1, ilayer+1));
                        geo.UpdateECALLayerSegmentation(layer_name, _gridSizeX, _stripSizeX, this->GetSegmentationConfig(1, ilayer+1));
                    }
                }
            }

            //Endcap
            for(int istave = 0; istave < 4; istave++)
            {
                for(int imodule = -1; imodule < 6; imodule++)
                {
                    for(unsigned int ilayer = 0; ilayer < 45; ilayer++)
                    {
                        std::string layer_name = std::string(TString::Format("EndcapECal_stave%02i_module%02i_layer_%02i", istave+1, imodule+1, ilayer+1));
                        geo.UpdateECALLayerSegmentation(layer_name, _gridSizeX, _stripSizeX, this->GetSegmentationConfig(2, ilayer+1));
                    }
                }
            }

            std::cout << " ######### gar::geo::seg::ECALSegmentationMultiGridStripXYAlg() " << std::endl ;
            std::cout << " ECALSegmentationMultiGridStripXYAlg::UpdateGeometrySegmentation() --- Geometry Map updated with the segmentation " << std::endl;
        }

        bool ECALSegmentationMultiGridStripXYAlg::GetSegmentationConfig(const unsigned int& det_id, const unsigned int& layer) const
        {
            bool isTile = false;

            if(det_id == 1)
            TokenizeLayerVectors(_gridBarrelLayers, _stripBarrelLayers);

            if(det_id == 2)
            TokenizeLayerVectors(_gridEndcapLayers, _stripEndcapLayers);

            //Check if it is tile configuration
            for(unsigned int i = 0; i < _gridFirst.size(); i++)
            {
                if( layer >= _gridFirst.at(i) && layer <= _gridLast.at(i) )
                {
                    isTile = true;
                    break;
                }
                else{
                    continue;
                }
            }

            return isTile;
        }
    } // geo
} //gar
