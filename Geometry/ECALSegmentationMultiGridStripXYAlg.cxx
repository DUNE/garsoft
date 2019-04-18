#include "Geometry/ECALSegmentationMultiGridStripXYAlg.h"

#include "fhiclcpp/ParameterSet.h"

#include "Geometry/GeometryCore.h"

#include <algorithm>
#include <sstream>
#include <iostream>
#include <cmath>

#include "TGeoBBox.h"
#include "TString.h"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

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

            _OnSameLayer = pset.get<bool>("strip_on_same_layer");

            this->PrintParameters();

            return;
        }

        void ECALSegmentationMultiGridStripXYAlg::Initialize(const gar::geo::GeometryCore& geo)
        {

        }

        G4ThreeVector ECALSegmentationMultiGridStripXYAlg::position(const gar::geo::GeometryCore& geo, const long64& cID) const
        {
            G4ThreeVector cellPosition;

            //Need to differentiate case tile and strips based on layer and slice
            bool isTile = this->isTile(cID);

            if(isTile){
                cellPosition.setX(binToPosition(_decoder->get(cID, _xId), _gridSizeX, _offsetX));
                cellPosition.setY(binToPosition(_decoder->get(cID, _yId), _gridSizeY, _offsetY));
                cellPosition.setZ(0.);
            }
            else
            {
                int cellIndexX = _decoder->get(cID,_xId);
                int cellIndexY = _decoder->get(cID,_yId);

                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
                {
                    //Segmentation in Y
                    int nCellsX = 1;
                    // int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    cellPosition.setX( ( cellIndexX + 0.5 ) * (_layer_dim_X / nCellsX ) - (_layer_dim_X / 2.) );
                    cellPosition.setY( ( cellIndexY + 0.5 ) * _stripSizeY );
                    cellPosition.setZ(0.);
                }

                if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) )
                {
                    //Segmentation in X
                    // int nCellsX = int(_layer_dim_X / _stripSizeX);
                    int nCellsY = 1;

                    cellPosition.setX( ( cellIndexX + 0.5 ) * _stripSizeX );
                    cellPosition.setY( ( cellIndexY + 0.5 ) * (_layer_dim_Y / nCellsY ) - (_layer_dim_Y / 2.) );
                    cellPosition.setZ(0.);
                }
            }

            return cellPosition;
        }

        /// determine the cell ID based on the position
        long64 ECALSegmentationMultiGridStripXYAlg::cellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const G4ThreeVector& localPosition) const
        {
            long64 cID = 0;

            _decoder->set(cID, "system", det_id);
            _decoder->set(cID, "module", module);
            _decoder->set(cID, "stave", stave);
            _decoder->set(cID, "layer", layer);
            _decoder->set(cID, "slice", slice);

            //Need to differentiate case tile and strips based on layer and slice
            bool isTile = this->isTile(cID);

            double localX = localPosition.x();
            double localY = localPosition.y();

            if(isTile)
            {
                _decoder->set(cID, _xId, positionToBin(localX, _gridSizeX, _offsetX));
                _decoder->set(cID, _yId, positionToBin(localY, _gridSizeY, _offsetY));
            }
            else
            {
                if( (_OnSameLayer && slice == 2) || (not _OnSameLayer && layer%2 == 0) )
                {
                    //Segmentation in Y
                    int nCellsX = 1;
                    int nCellsY = int(_layer_dim_Y / _stripSizeY);

                    int _cellIndexX = int ( localX / ( _layer_dim_X / nCellsX ) );
                    int _cellIndexY = int ( localY / ( _layer_dim_Y / nCellsY ) );

                    _decoder->set(cID, _xId, _cellIndexX);
                    _decoder->set(cID, _yId, _cellIndexY);
                }

                if( (_OnSameLayer && slice == 3) || (not _OnSameLayer && layer%2 != 0) )
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
            std::cout << "strip_size_x: " << _stripSizeX << " mm" << std::endl;
            std::cout << "strip_size_y: " << _stripSizeY << " mm" << std::endl;
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

            std::cout << "strip_on_same_layer: " << _OnSameLayer << std::endl;
        }

        std::array<std::vector<unsigned int>, 2> ECALSegmentationMultiGridStripXYAlg::TokenizeLayerVectors(std::vector<std::string> grid) const
        {
            std::vector<unsigned int> _gridFirst;
            std::vector<unsigned int> _gridLast;

            for(unsigned int i = 0; i < grid.size(); i++)
            {
                std::vector<std::string> descriptors ;
                Tokenizer t( descriptors ,':') ;

                std::string value = grid.at(i);
                std::for_each( value.begin(), value.end(), t ) ;

                _gridFirst.push_back(std::atoi(descriptors.at(0).c_str()));
                _gridLast.push_back(std::atoi(descriptors.at(1).c_str()));
            }

            std::array<std::vector<unsigned int>, 2> _list = {{_gridFirst, _gridLast}};
            return _list;
        }

        bool ECALSegmentationMultiGridStripXYAlg::isTile(const long long int& cID) const
        {
            bool isTile = false;

            unsigned int det_id = _decoder->get(cID, "system");
            unsigned int layer = _decoder->get(cID, _layerId);

            std::array<std::vector<unsigned int>, 2> _list;

            if(det_id == 1)
            _list = TokenizeLayerVectors(_gridBarrelLayers);
            if(det_id == 2)
            _list =  TokenizeLayerVectors(_gridEndcapLayers);

            //Check if it is tile configuration
            for(unsigned int i = 0; i < _list.at(0).size(); i++)
            {
                if( layer >= _list.at(0).at(i) && layer <= _list.at(1).at(i) )
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

        double ECALSegmentationMultiGridStripXYAlg::getStripLength(const gar::geo::GeometryCore& geo, const long64& cID) const
        {
            double stripLength = 0.;

            //Strip along X
            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
            stripLength = _layer_dim_X;
            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) //Strip along Y
            stripLength = _layer_dim_Y;

            return stripLength;
        }

        std::pair<float, float> ECALSegmentationMultiGridStripXYAlg::CalculateLightPropagation(const gar::geo::GeometryCore& geo, const std::array<double, 3U> &local, const long64& cID) const
        {
            //Propagate the light to the SiPM on the side
            //convert c to mm/ns
            // c   = 299.792458 mm/ns
            float c = CLHEP::c_light * CLHEP::mm / CLHEP::ns;
            float time1 = 0.;
            float time2 = 0.;

            //Strip along X
            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
            {
                time1 = ( _layer_dim_X / 2 - std::abs(local[0] * CLHEP::cm) ) / c;
                time2 = ( _layer_dim_X / 2 + std::abs(local[0] * CLHEP::cm) ) / c;
            }
            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) //Strip along Y
            {
                time1 = ( _layer_dim_Y / 2 - std::abs(local[1] * CLHEP::cm) ) / c;
                time2 = ( _layer_dim_Y / 2 + std::abs(local[1] * CLHEP::cm) ) / c;
            }

            return std::make_pair(time1, time2);
        }

        std::array<double, 3U> ECALSegmentationMultiGridStripXYAlg::ReconstructStripHitPosition(const gar::geo::GeometryCore& geo, const std::array<double, 3U> &local, const float &xlocal, const long64& cID) const
        {
            std::array<double, 3U> newlocal;

            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 2) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 == 0) )
            newlocal = { {xlocal, local[1], local[2]} };

            if( (_OnSameLayer && _decoder->get(cID, _sliceId) == 3) || (not _OnSameLayer && _decoder->get(cID, _layerId)%2 != 0) ) //Strip along Y
            newlocal = { {local[0], xlocal, local[2]} };

            return newlocal;
        }
    } // geo
} //gar
