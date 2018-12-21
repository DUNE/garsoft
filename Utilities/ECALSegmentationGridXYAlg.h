#ifndef ECALSEGMENTATIONGRIDXYALG_H
#define ECALSEGMENTATIONGRIDXYALG_H

#include "Utilities/ECALSegmentationAlg.h"

#include <string>
#include <vector>

namespace fhicl{
  class ParameterSet;
}

namespace util {

  class ECALSegmentationGridXYAlg: public ECALSegmentationAlg {

  public:
    ECALSegmentationGridXYAlg(fhicl::ParameterSet const& pset);

    ECALSegmentationGridXYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

    virtual ~ECALSegmentationGridXYAlg();

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    virtual ROOT::Math::XYZVector position(const long64& cID) const;

    virtual long64 cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const;

    virtual void PrintParameters() const;

    virtual void SetLayerParameters(TGeoShape *shape) { return; }

    double gridSizeX() const { return _gridSizeX; }

    double gridSizeY() const { return _gridSizeY; }

    double offsetX() const { return _offsetX; }

    double offsetY() const { return _offsetY; }

    const std::string& fieldNameX() const { return _xId; }

    const std::string& fieldNameY() const { return _yId; }

    void setGridSizeX(double cellSize) { _gridSizeX = cellSize; }

    void setGridSizeY(double cellSize) { _gridSizeY = cellSize; }

    void setOffsetX(double offset) { _offsetX = offset; }

    void setOffsetY(double offset) { _offsetY = offset; }

    void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

    void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

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

#endif
