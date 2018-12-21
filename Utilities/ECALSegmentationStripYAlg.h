#ifndef ECALSEGMENTATIONSTRIPYALG_H
#define ECALSEGMENTATIONSTRIPYALG_H

#include "Utilities/ECALSegmentationAlg.h"

#include <string>
#include <vector>

namespace fhicl{
  class ParameterSet;
}

namespace util {

  class ECALSegmentationStripYAlg: public ECALSegmentationAlg {

  public:
    ECALSegmentationStripYAlg(fhicl::ParameterSet const& pset);

    ECALSegmentationStripYAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

    virtual ~ECALSegmentationStripYAlg();

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    virtual ROOT::Math::XYZVector position(const long64& cID) const;

    virtual long64 cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const;

    virtual void PrintParameters() const;

    virtual void SetLayerParameters(TGeoShape *shape) { return; }

    double stripSizeY() const { return _stripSizeY; }

    double offsetY() const { return _offsetY; }

    const std::string& fieldNameY() const { return _yId; }

    void setStripSizeY(double cellSize) { _stripSizeY = cellSize; }

    void setOffsetY(double offset) { _offsetY = offset; }

    void setFieldNameY(const std::string& fieldName) { _yId = fieldName; }

  protected:
    /// the strip size in Y
    double _stripSizeY;
    /// the coordinate offset in Y
    double _offsetY;
    /// the field name used for Y
    std::string _yId;
  };
}

#endif
