#ifndef ECALSEGMENTATIONSTRIPXALG_H
#define ECALSEGMENTATIONSTRIPXALG_H

#include "Utilities/ECALSegmentationAlg.h"

#include <string>
#include <vector>

namespace fhicl{
  class ParameterSet;
}

namespace util {

  class ECALSegmentationStripXAlg: public ECALSegmentationAlg {

  public:
    ECALSegmentationStripXAlg(fhicl::ParameterSet const& pset);

    ECALSegmentationStripXAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

    virtual ~ECALSegmentationStripXAlg();

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    virtual ROOT::Math::XYZVector position(const long64& cID) const;

    virtual long64 cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const;

    virtual void PrintParameters() const;

    virtual void SetLayerParameters(TGeoShape *shape) { return; }

    double stripSizeX() const { return _stripSizeX; }

    double offsetX() const { return _offsetX; }

    const std::string& fieldNameX() const { return _xId; }

    void setStripSizeX(double cellSize) { _stripSizeX = cellSize; }

    void setOffsetX(double offset) { _offsetX = offset; }

    void setFieldNameX(const std::string& fieldName) { _xId = fieldName; }

  protected:
    /// the strip size in X
    double _stripSizeX;
    /// the coordinate offset in X
    double _offsetX;
    /// the field name used for X
    std::string _xId;
  };
}

#endif
