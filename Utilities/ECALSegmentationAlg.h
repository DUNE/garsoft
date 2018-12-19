#ifndef ECALSEGMENTATIONALG_H
#define ECALSEGMENTATIONALG_H

#include "Utilities/BitFieldCoder.h"

#include <vector>
#include <map>

#include "Math/Vector3D.h"

namespace fhicl{
  class ParameterSet;
}

namespace util {

  typedef long long long64 ;
  typedef unsigned long long ulong64 ;

  class ECALSegmentationAlg {

  public:
    virtual ~ECALSegmentationAlg();

    virtual const std::string& name() const {
      return _name;
    }

    virtual void setName(const std::string& value) {
      _name = value;
    }

    virtual const std::string& type() const {
      return _type;
    }

    virtual const std::string& description() const {
      return _description;
    }

    virtual const BitFieldCoder* decoder()  const {
      return _decoder;
    }

    virtual void setDecoder(const BitFieldCoder* decoder);

    virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

    virtual long64 cellID(const unsigned int det_id, const unsigned int stave, const unsigned int module, const unsigned int layer, unsigned int slice, const ROOT::Math::XYZVector& localPosition) const = 0;

    virtual ROOT::Math::XYZVector position(const long64& cID) const = 0;

    virtual void PrintParameters() const = 0;

  protected:
    ECALSegmentationAlg(fhicl::ParameterSet const& pset);

    ECALSegmentationAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

    static double binToPosition(long64 bin, double cellSize, double offset = 0);

    static int positionToBin(double position, double cellSize, double offset = 0);

    std::string _name;

    std::string _type;

    std::string _description;

    const BitFieldCoder* _decoder = 0;

  private:
    ECALSegmentationAlg(const ECALSegmentationAlg&);
  };

}

#endif
