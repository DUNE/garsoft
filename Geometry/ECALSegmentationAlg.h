#ifndef ECALSEGMENTATIONALG_H
#define ECALSEGMENTATIONALG_H

#include "Geometry/BitFieldCoder.h"

#include <vector>
#include <map>

#include "Geometry/GeometryCore.h"

#include "CLHEP/Vector/ThreeVector.h"

typedef CLHEP::Hep3Vector G4ThreeVector;

namespace fhicl{
    class ParameterSet;
}

namespace gar {
    namespace geo {

        class GeometryCore;

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

            virtual const double& gridSizeX() const {
                return _gridSizeX;
            }

            virtual const double& stripSizeX() const {
                return _stripSizeX;
            }

            virtual void setDecoder(const BitFieldCoder* decoder);

            virtual void reconfigure(fhicl::ParameterSet const& pset) = 0;

            virtual long64 cellID(const gar::geo::GeometryCore& geo, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const G4ThreeVector& localPosition) const = 0;

            virtual int getIDbyCellID(const long64& cID, const char* id) const = 0;

            virtual G4ThreeVector position(const gar::geo::GeometryCore& geo, const long64& cID) const = 0;

            virtual void PrintParameters() const = 0;

            virtual void Initialize(const gar::geo::GeometryCore & geo) = 0;

            virtual bool isTile(const unsigned int& det_id, const unsigned int& layer) const = 0;

        protected:
            ECALSegmentationAlg(fhicl::ParameterSet const& pset);

            ECALSegmentationAlg(const BitFieldCoder* decoder, fhicl::ParameterSet const& pset);

            static double binToPosition(long64 bin, double cellSize, double offset = 0);

            static int positionToBin(double position, double cellSize, double offset = 0);

            std::string _name;

            std::string _type;

            std::string _description;

            double _gridSizeX;

            double _stripSizeX;

            const BitFieldCoder* _decoder = 0;

        private:
            ECALSegmentationAlg(const ECALSegmentationAlg&);
        };
    }
}

#endif
