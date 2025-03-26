#include "SummaryDataProducts/GeometryConfigurationInfo.h"
// C/C++ standard library
#include <ostream>

// -----------------------------------------------------------------------------
std::ostream& gar::sumdata::operator<<(std::ostream& out,
                                       gar::sumdata::GeometryConfigurationInfo const& info)
{
  if (!info.isDataValid()) return out << "Invalid geometry configuration information" << std::endl;

  out << "Geometry information version: " << info.dataVersion;

  if (info.dataVersion >= gar::sumdata::GeometryConfigurationInfo::DataVersion_t{1})
    out << "\nDetector name:               '" << info.detectorName << "'";

  if (info.dataVersion >= gar::sumdata::GeometryConfigurationInfo::DataVersion_t{2}) {
    out << "\nFull configuration:"
        << "\n"
        << std::string(80, '-') << "\n"
        << info.geometryServiceConfiguration << "\n"
        << std::string(80, '-');
  }

  if (info.dataVersion > gar::sumdata::GeometryConfigurationInfo::DataVersion_t{2}) {
    out << "\n[this version of code can't fully decode further information]";
  }

  return out;
} // operator<< (GeometryConfigurationInfo const&)

// -----------------------------------------------------------------------------