#include "StandardRecord/SRMuIDCluster.h"

#include <limits>
#include <cstddef>

// This value will cause an error as soon as it is used, which should help find
// uninitialized variables.
const float kNaN = std::numeric_limits<float>::signaling_NaN();

namespace caf
{
  SRMuIDCluster::SRMuIDCluster() :
    id(-1), nhits(0), E(kNaN), t(kNaN), TimeDiffFirstLast(kNaN), x(kNaN), y(kNaN), z(kNaN), theta(kNaN), phi(kNaN), pid(kNaN), mcidx(-1), mcfrac(kNaN)
  {
  }
}
