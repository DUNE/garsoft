#ifndef CAFSRTRACKTRAJECTORY_H
#define CAFSRTRACKTRAJECTORY_H

#include "garsoft/StandardRecord/SRTrackPID.h"
#include "garsoft/StandardRecord/SRVector3D.h"

#include <vector>

namespace caf
{
  class SRTrackTrajectory
  {
  public:
    SRTrackTrajectory();

    float len;
    float chi2;
    float avgion; ///< average ionization

    std::vector<SRTrackPID> pid;

    std::vector<SRVector3D> traj;
  };
}

#endif
