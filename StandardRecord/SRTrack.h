#ifndef CAFSRTRACK_H
#define CAFSRTRACK_H

#include "StandardRecord/SRTrackTrajectory.h"
#include <cstddef>

namespace caf
{
  class SRTrack
  {
  public:
    SRTrack();

    // track data
    size_t id;
    SRVector3D start;
    SRVector3D startp;
    int startq;

    SRVector3D end;
    SRVector3D endp;
    int endq;

    int NTPCClustersOnTrack;

    SRTrackTrajectory fwd, bak;

    int   mcidx;  // Branch index (NOT the GEANT track ID) of MCParticle
    float mcfrac; // that best matches & fraction of ionization therefrom
  };
}

#endif
