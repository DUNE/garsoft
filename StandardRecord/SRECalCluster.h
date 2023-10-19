#ifndef CAFSRECALCLUSTER_H
#define CAFSRECALCLUSTER_H

#include "StandardRecord/SRVector3D.h"
#include <vector>
#include <cstddef>

namespace caf
{
  class SRECalCluster
  {
  public:
    SRECalCluster();

    size_t id;
    unsigned int nhits;
    float E;
    float t;
    float TimeDiffFirstLast;
    float x, y, z;
    float theta, phi;
    float pid;
    SRVector3D mainAxis;
    int   mcidx;  // Branch index (NOT the GEANT track ID) of MCPartice
    float mcfrac; // that best matches & fraction of ionization therefrom
  };
}

#endif
