#ifndef CAFSRMUIDCLUSTER_H
#define CAFSRMUIDCLUSTER_H

#include "StandardRecord/SRVector3D.h"
#include <cstddef>

#include <vector>

namespace caf
{
  class SRMuIDCluster
  {
  public:
    SRMuIDCluster();

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
