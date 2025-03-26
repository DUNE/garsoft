#ifndef GARSOFT_RECO_TRACKER2ALGS_H
#define GARSOFT_RECO_TRACKER2ALGS_H

#include <ReconstructionDataProducts/TPCCluster.h>
#include <TMath.h>
#include <TVector3.h>
#include <cstddef>
#include <iostream>
#include <vector>

namespace gar {
  namespace rec {

    float capprox(float x1, float y1, float x2, float y2, float x3, float y3, float& xc, float& yc);

    // from OPAL
    void ouchef(double* x,
                double* y,
                int np,
                double& xcirccent,
                double& ycirccent,
                double& rcirc,
                double& chisq,
                int& ifail);

    int initial_trackpar_estimate(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                  std::vector<int>& TPCClusterlist,
                                  float& curvature_init,
                                  float& lambda_init,
                                  float& phi_init,
                                  float& xpos,
                                  float& ypos,
                                  float& zpos,
                                  float& x_other_end,
                                  unsigned int initialtpnTPCClusters,
                                  int printlevel);

    void sort_TPCClusters_along_track(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                      std::vector<int>& hlf,
                                      std::vector<int>& hlb,
                                      int printlevel,
                                      float& lengthforwards,
                                      float& lengthbackwards,
                                      float sorttransweight,
                                      float sortdistback);

    void sort_TPCClusters_along_track2(const std::vector<gar::rec::TPCCluster>& TPCClusters,
                                       std::vector<int>& hlf,
                                       std::vector<int>& hlb,
                                       int printlevel,
                                       float& lengthforwards,
                                       float& lengthbackwards,
                                       float dcut);

    // utilities to get upper-triangular packed vector indices -- go forwards and backwards

    size_t ij2idxsort2(size_t nclus, size_t i, size_t j);

    void idx2ijsort2(size_t nclus, size_t idx, size_t& i, size_t& j);

    bool sort2_check_cyclic(std::vector<int>& link1, std::vector<int>& link2, int i, int j);

  }
}

#endif
