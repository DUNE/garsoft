#ifndef GARSOFT_RECO_TRACKER2ALGS_H
#define GARSOFT_RECO_TRACKER2ALGS_H

#include <vector>
#include <iostream>
#include <cstddef>
#include <TMath.h>
#include <TVector3.h>
#include <ReconstructionDataProducts/TPCCluster.h>

namespace gar {
  namespace rec {

    float capprox(float x1,float y1,
		  float x2,float y2,
		  float x3,float y3,
		  float &xc, float &yc);


    int initial_trackpar_estimate(const std::vector<gar::rec::TPCCluster>  &TPCClusters,
				  std::vector<int> &TPCClusterlist,
				  float &curvature_init,
				  float &lambda_init,
				  float &phi_init,
				  float &xpos,
				  float &ypos,
				  float &zpos,
				  float &x_other_end,
				  unsigned int initialtpnTPCClusters,
				  int printlevel);

    void sort_TPCClusters_along_track(const std::vector<gar::rec::TPCCluster>  &TPCClusters,
				      std::vector<int> &hlf,
				      std::vector<int> &hlb,
				      int printlevel,
				      float &lengthforwards,
				      float &lengthbackwards,
				      float sorttransweight,
				      float sortdistback);

    void sortTPCClusterAux1(std::vector<gar::rec::TPCCluster> &tcc, std::vector<int> &lss, int printlevel);

  }
}


#endif
