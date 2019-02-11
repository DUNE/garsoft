#ifndef GARSOFT_RECO_TRACKER2ALGS_H
#define GARSOFT_RECO_TRACKER2ALGS_H

#include <vector>
#include <iostream>
#include <cstddef>
#include <TMath.h>
#include <TVector3.h>
#include <ReconstructionDataProducts/Hit.h>

namespace gar {
  namespace rec {

    float capprox(float x1,float y1,
		  float x2,float y2,
		  float x3,float y3,
		  float &xc, float &yc);


    int initial_trackpar_estimate(const std::vector<Hit>  &hits,
				  std::vector<int> &hitlist,
				  float &curvature_init,
				  float &lambda_init,
				  float &phi_init,
				  float &xpos,
				  float &ypos,
				  float &zpos,
				  float &x_other_end,
				  unsigned int initialtpnhits,
				  int printlevel);

    void sort_hits_along_track(const std::vector<Hit>  &hits,
	 	               std::vector<int> &hlf,
			       std::vector<int> &hlb,
			       int printlevel,
			       float &lengthforwards,
			       float &lengthbackwards);

  }
}


#endif
