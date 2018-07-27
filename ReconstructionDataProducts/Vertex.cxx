//
//  Vertex.cxx
//  garsoft-mrb
//
//  Created by Tom Junk on July 12, 2018

#include "ReconstructionDataProducts/Vertex.h"

namespace gar {
  namespace rec {
    
    //--------------------------------------------------------------------------
    Vertex::Vertex(float       *pos,
                   float       *covmat,
		   ULong64_t time)
      : fTime(time)
    {
      size_t icounter = 0;
      for (size_t i=0;i<3;++i)
	{
	  fPosition[i] = pos[i];
	  for (size_t j=0; j<3; ++j)
	    {
	      fCovMat[i][j] = covmat[icounter];
	      ++icounter;
	    }
	}
    }

    //--------------------------------------------------------------------------
    // empty constructor

    Vertex::Vertex() 
      : fTime(0)
    {
      for (size_t i=0;i<3;++i)
	{
	  fPosition[i] = 0;
	  for (size_t j=0; j<3; ++j)
	    {
	      fCovMat[i][j] = 0;
	    }
	}
    }

  } // rec
} // gar