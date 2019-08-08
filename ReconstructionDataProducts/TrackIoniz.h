//
//  TrackIoniz.h
//
//  Created by Leo Bellantoni on 22 Feb 2019.
//  Track ionization data -- use associations with tracks to find them.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_TrackIoniz_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_TrackIoniz_h

#include <vector>



namespace gar {
  namespace rec {

    class TrackIoniz {

    public:
      TrackIoniz();
      inline const std::vector<std::pair<float,float>> getFWD_dSigdXs() {return fFWD_dSigdXs;};
      inline const std::vector<std::pair<float,float>> getBAK_dSigdXs() {return fBAK_dSigdXs;};

      // let the compiler provide the dtor

    private:
      std::vector<std::pair<float,float>> fFWD_dSigdXs;  ///< Segment values ordered by forward fit
      std::vector<std::pair<float,float>> fBAK_dSigdXs;  ///< Segment values ordered by backward fit



#ifndef __GCCXML__

    public:

      void setData(std::vector<std::pair<float,float>> dSigdX_FWD,std::vector<std::pair<float,float>> dSigdX_BAK);

#endif

    };






  } // rec
} // gar


#endif
