//
//  TPCCluster.h
//
//  Created by Brian Rebel on 10/6/16.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_TPCCluster_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_TPCCluster_h

#include <iostream>
#include "IDNumberGen.h"

namespace gar {
  namespace rec {

    class TPCCluster {
    
    public:
      TPCCluster();
      
      // let the compiler provide the dtor
      
    private:
      static gar::rec::IDNumber const FirstNumber = 100100000;
      gar::rec::IDNumber fIDnumero;
      
      float        fSignal;      ///< size of the signal for this TPCCluster  (integral of ADC values)
      float        fPosition[3]; ///< position of the TPCCluster
      float        fTime;        ///< time of TPCCluster charge arrival at the readout plane (ticks)
      float        fStartTime;   ///< start time of the TPCCluster (ticks)
      float        fEndTime;     ///< end time of the TPCCluster (ticks)
      float        fRMS;         ///< TPCCluster width calculated with RMS (in ticks)

#ifndef __GCCXML__
      
    public:
      
      TPCCluster(
          float        sig,
          float       *pos,
          float        startT,
          float        endT,
	      float        Time,
	      float        RMS);
      
      bool operator==(const TPCCluster& rhs) const;
      bool operator!=(const TPCCluster& rhs) const;
      gar::rec::IDNumber getIDNumber() const;
      
      const float*        Position()  const;
      float        const& Signal()    const;
      unsigned int        Channel()   const;
      float               StartTime() const;
      float               EndTime()   const;
      float               Time()      const;
      float               RMS()       const;

      void operator += (gar::rec::TPCCluster const& h);
      
      friend std::ostream& operator << (std::ostream & o, gar::rec::TPCCluster const& h);
      
#endif
      
    };
    
    inline const  float* TPCCluster::Position()  const { return &fPosition[0]; }
    inline float  const& TPCCluster::Signal()    const { return fSignal;       }
    inline float         TPCCluster::StartTime() const { return fStartTime;    }
    inline float         TPCCluster::EndTime()   const { return fEndTime;      }
    inline float         TPCCluster::RMS()       const { return fRMS;          }
    inline float         TPCCluster::Time()      const { return fTime;         }    
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_TPCCluster_h */
