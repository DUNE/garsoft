//
//  TPCCluster.cxx
//
//  Created by Brian Rebel on 10/6/16.
//
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "ReconstructionDataProducts/TPCCluster.h"
#include "TMath.h"

namespace gar {
  namespace rec {
   
    //--------------------------------------------------------------------------
    TPCCluster::TPCCluster()
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();
      return;
    }



    bool TPCCluster::operator==(const TPCCluster& rhs) const {
        return (this->fIDnumero == rhs.fIDnumero);
    }

    bool TPCCluster::operator!=(const TPCCluster& rhs) const {
        return (this->fIDnumero != rhs.fIDnumero);
    }

    gar::rec::IDNumber TPCCluster::getIDNumber() const {return fIDnumero;}



    //--------------------------------------------------------------------------
    TPCCluster::TPCCluster( float        sig,
         float       *pos,
         float        startT,
         float        endT,
	     float        Time,
	     float        RMS)
    : fSignal   (sig   )
    , fTime     (Time  )
    , fStartTime(startT)
    , fEndTime  (endT  )
    , fRMS      (RMS   )
    {
      IDNumberGen::create(FirstNumber);
      fIDnumero = IDNumberGen::create()->getNewOne();
      
      fPosition[0] = pos[0];
      fPosition[1] = pos[1];
      fPosition[2] = pos[2];
      
      return;
    }
    
    //--------------------------------------------------------------------------
    // Users should be careful with this method.  It will just add TPCClusters 
    // without checking if the TPCClusters are overlapping in time or close
    // to each other in time
    void TPCCluster::operator+=(gar::rec::TPCCluster const& h)
    {
      
      // add h to this TPCCluster.  Set the new position, etc to be the weighted
      // average of the individual TPCClusters
      float totSig = fSignal + h.Signal();
      
      if(totSig == 0.){
        LOG_WARNING("TPCCluster")
        << "attempting to add two TPCClusters and neithr has any signal, bail.";
        return;
      }
 
      for(size_t i = 0; i < 3; ++i)
        fPosition[i] = (fPosition[i] * fSignal + h.Position()[i] * h.Signal()) / totSig;

      // don't do a weighted average but just make a wider TPCCluster
      
      //fStartTime = (fStartTime * fSignal + h.StartTime() * h.Signal()) / totSig;
      //fEndTime   = (fEndTime   * fSignal + h.EndTime()   * h.Signal()) / totSig;
      
      fStartTime = TMath::Min(fStartTime,h.StartTime());
      fEndTime   = TMath::Max(fEndTime,h.EndTime());

      float avgtime = (fTime * fSignal + h.Time() * h.Signal()) / totSig;

      fRMS = TMath::Sqrt(  (fSignal*(TMath::Sq(fTime-avgtime)+TMath::Sq(fRMS))
			    + h.Signal()*(TMath::Sq(h.Time()-avgtime)+TMath::Sq(h.RMS())))/totSig );
                          
      fTime      = avgtime;

      fSignal += h.Signal();
      
      return;
    }


    //--------------------------------------------------------------------------
    std::ostream& operator<< (std::ostream& o, gar::rec::TPCCluster const& h)
    {
      
      o << "TPCCluster "
      << "\n\tID number = "
      << h.getIDNumber()
      << "\n\tposition = ("
      << h.Position()[0]
      << ", "
      << h.Position()[1]
      << ", "
      << h.Position()[2]
      << ")"
      << "\n\tsignal = "
      << h.Signal()
      << "\n\tstart time: "
      << h.StartTime()
      << " end time: "
      << h.EndTime();
      
      return o;
    }

  } // rec
} // gar
