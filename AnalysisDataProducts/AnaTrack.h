#ifndef GAR_ANATRACK_H
#define GAR_ANATRACK_H

#include "ReconstructionDataProducts/Track.h"
#include "garana/DataProducts/Track.h"
#include <vector>

using std::pair;
using std::vector;

namespace gar {
  namespace adp {
    class AnaTrack {

      private:
        rec::Track              fTrack;  // copy all reco::Track info (for now...)
        vector<pair<int,float>> fPidF; ///< PID: (PDG code, probability) front fit
        vector<pair<int,float>> fPidB; ///<  " backfit
        float                   fIonF; ///< average ionization along foward fit
        float                   fIonB; ///< average ionization along backward fit
        int                     fChgF; ///< charge (+/- 1e) if vertex at correct end
        int                     fChgB; ///< charge (+/- 1e) if end is actually vertex

      public:

        // constructors
        AnaTrack() {}
        AnaTrack(const rec::Track& trk, const vector<pair<int,float>>& pidf, 
                 const vector<pair<int,float>>& pidb, float ionf, float ionb) :
          fTrack(trk),
          fPidF(pidf),
          fPidB(pidb),
          fIonF(ionf),
          fIonB(ionb),
          fChgF(trk.ChargeBeg()),
          fChgB(trk.ChargeEnd())
          {}

        // getters
        const rec::Track* GetTrack() const { return &fTrack;      }
        garana::Track GaranaTrack();
        const rec::IDNumber TrackID() const { return fTrack.getIDNumber(); }
        const vector<pair<int,float>>*  PidF() const { return &fPidF;  }
        const vector<pair<int,float>>*  PidB() const { return &fPidB;  }
        //const float  PidProbF() const { return fPidF.second; }
        //const float  PidProbB() const { return fPidB.second; }
        const float  IonizF()   const { return fIonF;        }
        const float  IonizB()   const { return fIonB;        }
    };//class
  }//adp
}//gar

#endif