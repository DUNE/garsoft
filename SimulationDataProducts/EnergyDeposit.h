//
//  EnergyDeposit.h
//
//  Created by Brian Rebel on 1/30/17.
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h
#define GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h

#include <list>

namespace gar {
  namespace sdp {


    class EnergyDeposit{
    public:

      EnergyDeposit();

#ifndef __GCCXML__
      EnergyDeposit(int    trackID,
                    float t,
                    float  e,
                    float  xpos,
                    float  ypos,
                    float  zpos,
                    float  length,
                    bool   isPrimary)
      : fTrackID  (trackID)
      , fTime     (t)
      , fEnergy   (e)
      , fX        (xpos)
      , fY        (ypos)
      , fZ        (zpos)
      , fdX       (length)
      , fIsPrimary(isPrimary)
      {}

      int    const& TrackID()   const { return fTrackID;   }
      float  const& Time()      const { return fTime;      }
      float  const& Energy()    const { return fEnergy;    }
      float  const& X()         const { return fX;         }
      float  const& Y()         const { return fY;         }
      float  const& Z()         const { return fZ;         }
      float  const& dX()        const { return fdX;        }
      bool          IsPrimary() const { return fIsPrimary; }

      bool operator  <(gar::sdp::EnergyDeposit const& b) const;

#endif

    private:

      int    fTrackID;   ///< g4 track ID of particle making the deposit
      float fTime;      ///< time of the energy deposit
      float  fEnergy;    ///< energy deposited
      float  fX;         ///< x position of the energy deposit
      float  fY;         ///< y position of the energy deposit
      float  fZ;         ///< z position of the energy deposit
      float  fdX;        ///< the size of the step for this energy deposit
      bool   fIsPrimary; ///< is this from the primary particle, or a secondary (ie EM)
    };


  } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_EnergyDeposit_h */
