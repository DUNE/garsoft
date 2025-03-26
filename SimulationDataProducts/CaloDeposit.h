//
//  CaloDeposit.h
//
//  Created by Eldwan Brianne on 8/23/17.
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h
#define GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h

#include <list>
#include <vector>

#include "RawDataProducts/CaloRawDigit.h"

namespace gar {
  namespace sdp {

    class CaloDeposit {
    public:
      CaloDeposit();

#ifndef __GCCXML__
      CaloDeposit(int trackID, float t, float e, double pos[3], raw::CellID_t CellID)
        : fTrackID(trackID), fTime(t), fEnergy(e), fCellID(CellID), fStepLength(0.)
      {
        fPos[0] = pos[0];
        fPos[1] = pos[1];
        fPos[2] = pos[2];
      }

      CaloDeposit(int trackID, float t, float e, double pos[3], raw::CellID_t CellID, float step_l)
        : fTrackID(trackID), fTime(t), fEnergy(e), fCellID(CellID), fStepLength(step_l)
      {
        fPos[0] = pos[0];
        fPos[1] = pos[1];
        fPos[2] = pos[2];
      }

      CaloDeposit(double pos[3]) : fTrackID(0), fTime(0.), fEnergy(0.), fCellID(0), fStepLength(0.)
      {
        fPos[0] = pos[0];
        fPos[1] = pos[1];
        fPos[2] = pos[2];
      }

      int const& TrackID() const { return fTrackID; }
      float const& Time() const { return fTime; }
      float const& Energy() const { return fEnergy; }
      double const& X() const { return fPos[0]; }
      double const& Y() const { return fPos[1]; }
      double const& Z() const { return fPos[2]; }
      raw::CellID_t const& CellID() const { return fCellID; }
      double const* Pos() const { return &fPos[0]; }
      float const& StepLength() const { return fStepLength; }
      bool operator<(gar::sdp::CaloDeposit const& b) const;

      void operator+=(gar::sdp::CaloDeposit const& b);

#endif

    private:
      int fTrackID;   ///< g4 track ID of particle making the deposit
      float fTime;    ///< time of the energy deposit
      float fEnergy;  ///< energy deposited
      double fPos[3]; ///< position of the energy deposit
      raw::CellID_t
        fCellID; ///< cellID encoded in 64 bits containing det_id, stave, module, layer, slice, cellX and cellY, use Helper to access the values
      float fStepLength; ///< step length
    };

  } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h */
