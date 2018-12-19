//
//  CaloDeposit.h
//
//  Created by Eldwan Brianne on 8/23/17.
//

#ifndef GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h
#define GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h

#include <list>
#include <vector>

namespace gar {
    namespace sdp {


        class CaloDeposit{
        public:

            CaloDeposit();

            #ifndef __GCCXML__
            CaloDeposit(int trackID,
            float t,
            float e,
            double pos[3],
            long long int CellID)
            : fTrackID  (trackID)
            , fTime     (t)
            , fEnergy   (e)
            , fCellID   (CellID)
            {
                fPos[0] = pos[0];
                fPos[1] = pos[1];
                fPos[2] = pos[2];

                fContrib.clear();
            }

            CaloDeposit(double pos[3])
            : fTrackID(0),
            fTime(0.),
            fEnergy(0.),
            fCellID(0)
            {
                fPos[0] = pos[0];
                fPos[1] = pos[1];
                fPos[2] = pos[2];

                fContrib.clear();
            }

            int    const& TrackID()   const { return fTrackID;   }
            float  const& Time()      const { return fTime;      }
            float  const& Energy()    const { return fEnergy;    }
            double  const& X()         const { return fPos[0];         }
            double  const& Y()         const { return fPos[1];         }
            double  const& Z()         const { return fPos[2];         }
            long long int  const& CellID()      const { return fCellID;    }
            const double* Pos() const { return fPos; }
            std::vector<gar::sdp::CaloDeposit*> const& Contrib() const { return fContrib; }

            bool operator  <(gar::sdp::CaloDeposit const& b) const;

            #endif

            void setCellID(float cID) { fCellID = cID; }
            void setTrackID(float trkID) { fTrackID = trkID; }
            void setTime(float time) { fTime = time; }

            void AddEnergy(float energy) { fEnergy += energy; }
            void AddContrib(gar::sdp::CaloDeposit *contrib) { fContrib.push_back(contrib); }

        private:

            int fTrackID;   ///< g4 track ID of particle making the deposit
            float fTime;      ///< time of the energy deposit
            float fEnergy;    ///< energy deposited
            double fPos[3]; ///< position of the energy deposit
            long long int fCellID; ///< cellID encoded in 64 bits containing det_id, stave, module, layer, slice, cellX and cellY, use Helper to access the values
            std::vector<gar::sdp::CaloDeposit*> fContrib; ///<Vector of contributions
        };


    } // sdp
} // gar

#endif /* GAR_SIMULATIONDATAPRODUCTS_CaloDeposit_h */
