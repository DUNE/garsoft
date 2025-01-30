//
//  RecoParticle.h
//
//  Created by Francisco Martinez Lopez on 28/02/2024.
//  f.martinezlopez@qmul.ac.uk
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_RecoParticle_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_RecoParticle_h

#include <iostream>

#include "IDNumberGen.h"

namespace gar {
    namespace rec {

        class RecoParticle {

            public:
            RecoParticle();

            inline gar::rec::IDNumber getIDNumber()      const { return fIDnumber;         };

            inline float              Momentum()         const { return fMomentum;         };

            inline float              TotalCaloEnergy()  const { return fTotalCaloEnergy;  };
            inline float              MeanCaloEnergy()   const { return fMeanCaloEnergy;   };

            inline float              ProtondEdxScore()  const { return fProtondEdxScore;  };

            inline float              TotalECALEnergy()  const { return fTotalECALEnergy;  };
            inline int                NHitsECAL()        const { return fNHitsECAL;        };
            inline float              TotalMuIDEnergy()  const { return fTotalMuIDEnergy;  };
            inline int                NHitsMuID()        const { return fNHitsMuID;        };

            inline int                TrackEndECALed()   const { return fTrackEndECALed;   };

            inline float              MuonScore()        const { return fMuonScore;        };

            inline float              ECALToFTime()      const { return fECALToFTime;      };
            inline float              ECALToFBeta()      const { return fECALToFBeta;      };
            inline float              ECALToFMass()      const { return fECALToFMass;      };

            inline float              ProtonToFScore()   const { return fProtonToFScore;   };

            inline int                TrackEndVertexed() const { return fTrackEndVertexed; };

            inline int                Charge()           const { return fCharge;           };

            inline const float*       Start()            const { return fStart;            };
            inline const float*       End()              const { return fEnd;              };

            inline const float*       Direction()        const { return fDirection;        };

            // let the compiler provide the dtor

            private:

            static gar::rec::IDNumber const FirstNumber = 900000;
            gar::rec::IDNumber fIDnumber;

            float fMomentum = 0.0; ///< momentum from curvature in the B field [GeV/c]

            float fTotalCaloEnergy = 0.0; ///< total energy deposited in the TPC [MeV]
            float fMeanCaloEnergy  = 0.0; ///< mean dE/dx in TPC [keV/cm]

            float fProtondEdxScore = 0.0;

            float fTotalECALEnergy = 0.0;
            int   fNHitsECAL       = 0;
            float fTotalMuIDEnergy = 0.0;
            int   fNHitsMuID       = 0;

            int   fTrackEndECALed = -1;

            float fMuonScore = 0.0;

            float fECALToFTime = -1.0;  ///< arrival time to ECal [ns]
            float fECALToFBeta = -1.0;  ///< velocity from ToF and length
            float fECALToFMass = -1.0;  ///< mass from velocity and momentum

            float fProtonToFScore = 0.0;

            int   fTrackEndVertexed = -1;

            int   fCharge = 0;

            float fStart[3];    ///< particle start position [cm]
            float fEnd[3];      ///< particle end position [cm]

            float fDirection[3];      ///< particle direction


#ifndef __GCCXML__

            public:

            void setMomentum(float Momentum);

            void setTotalCaloEnergy(float TotalCaloEnergy);
            void setMeanCaloEnergy(float MeanCaloEnergy);

            void setProtondEdxScore(float ProtondEdxScore);

            void setTotalECALEnergy(float TotalECALEnergy);
            void setNHitsECAL(int NHitsECAL);
            void setTotalMuIDEnergy(float TotalMuIDEnergy);
            void setNHitsMuID(int NHitsMuID);

            void setTrackEndECALed(int TrackEndECALed);

            void setMuonScore(float MuonScore);

            void setECALToFTime(float ECALToFTime);
            void setECALToFBeta(float ECALToFBeta);
            void setECALToFMass(float ECALToFMass);

            void setProtonToFScore(float ProtonToFScore);

            void setTrackEndVertexed(int TrackEndVertexed);

            void setCharge(int Charge);

            void setStart(const float *Start);
            void setEnd(const float *End);

            void setDirection(const float *Direction);

#endif

        };

    } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_RecoParticle_h */
