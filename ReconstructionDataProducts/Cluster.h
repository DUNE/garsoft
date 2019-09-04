//
//  Cluster.h
//
//  Created by Eldwan Brianne on 08/29/18.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include <map>
#include <list>

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/IDNumberGen.h"

namespace gar {
    namespace rec {

        class Cluster {

        public:
            Cluster();

            // let the compiler provide the dtor

        private:
            static gar::rec::IDNumber const FirstNumber = 300000;
            gar::rec::IDNumber fIDnumero;

            float                           fEnergy{0};      ///< energy of the calo hit in GeV
            float                           fPosition[3] = {0, 0, 0}; ///< position of the cluster in cm
            float                           fShape[6] = {0, 0, 0, 0, 0, 0}; ///< cluster shape parameters (Ellipsoid r1, r2, r3, vol, width)
            float                           fTheta{0}; ///< intrasic direction of the cluster theta
            float                           fPhi{0}; ///< intrasic direction of the cluster phi
            float                           fEigenVector[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0}; ///< EigenVectors of the cluster corresponding to the 3 main axis sorted in ascending order (main principal axis with smallest inertial of mass) normalised to length of 1
            int                             fParticleId{0}; ///< particle id flag
            std::vector<gar::rec::Track*>   fTracks{}; ///< vector of tracks associated to the cluster
            std::vector<gar::rec::CaloHit*> fHits{}; ///< vector of hit contribution
            std::vector<float>              fWeights{}; ///< vector of energy contribution of the hits

            #ifndef __GCCXML__

        public:

            //Copy constructor
            Cluster(const gar::rec::Cluster &) = default;

            bool operator==(const Cluster& rhs) const;
            bool operator!=(const Cluster& rhs) const;
            gar::rec::IDNumber getIDNumber() const;

            void addHit(gar::rec::CaloHit* hit, float contribution);
            void addTrack(gar::rec::Track* trk);
            void setEnergy(float energy);
            void setPosition(const float* position);
            void setITheta(float theta);
            void setIPhi(float phi);
            void setEigenVectors(const float* eigenvectors);
            void setShape(const float* shape);
            void setParticleID(int pid);

            const float                                 Energy() const;
            const float*                                Position() const;
            const float                                 ITheta() const;
            const float                                 IPhi() const;
            const float*                                EigenVectors() const;
            const float*                                Shape() const;
            const int                                   ParticleID() const;
            const std::vector<gar::rec::Track*>&        Tracks() const;
            const std::vector<gar::rec::CaloHit*>&      CalorimeterHits() const;
            const std::vector<float>&                   HitContributions() const;

            friend std::ostream& operator << (std::ostream & o, gar::rec::Cluster const& h);

            #endif

        };

        inline const float               gar::rec::Cluster::Energy()              const { return fEnergy;          }
        inline const float*              gar::rec::Cluster::Position()            const { return fPosition;        }
        inline const float               gar::rec::Cluster::ITheta()              const { return fTheta;          }
        inline const float               gar::rec::Cluster::IPhi()                const { return fPhi;          }
        inline const float*              gar::rec::Cluster::EigenVectors()        const { return fEigenVector;     }
        inline const float*              gar::rec::Cluster::Shape()               const { return fShape;           }
        inline const int                 gar::rec::Cluster::ParticleID()          const { return fParticleId;      }
        inline const std::vector<gar::rec::Track*>&   gar::rec::Cluster::Tracks() const { return fTracks;      }
        inline const std::vector<gar::rec::CaloHit*>& gar::rec::Cluster::CalorimeterHits()  const { return fHits;    }
        inline const std::vector<float>&              gar::rec::Cluster::HitContributions() const { return fWeights; }
    } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h */
