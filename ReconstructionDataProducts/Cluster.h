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

namespace gar {
    namespace rec {

        class Cluster {

        public:
            Cluster();

            // let the compiler provide the dtor

        private:

            float                                  fEnergy;      ///< energy of the calo hit in GeV
            CLHEP::Hep3Vector                      fDirection; ///< main direction of the cluster in cm
            unsigned int                           fInnerLayer; ///< first layer of the cluster
            unsigned int                           fOuterLayer; ///< last layer of the cluster
            unsigned int                           fnCaloHits; ///< number of hits in the cluster
            float                                  fGoG[3]; ///< coordinates of the center of gravity of the cluster in cm
            CLHEP::Hep3Vector                      fEigenVector[3]; ///< EigenVectors of the cluster corresponding to the 3 main axis sorted in ascending order (main principal axis with smallest inertial of mass) normalised to length of 1
            int                                    fparticleId; ///< particle id flag

            #ifndef __GCCXML__

        public:

            Cluster(const float energy,
            const CLHEP::Hep3Vector direction,
            const unsigned int innerlayer,
            const unsigned int outerlayer,
            const unsigned int nCaloHits,
            const float* cog,
            const CLHEP::Hep3Vector* eigenvector,
            const int pid);

            const float                                 Energy() const;
            const CLHEP::Hep3Vector                     Direction() const;
            const unsigned int                          InnerLayer() const;
            const unsigned int                          OuterLayer() const;
            const unsigned int                          NCaloHits() const;
            const float*                                CenterOfGravity() const;
            const CLHEP::Hep3Vector*                    EigenVectors() const;
            const int                                   ParticleID() const;

            friend std::ostream& operator << (std::ostream & o, gar::rec::Cluster const& h);

            #endif

        };

        inline const float               gar::rec::Cluster::Energy()                           const { return fEnergy;      }
        inline const CLHEP::Hep3Vector   gar::rec::Cluster::Direction()                        const { return fDirection;      }
        inline const unsigned int        gar::rec::Cluster::InnerLayer()                       const { return fInnerLayer;      }
        inline const unsigned int        gar::rec::Cluster::OuterLayer()                       const { return fOuterLayer;      }
        inline const unsigned int        gar::rec::Cluster::NCaloHits()                        const { return fnCaloHits;      }
        inline const float*              gar::rec::Cluster::CenterOfGravity()                  const { return &fGoG[0];      }
        inline const CLHEP::Hep3Vector*  gar::rec::Cluster::EigenVectors()                     const { return &fEigenVector[0];      }
        inline const int                 gar::rec::Cluster::ParticleID()                       const { return fparticleId;      }
    } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h */
