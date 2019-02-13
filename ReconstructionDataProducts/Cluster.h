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

namespace gar {
    namespace rec {

        typedef std::list< const gar::rec::CaloHit* > CaloHitList;
        typedef std::map< unsigned int, CaloHitList* > OrderedCaloHitList;

        class Cluster {

        public:
            Cluster();

            // let the compiler provide the dtor

            void AddToCluster(const gar::rec::CaloHit *p);

        private:

            void ResetClusterParameters();
            void UpdateClusterParameters();
            void EigenVector();
            void AddtoOrderedList(const gar::rec::CaloHit *const pCaloHit, const unsigned int layer);

            float                                  fEnergy;      ///< energy of the calo hit in GeV
            CLHEP::Hep3Vector                      fInitialDirection; ///< initial direction of the cluster (can be used if track seed used) in cm
            CLHEP::Hep3Vector                      fDirection; ///< main direction of the cluster in cm
            OrderedCaloHitList                     fOrderedCaloHitList; ///< list of hits that make the cluster
            unsigned int                           fInnerLayer; ///< first layer of the cluster
            unsigned int                           fOuterLayer; ///< last layer of the cluster
            unsigned int                           fnCaloHits; ///< number of hits in the cluster
            float                                  fGoG[3]; ///< coordinates of the center of gravity of the cluster in cm
            CLHEP::Hep3Vector                      fEigenVector[3]; ///< EigenVectors of the cluster corresponding to the 3 main axis sorted in ascending order (main principal axis with smallest inertial of mass) normalised to length of 1

            #ifndef __GCCXML__

        public:

            const float                                 Energy()      const;
            const OrderedCaloHitList                    &OrderedCaloHitList() const;
            const CLHEP::Hep3Vector                     InitialDirection() const;
            const CLHEP::Hep3Vector                     Direction() const;
            const unsigned int                          NCaloHits() const;
            const float*                                CenterOfGravity()      const;
            const CLHEP::Hep3Vector*                    EigenVectors() const;

            friend std::ostream& operator << (std::ostream & o, gar::rec::Cluster const& h);

            #endif

        };

        inline const float                                  gar::rec::Cluster::Energy()       const { return fEnergy;      }
        inline const OrderedCaloHitList &gar::rec::Cluster::OrderedCaloHitList()              const { return fOrderedCaloHitList;      }
        inline const CLHEP::Hep3Vector  gar::rec::Cluster::InitialDirection()                 const { return fInitialDirection;      }
        inline const CLHEP::Hep3Vector  gar::rec::Cluster::Direction()                        const { return fDirection;      }
        inline const unsigned int       gar::rec::Cluster::NCaloHits()                        const { return fnCaloHits;      }
        inline const float*             gar::rec::Cluster::CenterOfGravity()                  const { return &fGoG[0];      }
        inline const CLHEP::Hep3Vector*  gar::rec::Cluster::EigenVectors()                    const { return &fEigenVector[0];      }
    } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h */
