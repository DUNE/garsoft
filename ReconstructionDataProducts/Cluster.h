//
//  Cluster.h
//
//  Created by Eldwan Brianne on 08/29/18.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include <list>

#include "ReconstructionDataProducts/CaloHit.h"

namespace gar {
  namespace rec {

    class Cluster {

    public:
      Cluster();

      // let the compiler provide the dtor

      void AddToCluster(const gar::rec::CaloHit *p);

    private:

      float                                  fEnergy;      ///< energy of the calo hit in GeV
      std::list< const gar::rec::CaloHit* >  fCaloHitList;

#ifndef __GCCXML__

    public:

      Cluster(float energy, std::list< const gar::rec::CaloHit* > hitlist);

      float                           Energy()    const;
      const std::list< const gar::rec::CaloHit* > CaloHitList() const;

      friend std::ostream& operator << (std::ostream & o, gar::rec::Cluster const& h);

#endif

    };

    inline float                            gar::rec::Cluster::Energy()                   const { return fEnergy;      }
    inline const std::list< const gar::rec::CaloHit* >  gar::rec::Cluster::CaloHitList()              const { return fCaloHitList;      }
  } // rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_Cluster_h */
