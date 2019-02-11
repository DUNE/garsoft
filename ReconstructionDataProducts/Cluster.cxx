#include "ReconstructionDataProducts/Cluster.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        Cluster::Cluster()
        : fEnergy(0.)
        {
            fCaloHitList.clear();
            return;
        }

        //--------------------------------------------------------------------------
        Cluster::Cluster(float energy, std::list< const gar::rec::CaloHit* > hitlist)
        : fEnergy(0.),
        fCaloHitList(hitlist.begin(), hitlist.end())
        {
            return;
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::Cluster const& h)
        {

            o << "Cluster "
            << "\n\tenergy = "
            << h.Energy();

            return o;
        }

        void Cluster::AddToCluster(const gar::rec::CaloHit *p)
        {
            fCaloHitList.push_back(p);
        }

    }
}
