#include "ReconstructionDataProducts/Cluster.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        //Default constructor
        Cluster::Cluster()
        : fEnergy(0.),
        fInnerLayer(0.),
        fOuterLayer(0.),
        fnCaloHits(0),
        fparticleId(0)
        {
            fDirection = CLHEP::Hep3Vector(0., 0., 0.);
            fGoG[0] = 0.;
            fGoG[1] = 0.;
            fGoG[2] = 0.;

            return;
        }

        //--------------------------------------------------------------------------
        //Constructor
        Cluster::Cluster(const float energy,
        const CLHEP::Hep3Vector direction,
        const unsigned int innerlayer,
        const unsigned int outerlayer,
        const unsigned int nCaloHits,
        const float* cog,
        const CLHEP::Hep3Vector* eigenvector,
        const int pid)
        : fEnergy(energy),
        fDirection(direction),
        fInnerLayer(innerlayer),
        fOuterLayer(outerlayer),
        fnCaloHits(nCaloHits),
        fparticleId(pid)
        {
            for(unsigned int i = 0; i < 3; i++)
            {
                fGoG[i] = cog[i];
                fEigenVector[i] = eigenvector[i];
            }
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::Cluster const& h)
        {
            o << "Cluster "
            << "\n\tEnergy = "
            << h.Energy()
            << "\n\tnHits = "
            << h.NCaloHits()
            << "\n\tPID = "
            << h.ParticleID()
            << "\n\tDirection = ("
            << h.Direction().x() << ", " <<  h.Direction().y() << ", " << h.Direction().z() << ")"
            << "\n\tCoG = "
            << h.CenterOfGravity()[0] << ", " <<  h.CenterOfGravity()[1] << ", " << h.CenterOfGravity()[2] << ")";

            return o;
        }

    }
}
