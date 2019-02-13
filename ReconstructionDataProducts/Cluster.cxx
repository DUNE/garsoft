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
        Cluster::Cluster()
        : fEnergy(0.),
        fInnerLayer(0.),
        fOuterLayer(0.),
        fnCaloHits(0)
        {
            fInitialDirection = CLHEP::Hep3Vector(0., 0., 0.);
            fDirection = CLHEP::Hep3Vector(0., 0., 0.);
            fGoG[0] = 0.;
            fGoG[1] = 0.;
            fGoG[2] = 0.;

            fOrderedCaloHitList.clear();
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

        void Cluster::AddToCluster(const gar::rec::CaloHit *pCaloHit)
        {
            this->ResetClusterParameters();

            const unsigned int layer(pCaloHit->GetLayer());
            this->AddtoOrderedList(pCaloHit, layer);

            this->UpdateClusterParameters();

            if (fInnerLayer == 0 || (layer < fInnerLayer))
            fInnerLayer = layer;

            if (fOuterLayer == 0 || (layer > fOuterLayer))
            fOuterLayer = layer;
        }

        void Cluster::ResetClusterParameters()
        {
            //Reset main cluster parameters
            fnCaloHits = 0;
            fEnergy = 0.;

            fInitialDirection.set(0., 0., 0.);
            fDirection.set(0., 0., 0.);
            fEigenVector[0].set(0., 0., 0.);
            fEigenVector[1].set(0., 0., 0.);
            fEigenVector[2].set(0., 0., 0.);

            fGoG[0] = 0.;
            fGoG[1] = 0.;
            fGoG[2] = 0.;
        }

        void Cluster::UpdateClusterParameters()
        {
            //Update energy, vector position, direction...
            float ePosSum[3] = {0., 0., 0.};

            for(OrderedCaloHitList::iterator it = fOrderedCaloHitList.begin(); it!= fOrderedCaloHitList.end(); ++it)
            {
                CaloHitList *const pCaloHitList(it->second);
                for(const CaloHit *const pCaloHit : *pCaloHitList)
                {
                    fEnergy += pCaloHit->Energy();

                    ePosSum[0] += pCaloHit->Energy() * pCaloHit->Position()[0];
                    ePosSum[1] += pCaloHit->Energy() * pCaloHit->Position()[1];
                    ePosSum[2] += pCaloHit->Energy() * pCaloHit->Position()[2];

                    ++fnCaloHits;
                }
            }

            fGoG[0] = ePosSum[0] / fEnergy;
            fGoG[1] = ePosSum[1] / fEnergy;
            fGoG[2] = ePosSum[2] / fEnergy;

            this->EigenVector();
        }

        void Cluster::EigenVector()
        {
            /* See https://github.com/iLCSoft/MarlinUtil/blob/master/source/src/ClusterShapes.cc */

            double aIne[3][3];
            float _ValAnalogInertia[3];
            float _VecAnalogInertia[9];

            for (int i(0); i < 3; ++i) {
                for (int j(0); j < 3; ++j) {
                    aIne[i][j] = 0.0;
                }
            }

            for(OrderedCaloHitList::iterator it = fOrderedCaloHitList.begin(); it!= fOrderedCaloHitList.end(); ++it)
            {
                CaloHitList *const pCaloHitList(it->second);
                for(const CaloHit *const pCaloHit : *pCaloHitList)
                {
                    float dX = pCaloHit->Position()[0] - fGoG[0];
                    float dY = pCaloHit->Position()[1] - fGoG[1];
                    float dZ = pCaloHit->Position()[2] - fGoG[2];

                    aIne[0][0] += pCaloHit->Energy() * (dY*dY+dZ*dZ);
                    aIne[1][1] += pCaloHit->Energy() * (dX*dX+dZ*dZ);
                    aIne[2][2] += pCaloHit->Energy() * (dX*dX+dY*dY);
                    aIne[0][1] -= pCaloHit->Energy() * dX*dY;
                    aIne[0][2] -= pCaloHit->Energy() * dX*dZ;
                    aIne[1][2] -= pCaloHit->Energy() * dY*dZ;
                }
            }

            for (int i(0); i < 2; ++i) {
                for (int j = i+1; j < 3; ++j) {
                    aIne[j][i] = aIne[i][j];
                }
            }

            //****************************************
            // analog Inertia
            //****************************************

            gsl_matrix_view aMatrix = gsl_matrix_view_array((double*)aIne, 3, 3);
            gsl_vector* aVector = gsl_vector_alloc(3);
            gsl_matrix* aEigenVec = gsl_matrix_alloc(3, 3);
            gsl_eigen_symmv_workspace* wa = gsl_eigen_symmv_alloc(3);
            gsl_eigen_symmv(&aMatrix.matrix, aVector, aEigenVec, wa);
            gsl_eigen_symmv_free(wa);
            gsl_eigen_symmv_sort(aVector, aEigenVec, GSL_EIGEN_SORT_ABS_ASC);

            for (int i(0); i < 3; i++) {
                _ValAnalogInertia[i] = gsl_vector_get(aVector, i);
                for (int j(0); j < 3; j++) {
                    _VecAnalogInertia[i+3*j] = gsl_matrix_get(aEigenVec, i, j);
                }
            }

            // Main principal points away from IP
            float radius = 0.;
            float radius2 = 0.;

            for (int i(0); i < 3; ++i) {
                radius += fGoG[i] * fGoG[i];
                radius2 += (fGoG[i] + _VecAnalogInertia[i]) * (fGoG[i] + _VecAnalogInertia[i]);
            }

            if ( radius2 < radius ) {
                for (int i(0); i < 3; ++i)
                _VecAnalogInertia[i] = - _VecAnalogInertia[i];
            }

            //Direction
            fDirection.set(_ValAnalogInertia[0], _ValAnalogInertia[1], _ValAnalogInertia[2]);
            //EigenVectors for the 3 cluster axis
            fEigenVector[0].set(_VecAnalogInertia[0], _VecAnalogInertia[1], _VecAnalogInertia[2]);
            fEigenVector[1].set(_VecAnalogInertia[3], _VecAnalogInertia[4], _VecAnalogInertia[5]);
            fEigenVector[2].set(_VecAnalogInertia[6], _VecAnalogInertia[7], _VecAnalogInertia[8]);

            gsl_vector_free(aVector);
            gsl_matrix_free(aEigenVec);
        }

        void Cluster::AddtoOrderedList(const gar::rec::CaloHit *const pCaloHit, const unsigned int layer)
        {
            OrderedCaloHitList::iterator iter = fOrderedCaloHitList.find(layer);

            if (fOrderedCaloHitList.end() == iter)
            {
                CaloHitList *const pCaloHitList = new CaloHitList;
                pCaloHitList->push_back(pCaloHit);

                if (!(fOrderedCaloHitList.insert(OrderedCaloHitList::value_type(layer, pCaloHitList)).second))
                {
                    delete pCaloHitList;
                    return;
                }
            }
            else
            {
                if (iter->second->end() != std::find(iter->second->begin(), iter->second->end(), pCaloHit))
                return;

                iter->second->push_back(pCaloHit);
            }
        }

    }
}
