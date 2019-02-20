#include "RecoAlg/Cluster.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>

namespace gar {
    namespace rec {
        namespace alg {

            //--------------------------------------------------------------------------
            Cluster::Cluster()
            : m_Energy(0.),
            m_pTrackSeed(nullptr),
            m_InnerLayer(0.),
            m_OuterLayer(0.),
            m_nCaloHits(0),
            m_particleId(0)
            {
                m_InitialDirection = CLHEP::Hep3Vector(0., 0., 0.);
                m_Direction = CLHEP::Hep3Vector(0., 0., 0.);
                m_GoG[0] = 0.;
                m_GoG[1] = 0.;
                m_GoG[2] = 0.;
                m_eGoG[0] = 0.;
                m_eGoG[1] = 0.;
                m_eGoG[2] = 0.;

                m_OrderedCaloHitList.clear();
                m_TrackList.clear();
                return;
            }

            //--------------------------------------------------------------------------
            Cluster::~Cluster()
            {

            }

            //--------------------------------------------------------------------------
            void Cluster::AddToCluster(const gar::rec::CaloHit *pCaloHit)
            {
                this->ResetClusterParameters();

                const unsigned int layer(pCaloHit->GetLayer());
                this->AddtoOrderedList(pCaloHit, layer);

                const float x(pCaloHit->GetPositionVector().x());
                const float y(pCaloHit->GetPositionVector().y());
                const float z(pCaloHit->GetPositionVector().z());

                OrderedCaloHitList::const_iterator iter = m_OrderedCaloHitList.find(layer);
                if ((m_OrderedCaloHitList.end() != iter) && (iter->second->size() > 1))
                {
                    SimplePoint &mypoint = m_sumXYZByLayer[layer];
                    mypoint.m_xyzPositionSums[0] += x;
                    mypoint.m_xyzPositionSums[1] += y;
                    mypoint.m_xyzPositionSums[2] += z;
                    ++mypoint.m_nHits;
                }
                else
                {
                    SimplePoint &mypoint = m_sumXYZByLayer[layer];
                    mypoint.m_xyzPositionSums[0] = x;
                    mypoint.m_xyzPositionSums[1] = y;
                    mypoint.m_xyzPositionSums[2] = z;
                    mypoint.m_nHits = 1;
                }

                this->UpdateClusterParameters();

                if (m_InnerLayer == 0 || (layer < m_InnerLayer))
                m_InnerLayer = layer;

                if (m_OuterLayer == 0 || (layer > m_OuterLayer))
                m_OuterLayer = layer;
            }

            //--------------------------------------------------------------------------
            void Cluster::AddToCluster(const gar::rec::Track *t)
            {
                if(nullptr == m_pTrackSeed)
                m_pTrackSeed = t;

                //push back the track into the list of associated tracks
                m_TrackList.push_back(t);

                //Get track state at the calorimeter and use it to set the initial direction
                m_InitialDirection = CLHEP::Hep3Vector(0., 0., 0.);
            }

            //--------------------------------------------------------------------------
            void Cluster::ResetClusterParameters()
            {
                //Reset main cluster parameters
                m_nCaloHits = 0;
                m_Energy = 0.;

                m_InitialDirection.set(0., 0., 0.);
                m_Direction.set(0., 0., 0.);
                m_EigenVector[0].set(0., 0., 0.);
                m_EigenVector[1].set(0., 0., 0.);
                m_EigenVector[2].set(0., 0., 0.);

                m_GoG[0] = 0.;
                m_GoG[1] = 0.;
                m_GoG[2] = 0.;

                m_eGoG[0] = 0.;
                m_eGoG[1] = 0.;
                m_eGoG[2] = 0.;
            }

            //--------------------------------------------------------------------------
            void Cluster::UpdateClusterParameters()
            {
                //Update energy, vector position, direction...
                float ePosSum[3] = {0., 0., 0.};
                float PosSum[3] = {0., 0., 0.};

                for(OrderedCaloHitList::iterator it = m_OrderedCaloHitList.begin(); it!= m_OrderedCaloHitList.end(); ++it)
                {
                    CaloHitList *const pCaloHitList(it->second);
                    for(const CaloHit *const pCaloHit : *pCaloHitList)
                    {
                        m_Energy += pCaloHit->Energy();

                        PosSum[0] += pCaloHit->Position()[0];
                        PosSum[1] += pCaloHit->Position()[1];
                        PosSum[2] += pCaloHit->Position()[2];

                        ePosSum[0] += pCaloHit->Energy() * pCaloHit->Position()[0];
                        ePosSum[1] += pCaloHit->Energy() * pCaloHit->Position()[1];
                        ePosSum[2] += pCaloHit->Energy() * pCaloHit->Position()[2];

                        ++m_nCaloHits;
                    }
                }

                //Center of gravity
                m_GoG[0] = PosSum[0] / m_nCaloHits;
                m_GoG[1] = PosSum[1] / m_nCaloHits;
                m_GoG[2] = PosSum[2] / m_nCaloHits;

                //Energy weighted center of gravity
                m_eGoG[0] = ePosSum[0] / m_Energy;
                m_eGoG[1] = ePosSum[1] / m_Energy;
                m_eGoG[2] = ePosSum[2] / m_Energy;

                this->EigenVector();
            }

            //--------------------------------------------------------------------------
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

                for(OrderedCaloHitList::iterator it = m_OrderedCaloHitList.begin(); it!= m_OrderedCaloHitList.end(); ++it)
                {
                    CaloHitList *const pCaloHitList(it->second);
                    for(const CaloHit *const pCaloHit : *pCaloHitList)
                    {
                        float dX = pCaloHit->Position()[0] - m_eGoG[0];
                        float dY = pCaloHit->Position()[1] - m_eGoG[1];
                        float dZ = pCaloHit->Position()[2] - m_eGoG[2];

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
                    radius += m_eGoG[i] * m_eGoG[i];
                    radius2 += (m_eGoG[i] + _VecAnalogInertia[i]) * (m_eGoG[i] + _VecAnalogInertia[i]);
                }

                if ( radius2 < radius ) {
                    for (int i(0); i < 3; ++i)
                    _VecAnalogInertia[i] = - _VecAnalogInertia[i];
                }

                //Direction
                m_Direction.set(_ValAnalogInertia[0], _ValAnalogInertia[1], _ValAnalogInertia[2]);
                //EigenVectors for the 3 cluster axis
                m_EigenVector[0].set(_VecAnalogInertia[0], _VecAnalogInertia[1], _VecAnalogInertia[2]);
                m_EigenVector[1].set(_VecAnalogInertia[3], _VecAnalogInertia[4], _VecAnalogInertia[5]);
                m_EigenVector[2].set(_VecAnalogInertia[6], _VecAnalogInertia[7], _VecAnalogInertia[8]);

                gsl_vector_free(aVector);
                gsl_matrix_free(aEigenVec);
            }

            //--------------------------------------------------------------------------
            void Cluster::AddtoOrderedList(const gar::rec::CaloHit *const pCaloHit, const unsigned int layer)
            {
                OrderedCaloHitList::iterator iter = m_OrderedCaloHitList.find(layer);

                if (m_OrderedCaloHitList.end() == iter)
                {
                    CaloHitList *const pCaloHitList = new CaloHitList;
                    pCaloHitList->push_back(pCaloHit);

                    if (!(m_OrderedCaloHitList.insert(OrderedCaloHitList::value_type(layer, pCaloHitList)).second))
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

            //--------------------------------------------------------------------------
            const CLHEP::Hep3Vector Cluster::GetCentroid(const unsigned int layer) const
            {
                std::map<unsigned int, SimplePoint>::const_iterator pointValueIter = m_sumXYZByLayer.find(layer);

                if (m_sumXYZByLayer.end() == pointValueIter)
                return CLHEP::Hep3Vector(0., 0., 0.);

                const SimplePoint &mypoint = pointValueIter->second;

                if (0 == mypoint.m_nHits)
                return CLHEP::Hep3Vector(0., 0., 0.);

                return CLHEP::Hep3Vector(static_cast<float>(mypoint.m_xyzPositionSums[0] / static_cast<float>(mypoint.m_nHits)),
                static_cast<float>(mypoint.m_xyzPositionSums[1] / static_cast<float>(mypoint.m_nHits)),
                static_cast<float>(mypoint.m_xyzPositionSums[2] / static_cast<float>(mypoint.m_nHits)));
            }

        }
    }
}
