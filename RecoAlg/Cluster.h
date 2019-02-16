#ifndef GAR_RECOALG_Cluster_h
#define GAR_RECOALG_Cluster_h

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

#include <map>
#include <list>

#include "ReconstructionDataProducts/CaloHit.h"
#include "ReconstructionDataProducts/Track.h"

namespace gar {
    namespace rec {
        namespace alg{

            typedef std::list< const gar::rec::Track* > TrackList;
            typedef std::list< const gar::rec::CaloHit* > CaloHitList;
            typedef std::map< unsigned int, CaloHitList* > OrderedCaloHitList;

            class Cluster {

            public:
                Cluster();

                virtual ~Cluster();

                void AddToCluster(const gar::rec::CaloHit *p);

                void AddToCluster(const gar::rec::Track *t);

                float getEnergy() const;

                unsigned int getNCaloHits() const;

                const OrderedCaloHitList &getOrderedCaloHitList() const;

                const TrackList &getTrackList() const;

                const CLHEP::Hep3Vector getInitialDirection() const;

                const CLHEP::Hep3Vector getDirection() const;

                const float* getCenterOfGravity() const;

                const float* getEnergyWeightedCenterOfGravity() const;

                const CLHEP::Hep3Vector* getEigenVectors() const;

                const gar::rec::Track* getTrackSeed() const;

                bool isTrackSeeded() const;

                unsigned int getInnerLayer() const;

                unsigned int getOuterLayer() const;

                int getParticleId() const;

                const CLHEP::Hep3Vector GetCentroid(const unsigned int layer) const;

            private:

                void ResetClusterParameters();

                void UpdateClusterParameters();

                void EigenVector();

                void AddtoOrderedList(const gar::rec::CaloHit *const pCaloHit, const unsigned int layer);

                //Defines a point used for each layer
                class SimplePoint
                {
                public:
                    double                  m_xyzPositionSums[3];           ///< The sum of the x, y and z hit positions in the pseudo layer
                    unsigned int            m_nHits;                        ///< The number of hits in the pseudo layer
                };

                std::map<unsigned int, SimplePoint>    m_sumXYZByLayer; ///< sum of the position by layer to simplify centroid
                float                                  m_Energy;      ///< energy of the calo hit in GeV
                CLHEP::Hep3Vector                      m_InitialDirection; ///< initial direction of the cluster (can be used if track seed used) in cm
                CLHEP::Hep3Vector                      m_Direction; ///< main direction of the cluster in cm
                OrderedCaloHitList                     m_OrderedCaloHitList; ///< list of hits that make the cluster
                const gar::rec::Track*                 m_pTrackSeed; ///< Address of the track used as seed for the cluster
                unsigned int                           m_InnerLayer; ///< first layer of the cluster
                unsigned int                           m_OuterLayer; ///< last layer of the cluster
                unsigned int                           m_nCaloHits; ///< number of hits in the cluster
                float                                  m_GoG[3]; ///< coordinates of the center of gravity of the cluster in cm
                float                                  m_eGoG[3]; ///< coordinates of the energy weighted center of gravity of the cluster in cm
                CLHEP::Hep3Vector                      m_EigenVector[3]; ///< EigenVectors of the cluster corresponding to the 3 main axis sorted in ascending order (main principal axis with smallest inertial of mass) normalised to length of 1
                int                                    m_particleId; ///< Particle id flag
                TrackList                              m_TrackList; ///< list of associated tracks

            };

            //---------------------------------------------------------------------------
            inline float Cluster::getEnergy() const
            {
                return m_Energy;
            }

            //---------------------------------------------------------------------------
            inline const OrderedCaloHitList &Cluster::getOrderedCaloHitList() const
            {
                return m_OrderedCaloHitList;
            }

            //---------------------------------------------------------------------------
            inline const TrackList &Cluster::getTrackList() const
            {
                return m_TrackList;
            }

            //---------------------------------------------------------------------------
            inline const CLHEP::Hep3Vector Cluster::getInitialDirection() const
            {
                return m_InitialDirection;
            }

            //---------------------------------------------------------------------------
            inline const CLHEP::Hep3Vector Cluster::getDirection() const
            {
                return m_Direction;
            }

            //---------------------------------------------------------------------------
            inline unsigned int Cluster::getNCaloHits() const
            {
                return m_nCaloHits;
            }

            //---------------------------------------------------------------------------
            inline const float* Cluster::getCenterOfGravity() const
            {
                return &m_GoG[0];
            }

            //---------------------------------------------------------------------------
            inline const float* Cluster::getEnergyWeightedCenterOfGravity() const
            {
                return &m_eGoG[0];
            }

            //---------------------------------------------------------------------------
            inline const CLHEP::Hep3Vector* Cluster::getEigenVectors() const
            {
                return &m_EigenVector[0];
            }

            //---------------------------------------------------------------------------
            inline const gar::rec::Track* Cluster::getTrackSeed() const
            {
                return m_pTrackSeed;
            }

            //----------------------------------------------------------------------------
            inline bool Cluster::isTrackSeeded() const
            {
                return (nullptr != m_pTrackSeed);
            }

            //----------------------------------------------------------------------------
            inline unsigned int Cluster::getInnerLayer() const
            {
                return m_InnerLayer;
            }

            //----------------------------------------------------------------------------
            inline unsigned int Cluster::getOuterLayer() const
            {
                return m_OuterLayer;
            }

            //----------------------------------------------------------------------------
            inline int Cluster::getParticleId() const
            {
                return m_particleId;
            }

        }
    }
}


#endif
