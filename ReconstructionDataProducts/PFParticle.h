//
//  PFParticle.h
//
//  Created by Eldwan Brianne on 07/06/2020.
//

#ifndef GAR_RECONSTRUCTIONDATAPRODUCTS_PFParticle_h
#define GAR_RECONSTRUCTIONDATAPRODUCTS_PFParticle_h

#include <iostream>

#include <map>
#include <list>
#include <vector>

#include "ReconstructionDataProducts/IDNumberGen.h"

namespace gar {
    namespace rec {

        class PFParticle {

        public:
            PFParticle();

            PFParticle(int type, float energy, float pos[3], float mom[3], float charge, int pdg, float goodness, size_t parent, std::vector<size_t> daughters);

            // let the compiler provide the dtor

            struct ParticleIDData {
                int fPdg;
                float fGoodness;
            };

        private:

            static gar::rec::IDNumber const FirstNumber = 800000;
            gar::rec::IDNumber fIDnumero;

            int                             fType{0};          ///< type of the PFParticle
            float                           fEnergy{0.};       ///< energy of the PFParticle in GeV
            float                           fPos[3];           ///< position of the PFParticle (cluster CoG)
            float                           fMom[3];           ///< momentum of the PFParticle in GeV
            float                           fMass{0};          ///< mass of the PFParticle
            float                           fCharge{0};        ///< charge of the PFParticle
            int                             fPdg;              ///< PDG of the PFParticle
            float                           fGoodness;         ///< Goodness of the PDG of the PFParticle
            size_t                          fParent{0};        ///< Parent of the PFParticle
            std::vector<size_t>             fDaughters;        ///< Daughters of the PFParticle

            #ifndef __GCCXML__

        public:

            bool operator==(const PFParticle& rhs) const;
            bool operator!=(const PFParticle& rhs) const;
            gar::rec::IDNumber getIDNumber() const;

            void setType(int type);
            void setEnergy(float energy);
            void setPosition(const float pos[3]);
            void setMomentum(const float mom[3]);
            void setMass(float mass);
            void setCharge(float charge);
            void setParticleID(int pdg, float goodness);
            void setParent(size_t parent);
            void setDaughters(std::vector<size_t> daughters);

            int                                  Type() const;
            float                                Energy() const;
            const float*                         Position() const;
            const float*                         Momentum() const;
            float                                Mass() const;
            float                                Charge() const;
            int                                  Pdg() const;
            float                                GoodnessOfPdg() const;
            size_t                               Parent() const;
            const std::vector<size_t>            Daughters() const;
            size_t                               NDaughters() const;

            friend std::ostream& operator << (std::ostream & o, gar::rec::PFParticle const& h);

            #endif

        };

        inline int                            gar::rec::PFParticle::Type()                const { return fType;         }
        inline float                          gar::rec::PFParticle::Energy()              const { return fEnergy;       }
        inline const float*                   gar::rec::PFParticle::Position()            const { return &fPos[0];      }
        inline const float*                   gar::rec::PFParticle::Momentum()            const { return &fMom[0];      }
        inline float                          gar::rec::PFParticle::Mass()                const { return fMass;         }
        inline float                          gar::rec::PFParticle::Charge()              const { return fCharge;       }
        inline int                            gar::rec::PFParticle::Pdg()                 const { return fPdg;          }
        inline float                          gar::rec::PFParticle::GoodnessOfPdg()       const { return fGoodness;     }
        inline size_t                         gar::rec::PFParticle::Parent()              const { return fParent;       }
        inline const std::vector<size_t>      gar::rec::PFParticle::Daughters()           const { return fDaughters;    }
        inline size_t                         gar::rec::PFParticle::NDaughters()          const { return fDaughters.size();    }
    }//rec
} // gar


#endif /* GAR_RECONSTRUCTIONDATAPRODUCTS_PFParticle_h */
