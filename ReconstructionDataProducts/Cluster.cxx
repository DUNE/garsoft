#include "ReconstructionDataProducts/Cluster.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        //Default constructor
        Cluster::Cluster(){
            // The default constructor is used e.g. by art::DataViewImpl::getByLabel
            // Make sure all Cluster objects are numbered, lest art deep-copy uninitialized
            // instances and then operator==() evaluates meaninglessly true.
            IDNumberGen::create(FirstNumber);
            fIDnumero = IDNumberGen::create()->getNewOne();
            return;
        }

        //--------------------------------------------------------------------------
        void Cluster::setEnergy(float energy ) {
            fEnergy = energy ;
        }

        //--------------------------------------------------------------------------
        void Cluster::setPosition(const float* position) {
            for(int i=0;i<3;i++) { fPosition[i] = position[i]; }
        }

        //--------------------------------------------------------------------------
        void Cluster::setITheta(float theta){
            fTheta = theta;
        }

        //--------------------------------------------------------------------------
        void Cluster::setIPhi(float phi){
            fPhi = phi;
        }

        //--------------------------------------------------------------------------
        void Cluster::setEigenVectors(const float* eigenvectors){
            for(int i=0;i<9;i++) { fEigenVector[i] = eigenvectors[i]; }
        }

        //--------------------------------------------------------------------------
        void Cluster::setShape(const float* shape) {
            for(int i=0;i<6;i++) { fShape[i] = shape[i]; }
        }

        //--------------------------------------------------------------------------
        void Cluster::setParticleID(int pid) {
            fParticleId = pid;
        }

        //--------------------------------------------------------------------------
        void Cluster::addHit(gar::rec::CaloHit* hit, float contribution) {
            fHits.push_back( hit ) ;
            fWeights.push_back( contribution ) ;
        }

        //--------------------------------------------------------------------------
        void Cluster::addTrack(gar::rec::Track* trk) {
            fTracks.push_back( trk ) ;
        }

        //--------------------------------------------------------------------------
        std::ostream& operator<< (std::ostream& o, gar::rec::Cluster const& h)
        {
            o << "Cluster "
            << "\n\tEnergy = "
            << h.Energy()
            << "\n\tPID = "
            << h.ParticleID()
            << "\n\tMain EigenVector = ("
            << h.EigenVectors()[0] << ", " <<  h.EigenVectors()[1] << ", " << h.EigenVectors()[2] << ")"
            << "\n\tPosition = "
            << h.Position()[0] << ", " <<  h.Position()[1] << ", " << h.Position()[2] << ")";

            return o;
        }



        //--------------------------------------------------------------------------
        // ID number methods
        bool Cluster::operator==(const Cluster& rhs) const {
            return (this->fIDnumero == rhs.fIDnumero);
        }
    
        bool Cluster::operator!=(const Cluster& rhs) const {
	        return (this->fIDnumero != rhs.fIDnumero);
        }
    
        IDNumberGen::IDNumber Cluster::getIDNumber() const {return fIDnumero;}



    }
}
