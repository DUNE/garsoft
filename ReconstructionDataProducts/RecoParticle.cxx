#include "ReconstructionDataProducts/RecoParticle.h"

namespace gar {
    namespace rec {

        //--------------------------------------------------------------------------
        //Default constructor
        RecoParticle::RecoParticle(){
            IDNumberGen::create(FirstNumber);
            fIDnumber = IDNumberGen::create()->getNewOne();

            return;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setMomentum(float Momentum){
            fMomentum = Momentum;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setTotalCaloEnergy(float TotalCaloEnergy){
            fTotalCaloEnergy = TotalCaloEnergy;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setMeanCaloEnergy(float MeanCaloEnergy){
            fMeanCaloEnergy = MeanCaloEnergy;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setProtondEdxScore(float ProtondEdxScore){
            fProtondEdxScore = ProtondEdxScore;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setTotalECALEnergy(float TotalECALEnergy){
            fTotalECALEnergy = TotalECALEnergy;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setNHitsECAL(int NHitsECAL){
            fNHitsECAL = NHitsECAL;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setTotalMuIDEnergy(float TotalMuIDEnergy){
            fTotalMuIDEnergy = TotalMuIDEnergy;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setNHitsMuID(int NHitsMuID){
            fNHitsMuID = NHitsMuID;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setTrackEndECALed(int TrackEndECALed){
            fTrackEndECALed = TrackEndECALed;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setMuonScore(float MuonScore){
            fMuonScore = MuonScore;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setECALToFTime(float ECALToFTime){
            fECALToFTime = ECALToFTime;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setECALToFBeta(float ECALToFBeta){
            fECALToFBeta = ECALToFBeta;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setECALToFMass(float ECALToFMass){
            fECALToFMass = ECALToFMass;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setProtonToFScore(float ProtonToFScore){
            fProtonToFScore = ProtonToFScore;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setTrackEndVertexed(int TrackEndVertexed){
            fTrackEndVertexed = TrackEndVertexed;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setCharge(int Charge){
            fCharge = Charge;
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setStart(const float *Start){
            for (int i=0; i<3; ++i) fStart[i] = Start[i];
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setEnd(const float *End){
            for (int i=0; i<3; ++i) fEnd[i] = End[i];
        }

        //--------------------------------------------------------------------------
        void RecoParticle::setDirection(const float *Direction){
            for (int i=0; i<3; ++i) fDirection[i] = Direction[i];
        }
    }
}
