////////////////////////////////////////////////////////////////////////
// Class:       tpctrackfit2
// Plugin Type: producer (art v3_00_00)
// File:        tpctrackfit2_module.cc
//
// Generated at Tue Feb  5 11:34:54 2019 by Thomas Junk using cetskelgen
// from cetlib version v3_04_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include <memory>
#include <set>

// ROOT includes

#include "TVectorF.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"

// GArSoft Includes
#include "DetectorInfo/GArMagneticField.h"
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "ReconstructionDataProducts/TrackIoniz.h"
#include "ReconstructionDataProducts/TrackTrajectory.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "Geometry/GeometryGAr.h"
#include "CoreUtils/ServiceUtil.h"
#include "RecoAlg/fastSimulation.h"
#include "RecoAlg/garutils.h"
#include "nug4/MagneticFieldServices/MagneticFieldService.h"
#include "Geant4/G4ThreeVector.hh"

namespace gar {
  namespace rec {

    class tpctrackfit2 : public art::EDProducer {
    public:
      explicit tpctrackfit2(fhicl::ParameterSet const& p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tpctrackfit2(tpctrackfit2 const&) = delete;
      tpctrackfit2(tpctrackfit2&&) = delete;
      tpctrackfit2& operator=(tpctrackfit2 const&) = delete;
      tpctrackfit2& operator=(tpctrackfit2&&) = delete;

      // Required functions.
      void produce(art::Event& e) override;

    private:

      // Declare member data here.

      std::string fPatRecLabel;            ///< input patrec tracks and associations
      int fPrintLevel;                     ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      float  fKalCurvStepUncSq;            ///< constant uncertainty term on each step of the Kalman fit -- squared, for curvature
      float  fKalPhiStepUncSq;             ///< constant uncertainty term on each step of the Kalman fit -- squared, for phi
      float  fKalLambdaStepUncSq;          ///< constant uncertainty term on each step of the Kalman fit -- squared, for lambda
      float  fKalCovZYMeasure;             ///< constant uncertainty term on measurement in Kalman (the R matrix)
      float  fTPCClusterResolYZ;           ///< pad size in cm in YZ to determine step size
      float  fTPCClusterResolX;            ///< drift direction contribution to determine step size (resolution of a TPCCluster)
      unsigned int fInitialTPNTPCClusters; ///< number of TPCClusters to use for initial trackpar estimate, if present
      unsigned int fMinNumTPCClusters;     ///< minimum number of TPCClusters to define a track
      int fDumpTracks;                     ///< 0: do not print out tracks, 1: print out tracks
      float fRoadYZinFit;                  ///< cut in cm for dropping TPCClusters from tracks in fit
      int   fSortAlg;                      ///< which hit sorting alg to use.  1: old, 2: greedy distance sort
      float fSortDistCut;                  ///< distance cut to pass to hit sorter #2
      float fSortTransWeight;              ///< for use in hit sorting algorithm #1 -- transverse distance weight factor
      float fSortDistBack;                 ///< for use in hit sorting algorithm #1 -- how far to go back before raising the distance figure of merit
      float fMinIonizGapCut;               ///< Don't compute dEdx for this dx or larger

      float fTPCClusterResid__CROC_b;      ///< parameters to estimate residuals in YZ plane
      float fTPCClusterResid__CROC_m;
      float fTPCClusterResid__IROC_b;
      float fTPCClusterResid__IROC_m;
      float fTPCClusterResid_IOROC_b;
      float fTPCClusterResid_IOROC_m;
      float fTPCClusterResid_OOROC_b;
      float fTPCClusterResid_OOROC_m;



      int KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
                            std::vector<TrackPar> &trackpars,  std::vector<TrackIoniz> &trackions, std::vector<TrackTrajectory> &tracktrajs);
      
      double CalculateChi2(RVec<AliExternalTrackParam4D> ParamMC, RVec<AliExternalTrackParam4D> ParamIn, double xstart, double ystart, double * GArCenter);

      void CalculateTrackInfo(float & length, std::vector<std::pair<float,float>>& dSigdX, std::vector<TVector3>& trajpts,
                                          std::vector<int> TPCClusterList, std::vector<TPCCluster> TPCClusters,
                                          RVec<AliExternalTrackParam4D> ParamIn, double xstart, double ystart, double * GArCenter);
      
      art::ServiceHandle<geo::GeometryGAr> euclid;

    };


    tpctrackfit2::tpctrackfit2(fhicl::ParameterSet const& p) : EDProducer{p}
    {
      // Call appropriate produces<>() functions here.
      // Call appropriate consumes<>() for any products to be retrieved by this module.

      fPatRecLabel       = p.get<std::string>("PatRecLabel","patrec");
      fPrintLevel        = p.get<int>("PrintLevel",0);
      fMinNumTPCClusters        = p.get<unsigned int>("MinNumTPCClusters",20);
      fKalCurvStepUncSq  = p.get<float>("KalCurvStepUncSq",1.0E-9);
      fKalPhiStepUncSq   = p.get<float>("KalPhiStepUncSq",1.0E-9);
      fKalLambdaStepUncSq = p.get<float>("KalLambdaStepUncSq",1.0E-9);
      fKalCovZYMeasure   = p.get<float>("KalCovZYMeasure", 4.0);
      fInitialTPNTPCClusters    = p.get<unsigned int>("InitialTPNTPCClusters",100);
      fDumpTracks        = p.get<int>("DumpTracks",0);
      fTPCClusterResolYZ        = p.get<float>("TPCClusterResolYZ",1.0); // TODO -- think about what this value is
      fTPCClusterResolX         = p.get<float>("TPCClusterResolX",0.5);  // this is probably much better
      fRoadYZinFit       = p.get<float>("RoadYZinFit",1.0);
      fSortTransWeight   = p.get<float>("SortTransWeight",0.1);
      fSortDistBack      = p.get<float>("SortDistBack",2.0);
      fMinIonizGapCut    = p.get<float>("MinIonizGapCut",5.0);
      fSortDistCut       = p.get<float>("SortDistCut",10.0);
      fSortAlg           = p.get<int>("SortAlg",2);

      fTPCClusterResid__CROC_b = p.get<float>("TPCClusterResid__CROC_b", 0.2);
      fTPCClusterResid__CROC_m = p.get<float>("TPCClusterResid__CROC_m", 0.1);
      fTPCClusterResid__IROC_b = p.get<float>("TPCClusterResid__IROC_b", 0.1);
      fTPCClusterResid__IROC_m = p.get<float>("TPCClusterResid__IROC_m", 0.2);
      fTPCClusterResid_IOROC_b = p.get<float>("TPCClusterResid__CROC_b", 0.1);
      fTPCClusterResid_IOROC_m = p.get<float>("TPCClusterResid__CROC_m", 0.3);
      fTPCClusterResid_OOROC_b = p.get<float>("TPCClusterResid__CROC_b", 0.1);
      fTPCClusterResid_OOROC_m = p.get<float>("TPCClusterResid__CROC_m", 0.9);

      art::InputTag patrecTag(fPatRecLabel);
      consumes<std::vector<gar::rec::Track>>(patrecTag);
      consumes<art::Assns<gar::rec::TPCCluster, gar::rec::Track>>(patrecTag);

      // probably don't need the vector hits at this point if we have the TPCClusters
      //consumes< std::vector<gar::rec::VecHit> >(patrecTag);
      //consumes< art::Assns<gar::rec::VecHit, gar::rec::Track> >(patrecTag);

      produces<std::vector<gar::rec::Track>>();
      produces<art::Assns<gar::rec::TPCCluster, gar::rec::Track>>();
      produces<std::vector<gar::rec::TrackIoniz>>();
      produces<art::Assns<rec::TrackIoniz, rec::Track>>();
      produces<std::vector<gar::rec::TrackTrajectory>>();
      produces<art::Assns<rec::TrackTrajectory, rec::Track>>();
    }



    void tpctrackfit2::produce(art::Event& e)
    {
      // output collections

      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);
      std::unique_ptr< std::vector<rec::TrackIoniz> > ionCol(new std::vector<rec::TrackIoniz>);
      std::unique_ptr< art::Assns<rec::TrackIoniz,rec::Track> > ionTrkAssns(new ::art::Assns<rec::TrackIoniz,rec::Track>);
      std::unique_ptr< std::vector<rec::TrackTrajectory> > trajCol(new std::vector<rec::TrackTrajectory>);
      std::unique_ptr< art::Assns<rec::TrackTrajectory,rec::Track> > trajTrkAssns(new ::art::Assns<rec::TrackTrajectory,rec::Track>);

      // inputs

      auto patrecTrackHandle = e.getValidHandle< std::vector<gar::rec::Track> >(fPatRecLabel);
      auto const& patrecTracks = *patrecTrackHandle;

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      auto const ionizPtrMaker = art::PtrMaker<rec::TrackIoniz>(e);
      auto const trajPtrMaker  = art::PtrMaker<rec::TrackTrajectory>(e);
      //auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());

      auto const *magFieldService = gar::providerFrom<mag::MagneticFieldService>();
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      double xtpccent = euclid->TPCXCent();
      double ytpccent = euclid->TPCYCent();
      double ztpccent = euclid->TPCZCent();
      TVector3 tpccent(xtpccent,ytpccent,ztpccent);
      TVector3 xhat(1,0,0);

      const art::FindManyP<gar::rec::TPCCluster> TPCClustersFromPatRecTracks(patrecTrackHandle,e,fPatRecLabel);

      for (size_t itrack = 0; itrack < patrecTracks.size(); ++itrack)
        {
          std::vector<gar::rec::TPCCluster> TPCClusters;
          for (size_t iTPCCluster=0; iTPCCluster < TPCClustersFromPatRecTracks.at(itrack).size(); ++iTPCCluster)
            {
              TPCClusters.push_back(*TPCClustersFromPatRecTracks.at(itrack).at(iTPCCluster));  // make our own local copy of TPCClusters.  Maybe we can skip this?
            }
          std::vector<TrackPar> trackparams;
          std::vector<TrackIoniz> trackions;
          std::vector<TrackTrajectory> tracktrajs;
          if (KalmanFitBothWays(TPCClusters,trackparams,trackions,tracktrajs) == 0)   // to think about -- unused TPCClusters?  Or just ignore them in the fit?
            {
              gar::rec::IDNumber ID = 0;
              for(size_t ireco = 0; ireco < trackparams.size(); ++ireco)
                {
                 Track tr = trackparams[ireco].CreateTrack();
                 if(ireco==0) ID = tr.getIDNumber();
                 tr.setIDNumber(ID);
                 //std::cout<<"Track ID: "<<tr.getIDNumber()<<std::endl;
                 trkCol->push_back(tr);
                 ionCol->push_back(trackions[ireco]);
                 trajCol->push_back(tracktrajs[ireco]);
                 auto const trackpointer = trackPtrMaker(trkCol->size()-1);
                 auto const ionizpointer = ionizPtrMaker(ionCol->size()-1);
                 auto const trajpointer  = trajPtrMaker(trajCol->size()-1);
                 for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
                    {
                     TPCClusterTrkAssns->addSingle(TPCClustersFromPatRecTracks.at(itrack).at(iTPCCluster),trackpointer);
                    }
              	 ionTrkAssns->addSingle(ionizpointer, trackpointer);
              	 trajTrkAssns->addSingle(trajpointer, trackpointer);
                }
            }
        }

      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterTrkAssns));
      e.put(std::move(ionCol));
      e.put(std::move(ionTrkAssns));
      e.put(std::move(trajCol));
      e.put(std::move(trajTrkAssns));
    }

    double tpctrackfit2::CalculateChi2(RVec<AliExternalTrackParam4D> ParamMC, RVec<AliExternalTrackParam4D> ParamIn, double xstart, double ystart, double * GArCenter)
    {
      double chisquared = 0;
      for(size_t t = 0; t<ParamIn.size(); ++t)
         {
          double xyz[3];
          ParamIn[t].GetXYZ(xyz);
          double xyzMC[3];
          ParamMC[t].GetXYZ(xyzMC);

          if(xyz[0]==0 || xyz[1]==0 || xyz[2]==0) continue;
          if(xyzMC[0]==0 || xyzMC[1]==0 || xyzMC[2]==0) continue;
          
          TVectorF ytilde(2);
               
          ytilde[0] = xyzMC[1] - xyz[1];
          ytilde[1] = xyzMC[0] - xyz[0];


          float impactAngle;
          TVector3 trajPerp(0.0, xyz[1]+(GArCenter[1]-ystart),xyz[0]+(GArCenter[2]-xstart));
          float rTrj = trajPerp.Mag();

          Double_t cb=TMath::Cos(-ParamIn[t].GetAlpha());
          Double_t sb=TMath::Sin(-ParamIn[t].GetAlpha());
          Double_t sfb=ParamIn[t].GetParameter()[2];
          Double_t cfb=TMath::Sqrt((1.- sfb)*(1.+sfb));
          Double_t sfrotb = sfb*cb - cfb*sb;
          
          TVector3 trajStepPerp(0.0, sfrotb, TMath::Sqrt(1-sfrotb*sfrotb));
          impactAngle = trajPerp.Dot(trajStepPerp) / rTrj;
          impactAngle = acos(abs(impactAngle));
          float IROC_OROC_boundary = (euclid->GetIROCOuterRadius() +euclid->GetOROCInnerRadius())/2.0;
          bool In_CROC =                                                  rTrj <= euclid->GetIROCInnerRadius();
          bool In_IROC = euclid->GetIROCInnerRadius() < rTrj 		   && rTrj <= IROC_OROC_boundary;
          bool InIOROC = IROC_OROC_boundary < rTrj 		               && rTrj <= euclid->GetOROCPadHeightChangeRadius();
          bool InOOROC = euclid->GetOROCPadHeightChangeRadius() < rTrj;
          float typicalResidual = 1.0;	// Shaddup, compiler
          if (In_CROC) {
            typicalResidual = fTPCClusterResid__CROC_m*impactAngle +fTPCClusterResid__CROC_b;
          } else if (In_IROC) {
            typicalResidual = fTPCClusterResid__IROC_m*impactAngle +fTPCClusterResid__IROC_b;
          } else if (InIOROC) {
            typicalResidual = fTPCClusterResid_IOROC_m*impactAngle +fTPCClusterResid_IOROC_b;
          } else if (InOOROC) {
            typicalResidual = fTPCClusterResid_OOROC_m*impactAngle +fTPCClusterResid_OOROC_b;
          }

          chisquared += ytilde.Norm2Sqr()/TMath::Sq(typicalResidual);
         }
      
      return chisquared;
    }

    void tpctrackfit2::CalculateTrackInfo(float &length, std::vector<std::pair<float,float>>& dSigdXs, std::vector<TVector3>& trajpts, 
                                          std::vector<int> TPCClusterList, std::vector<TPCCluster> TPCClusters,
                                          RVec<AliExternalTrackParam4D> ParamIn, double xstart, double ystart, double * GArCenter)
    {     
      for(size_t t = 0; t < ParamIn.size(); t++)
         { 
          double xyz[3];
          ParamIn[t].GetXYZ(xyz);
                
          if(xyz[0]==0 || xyz[1]==0 || xyz[2]==0) continue;
          trajpts.emplace_back(xyz[2], xyz[1]+(GArCenter[1]-ystart),xyz[0]+(GArCenter[2]-xstart));

          if(t!=0)
            {
              double xyz_prev[3];
              ParamIn[t-1].GetXYZ(xyz_prev);

              if(xyz[0]==0 || xyz[1]==0 || xyz[2]==0) continue;
              float d_length = TMath::Sqrt( TMath::Sq(xyz[0]-xyz_prev[0]) + TMath::Sq(xyz[1]-xyz_prev[1]) + TMath::Sq(xyz[2]-xyz_prev[2]) );
              length += d_length;
              
              float valSig = TPCClusters[TPCClusterList[t]].Signal();
              //std::cout<<"TPCClusterPos: "<< TPCClusters[TPCClusterList[t]].Position()[0]<<"\n";
              //std::cout<<"Param Position: "<< xyz[2]<<"\n";
              if (d_length < fMinIonizGapCut)
               {
                 std::pair pushme = std::make_pair(valSig,d_length);
                 dSigdXs.push_back( pushme );
               }
              else
               {
                 if (dSigdXs.size()>0) dSigdXs.pop_back();
               }
            }
         }
          
    }
    int tpctrackfit2::KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
                                        std::vector<TrackPar> &trackpars, std::vector<TrackIoniz> &trackions, std::vector<TrackTrajectory> &tracktrajs)

    {
      // For Garsoft Implementation:
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end
      //
      // For actual ALICE KF reconstruction
      // variables: z is independent variable with rotating frame in radial coordinates
      // 0: y
      // 1: x
      // 2: sinphi
      // 3: tanlambda
      // 4: q/pT


      std::vector<int> hlf;
      std::vector<int> hlb;
      float lftmp=0;  // need these dummy arguments so we can share code with the patrec sorter
      float lbtmp=0;
      if (fSortAlg == 1)
	{
          gar::rec::sort_TPCClusters_along_track(TPCClusters,hlf,hlb,fPrintLevel,lftmp,lbtmp,fSortTransWeight,fSortDistBack);
	}
      else if (fSortAlg == 2)
	{
	  gar::rec::sort_TPCClusters_along_track2(TPCClusters,hlf,hlb,fPrintLevel,lftmp,lbtmp,fSortDistCut);
	}
      else
	{
	  throw cet::exception("tpctrackfit2_module") << "Sort Algorithm switch not understood: " << fSortAlg;
	}
      if (hlf.size() == 0) return 1;

      /////////////////////////////////////////////////////////////////////////////////////////////New Code
     
      //////////////Hacky implementation of global params, to be done more properly with fcl files
      float resol[2]={0.0001,0.0001};
      resol[0]=0.4;                   
      resol[1]=0.3;
      const Int_t   nLayerTPC=278;
      const Float_t xx0=8.37758e-04; //1/X0 cm^-1 for ArCH4 at 10 atm
      const Float_t xrho=0.016770000; //rho g/cm^3 for ArCH4 at 10 atm
      double GArCenter[3]={0,-150.473,1486};
      int PDGcode = 11;

      //////////////Building fast geometry with TPC properties
      fastGeometry geom(nLayerTPC+1);
      geom.fBz=-5;
      geom.setLayerRadiusPower(0,nLayerTPC,1,nLayerTPC,1.0,xx0,xrho,resol);
      for (size_t iLayer=0; iLayer<geom.fLayerX0.size();iLayer++) 
        {
          geom.fLayerX0[iLayer] = xx0;
          geom.fLayerRho[iLayer] = xrho;
          geom.fLayerResolRPhi[iLayer] = resol[0];
          geom.fLayerResolZ[iLayer] = resol[1];
        }
 
      ////////////Creating Position vectors to be used in reconstruction
      std::vector<TVector3> TrkClusterXYZb_NDGAr;
      std::vector<TVector3> TrkClusterXYZf_NDGAr;
      for(size_t k=0;k<hlb.size();k++) 
        {
          TVector3 pos(TPCClusters[hlb[k]].Position()[0],TPCClusters[hlb[k]].Position()[1],TPCClusters[hlb[k]].Position()[2]);
          TrkClusterXYZb_NDGAr.push_back(pos);
        }
      for(size_t k=0;k<hlf.size();k++) 
        {
           TVector3 pos(TPCClusters[hlf[k]].Position()[0],TPCClusters[hlf[k]].Position()[1],TPCClusters[hlf[k]].Position()[2]);
           TrkClusterXYZf_NDGAr.push_back(pos);
        }
      
      //////////////////Forward trajectory reconstruction
      fastParticle particle_f(hlf.size()+1);
      particle_f.fAddMSsmearing=true;
      particle_f.fAddPadsmearing=false;
      particle_f.fUseMCInfo=false;

      ///////////////Choosing frame of reference
      double xstart=0;
      double ystart=0;          
      size_t size_dis = 30;
      if (int(TrkClusterXYZf_NDGAr.size()-1)<30) size_dis = int(TrkClusterXYZf_NDGAr.size()-1);
      double displacex = (TrkClusterXYZf_NDGAr.at(size_dis).Z()-TrkClusterXYZf_NDGAr.at(0).Z());
      double displacey = (TrkClusterXYZf_NDGAr.at(size_dis).Y()-TrkClusterXYZf_NDGAr.at(0).Y());
      double displace_mod = sqrt(displacex*displacex+displacey*displacey);
      xstart = -(TrkClusterXYZf_NDGAr.at(0).Z()-GArCenter[2])+20*displacex/displace_mod;
      ystart = -(TrkClusterXYZf_NDGAr.at(0).Y()-GArCenter[1])+20*displacey/displace_mod; 
      BuildParticlePoints(particle_f,GArCenter,geom,TrkClusterXYZf_NDGAr,PDGcode,xstart,ystart);
      
      //particle_f.reconstructParticleFullOut(geom,PDGcode,10000); ///outwards reconstruction
      //particle_f.reconstructParticleFull(geom,PDGcode,10000);
      //std::cout<<"NTPCClusters "<<particle_f.fChi2.size()<<std::endl;
      //std::cout<<"Chi2 Forward In/Out "<<particle_f.fChi2[0]<<" "<<particle_f.fChi2Out[particle_f.fChi2Out.size()-1]<<std::endl;


      //////////////////Backward Trajectory reconstruction
      fastParticle particle_b(hlb.size()+1);
      particle_b.fAddMSsmearing=true;
      particle_b.fAddPadsmearing=false;
      particle_b.fUseMCInfo=false;


      //////////////Choosing frame of reference
      double xend=0;
      double yend=0;
      size_t size_dis_b = 30;
      if (int(TrkClusterXYZb_NDGAr.size()-1)<30) size_dis_b = int(TrkClusterXYZb_NDGAr.size()-1);
      double displacex_b = (TrkClusterXYZb_NDGAr.at(size_dis_b).Z()-TrkClusterXYZb_NDGAr.at(0).Z());
      double displacey_b = (TrkClusterXYZb_NDGAr.at(size_dis_b).Y()-TrkClusterXYZb_NDGAr.at(0).Y());
      double displace_mod_b = sqrt(displacex_b*displacex_b+displacey_b*displacey_b);
      xend = -(TrkClusterXYZb_NDGAr.at(0).Z()-GArCenter[2])+20*displacex_b/displace_mod_b;
      yend = -(TrkClusterXYZb_NDGAr.at(0).Y()-GArCenter[1])+20*displacey_b/displace_mod_b;
      BuildParticlePoints(particle_b,GArCenter,geom,TrkClusterXYZb_NDGAr,PDGcode,xend,yend);
    
      //particle_b.reconstructParticleFull(geom,PDGcode,10000);
      //particle_b.reconstructParticleFullOut(geom,PDGcode,10000);
      //std::cout<<"Chi2 Backwards In/Out "<<particle_b.fChi2[0]<<" "<<particle_b.fChi2Out[particle_b.fChi2Out.size()-1]<<std::endl<<std::endl;


      ///Test Combinations
      int PDGcodes[]={11,13,211,2212};
      //double bestcombinedchi2=100000;
      //int bestPDG=0;
      //bool forward=1;
      std::vector<fastParticle> particle_h;

      for (auto p: PDGcodes){
        fastParticle pf = particle_f;
        fastParticle pb = particle_b;

        pf.reconstructParticleFull(geom,p,10000);
        pf.reconstructParticleFullOut(geom,p,10000);
        pb.reconstructParticleFull(geom,p,10000);
        pb.reconstructParticleFullOut(geom,p,10000);
 
        particle_h.push_back(pf);
        particle_h.push_back(pb);
       
     //   double chi2fIn = (particle_f.fChi2[0]<0.000001)? 100000 : particle_f.fChi2[0];
     //   double chi2fOut = (particle_f.fChi2Out[particle_f.fChi2Out.size()-1]<0.000001)? 100000 : particle_f.fChi2Out[particle_f.fChi2Out.size()-1];
     //   double combinedchi2f = chi2fIn+chi2fOut;

     //   double chi2bIn = (particle_b.fChi2[0]<0.000001)? 100000 : particle_b.fChi2[0];
     //   double chi2bOut = (particle_b.fChi2Out[particle_b.fChi2Out.size()-1]<0.000001)? 100000 : particle_b.fChi2Out[particle_b.fChi2Out.size()-1];
     //   double combinedchi2b = chi2bIn+chi2bOut;
     //      
     //   double bestchi2fb = 0;
     //   bool f = 0;
     //   if(combinedchi2f<combinedchi2b){
     //     bestchi2fb = combinedchi2f;
     //     f = 1; 
     //   }else{
     //     bestchi2fb = combinedchi2b;
     //     f = 0;
     //   }
     //   std::cout<<"combined chi2 : PDG "<<bestchi2fb<<" "<<p<<std::endl;
     //   if(bestchi2fb<bestcombinedchi2){
     //     forward=f;
     //     bestPDG=p;
     //     bestcombinedchi2 = bestchi2fb;
     //   }
       }

     // std::cout<<"Best Results (chi2, fb, PDG): "<<bestcombinedchi2<<" "<<forward<<" "<<bestPDG<<std::endl<<std::endl;
     // fastParticle particle;

     // if(forward){
     //   particle = particle_f;
     //   particle.reconstructParticleFull(geom,bestPDG,10000);
     //   particle.reconstructParticleFullOut(geom,bestPDG,10000);
     // }else{
     //   particle = particle_b;
     //   particle.reconstructParticleFull(geom,bestPDG,10000);
     //   particle.reconstructParticleFullOut(geom,bestPDG,10000);
     //   xstart=xend;
     //   ystart=yend;
     // }
     // if(bestPDG==0) bestPDG=11;
      ////////////////////////////////////////////////////////////////////////////////////////////////////

      for(size_t h=0; h<particle_h.size(); ++h)
        { 

          double xinit, yinit = 0;
          std::vector<int> TPClist;
          if(h % 2 == 0)
            {
             xinit = xstart;
             yinit = ystart;
             TPClist = hlf;
            }
          else
            {
             xinit = xend;
             yinit = yend;
             TPClist = hlb;
            }

          std::vector<float> tparbeg(6,0);
          float covmatbeg[25];
          float chisqbackwards = 0;
          float lengthbackwards = 0;
          std::vector<std::pair<float,float>> dSigdXs_BAK;
          std::vector<TVector3> trajpts_BAK;

          Double_t xyz_start[3];
          particle_h[h].fParamIn[0].GetXYZ(xyz_start);
          Double_t cb=TMath::Cos(-particle_h[h].fParamIn[0].GetAlpha());
          Double_t sb=TMath::Sin(-particle_h[h].fParamIn[0].GetAlpha());
          Double_t sfb=particle_h[h].fParamIn[0].GetParameter()[2];
          Double_t cfb=TMath::Sqrt((1.- sfb)*(1.+sfb));
          Double_t sfrotb = sfb*cb - cfb*sb;

          tparbeg[0]=xyz_start[1]+(GArCenter[1]-yinit);
          tparbeg[1]=xyz_start[0]+(GArCenter[2]-xinit);
          tparbeg[2]=particle_h[h].fParamIn[0].GetParameter()[4]*(5*0.299792458e-3);
          tparbeg[3]=TMath::ASin(sfrotb);
          tparbeg[4]=TMath::ATan(particle_h[h].fParamIn[particle_h[h].fParamIn.size()-1].GetParameter()[3]);
          tparbeg[5]=xyz_start[2];

          chisqbackwards=CalculateChi2(particle_h[h].fParamMC,particle_h[h].fParamIn,xinit,yinit,GArCenter);
          CalculateTrackInfo(lengthbackwards,dSigdXs_BAK,trajpts_BAK,TPClist,TPCClusters,particle_h[h].fParamIn,xinit,yinit,GArCenter);
      
          std::vector<float> tparend(6,0);
          float covmatend[25];
          float chisqforwards = 0;
          float lengthforwards = 0;
          std::set<int> unused_TPCClusters;
          std::vector<std::pair<float,float>> dSigdXs_FWD;
          std::vector<TVector3> trajpts_FWD;


          Double_t xyz_end[3];
          particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetXYZ(xyz_end);
          Double_t ca=TMath::Cos(-particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetAlpha());
          Double_t sa=TMath::Sin(-particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetAlpha());
          Double_t sf=particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetParameter()[2];
          Double_t cf=TMath::Sqrt((1.- sf)*(1.+sf));
          Double_t sfrot = sf*ca - cf*sa;
          tparend[0]=xyz_end[1]+(GArCenter[1]-yinit);
          tparend[1]=xyz_end[0]+(GArCenter[2]-xinit);
          tparend[2]=particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetParameter()[4]*(5*0.299792458e-3);
          tparend[3]=TMath::ASin(sfrot);
          tparend[4]=TMath::ATan(particle_h[h].fParamOut[particle_h[h].fParamOut.size()-1].GetParameter()[3]);
          tparend[5]=xyz_end[2];

          chisqforwards=CalculateChi2(particle_h[h].fParamMC,particle_h[h].fParamOut,xinit,yinit,GArCenter);            
          CalculateTrackInfo(lengthforwards,dSigdXs_FWD,trajpts_FWD,TPClist,TPCClusters,particle_h[h].fParamOut,xinit,yinit,GArCenter);

          size_t nTPCClusters=0;
          if (TPCClusters.size()>unused_TPCClusters.size())
           { nTPCClusters = TPCClusters.size()-unused_TPCClusters.size(); }
          TrackPar trackpar;
          trackpar.setNTPCClusters(nTPCClusters);
          trackpar.setTime(0);
          trackpar.setChisqForwards(chisqforwards);
          trackpar.setChisqBackwards(chisqbackwards);
          trackpar.setLengthForwards(lengthforwards);
          trackpar.setLengthBackwards(lengthbackwards);
          trackpar.setCovMatBeg(covmatbeg);
          trackpar.setCovMatEnd(covmatend);
          trackpar.setTrackParametersBegin(tparbeg.data());
          trackpar.setXBeg(tparbeg[5]);
          trackpar.setTrackParametersEnd(tparend.data());
          trackpar.setXEnd(tparend[5]);

          TrackIoniz trackion;
          TrackTrajectory tracktraj;     
          trackion.setData(dSigdXs_FWD,dSigdXs_BAK);
          tracktraj.setData(trajpts_FWD,trajpts_BAK);
         
          trackpars.push_back(trackpar);
          trackions.push_back(trackion);
          tracktrajs.push_back(tracktraj);
         }
      return 0;
    }


    //--------------------------------------------------------------------------------------------------------------


    DEFINE_ART_MODULE(tpctrackfit2)

  } // namespace rec
} // namespace gar
