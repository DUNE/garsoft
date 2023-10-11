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



      int KalmanFit( std::vector<TPCCluster> &TPCClusters,
                     std::vector<int> &TPCClusterlist,
                     std::vector<float> &trackparatend,
                     float &chisquared,
                     float &length,
                     float *covmat,    // 5x5 covariance matrix
                     std::set<int> &unused_TPCClusters,
                     std::vector<std::pair<float,float>>& dSigdX,
                     std::vector<TVector3>& trajpts);

      int KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
                            TrackPar &trackpar,  TrackIoniz &trackions, TrackTrajectory &tracktraj);

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
          TrackPar trackparams;
          TrackIoniz trackions;
          TrackTrajectory tracktraj;
          if (KalmanFitBothWays(TPCClusters,trackparams,trackions,tracktraj) == 0)   // to think about -- unused TPCClusters?  Or just ignore them in the fit?
            {
              trkCol->push_back(trackparams.CreateTrack());
              ionCol->push_back(trackions);
              trajCol->push_back(tracktraj);
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

      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterTrkAssns));
      e.put(std::move(ionCol));
      e.put(std::move(ionTrkAssns));
      e.put(std::move(trajCol));
      e.put(std::move(trajTrkAssns));
    }

    int tpctrackfit2::KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
                                        TrackPar &trackpar, TrackIoniz &trackions, TrackTrajectory &tracktraj)

    {
      // variables:  x is the independent variable
      // 0: y
      // 1: z
      // 2: curvature
      // 3: phi
      // 4: lambda = angle from the cathode plane
      // 5: x   /// added on to the end


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
     
      float resol[2]={0.0001,0.0001};
      resol[0]=0.4;                   
      resol[1]=0.3;
      const Int_t   nLayerTPC=278;
      const Float_t xx0=8.37758e-04; //1/X0 cm^-1 for ArCH4 at 10 atm
      const Float_t xrho=0.016770000; //rho g/cm^3 for ArCH4 at 10 atm
      double GArCenter[3]={0,-150.473,1486};
      int PDGcode = 211;

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
      
      
      fastParticle particle_f(hlf.size()+1);
      particle_f.fAddMSsmearing=true;
      particle_f.fAddPadsmearing=false;
      particle_f.fUseMCInfo=false;

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
      std::cout<<"NTPCClusters "<<TrkClusterXYZf_NDGAr.size()<<std::endl;
      size_t NTPCClusters = TrkClusterXYZf_NDGAr.size();
      particle_f.reconstructParticleFull(geom,PDGcode,10000);
      particle_f.reconstructParticleFullOut(geom,PDGcode,10000);

      Double_t xyz_beg[3];
      particle_f.fParamIn[0].GetXYZ(xyz_beg);
      Double_t xyz_end_out[3];
      particle_f.fParamOut[particle_f.fParamOut.size()-1].GetXYZ(xyz_end_out);
      //std::cout<<"Starting par ALICE (x,y,z,sinPhi,tanlambda,q/pT): "<<xyz_beg[2]<<" "<<xyz_beg[1]+(GArCenter[1]-ystart)<<" "<<xyz_beg[0]+(GArCenter[2]-xstart)<<" ";
      //std::cout<<particle_f.fParamIn[0].GetParameter()[2]<<" "<<particle_f.fParamIn[0].GetParameter()[3]<<" "<<particle_f.fParamIn[0].GetParameter()[4]<<"\n"; 
      Double_t pxyz_beg[3];
      particle_f.fParamIn[0].GetPxPyPz(pxyz_beg);
      //float sign = +1;
      //if ( pxyz_beg[2]>0 && particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]<particle_f.fParamMC[0].GetParameter()[1]) sign = -1;
      //if ( pxyz_beg[2]<0 && particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]>particle_f.fParamMC[0].GetParameter()[1]) sign = -1;
      //std::cout<<"xstart ystart "<<xstart<<" "<<ystart<<std::endl;
      //std::cout<<"ALICE end/start x: "<<particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]<<" "<<particle_f.fParamMC[0].GetParameter()[1]<<std::endl;
      //Double_t mod_p = sign*particle_f.fParamIn[0].GetP();
      Double_t ca=TMath::Cos(-particle_f.fParamIn[0].GetAlpha()), sa=TMath::Sin(-particle_f.fParamIn[0].GetAlpha());
      Double_t sf=particle_f.fParamIn[0].GetParameter()[2];
      Double_t cf=TMath::Sqrt((1.- sf)*(1.+sf));
      Double_t sfrot = sf*ca - cf*sa;
      //std::cout<<"Starting dir ALICE: "<<pxyz_beg[2]/mod_p<<" "<<pxyz_beg[1]/mod_p<<" "<<pxyz_beg[0]/mod_p<<"\n";
      //std::cout<<"Starting par ALICE conv (x,y,z,q/r,sinphi,tanlambda): "<<xyz_beg[2]<<" "<<xyz_beg[1]+(GArCenter[1]-ystart)<<" "<<xyz_beg[0]+(GArCenter[2]-xstart)<<" ";
      //std::cout<<" "<<particle_f.fParamIn[0].GetParameter()[4]*(5*0.299792458e-3)<<" "<<sfrot<<" "<<particle_f.fParamIn[0].GetParameter()[3]<<"\n";
      
      ca=TMath::Cos(-particle_f.fParamOut[particle_f.fParamOut.size()-1].GetAlpha()), sa=TMath::Sin(-particle_f.fParamOut[particle_f.fParamOut.size()-1].GetAlpha());
      sf=particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[2];
      cf=TMath::Sqrt((1.- sf)*(1.+sf));
      sfrot = sf*ca - cf*sa;
      //std::cout<<"Starting dir ALICE: "<<pxyz_beg[2]/mod_p<<" "<<pxyz_beg[1]/mod_p<<" "<<pxyz_beg[0]/mod_p<<"\n";
      if(NTPCClusters>100) std::cout<<"Ending par ALICE Out conv (x,y,z,q/r,sinphi,tanlambda): "<<xyz_end_out[2]<<" "<<xyz_end_out[1]+(GArCenter[1]-ystart)<<" "<<xyz_end_out[0]+(GArCenter[2]-xstart)<<" ";
      if(NTPCClusters>100) std::cout<<" "<<particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[4]*(5*0.299792458e-3)<<" "<<sfrot<<" "<<particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[3]<<"\n";

      double dira[3];
      dira[0] = TMath::Tan(-TMath::ATan(particle_f.fParamIn[0].GetParameter()[3]));
      dira[1] = TMath::Sin(TMath::ASin(sfrot));
      dira[2] = TMath::Cos(TMath::ASin(sfrot));
      float signa = +1;
      if ( dira[0]>0 && particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]<particle_f.fParamMC[0].GetParameter()[1]) signa = -1;
      if ( dira[0]<0 && particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]>particle_f.fParamMC[0].GetParameter()[1]) signa = -1;
      float norma = signa * TMath::Sqrt( 1.0 + dira[0]*dira[0]);
      dira[0] /= norma;
      dira[1] /= norma;
      dira[2] /= norma;
      //std::cout<<"Starting dir ALICE conv (xyz): "<<dira[0]<<" "<<dira[1]<<" "<<dira[2]<<std::endl;
      
      
      fastParticle particle_b(hlb.size()+1);
      particle_b.fAddMSsmearing=true;
      particle_b.fAddPadsmearing=false;
      particle_b.fUseMCInfo=false;

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
      //std::cout<<"NTPCClusters "<<TrkClusterXYZb_NDGAr.size()<<std::endl;
      particle_b.reconstructParticleFull(geom,PDGcode,10000);
      particle_b.reconstructParticleFullOut(geom,PDGcode,10000);


      Double_t xyz_end[3];
      particle_b.fParamIn[0].GetXYZ(xyz_end);
      Double_t xyz_start_out[3];
      particle_b.fParamOut[particle_b.fParamOut.size()-1].GetXYZ(xyz_start_out);
      //std::cout<<"Ending par ALICE (x,y,z,sinPhi,tanlambda,q/pT): "<<xyz_end[2]<<" "<<xyz_end[1]+(GArCenter[1]-yend)<<" "<<xyz_end[0]+(GArCenter[2]-xend)<<" ";
      //std::cout<<particle_b.fParamIn[0].GetParameter()[2]<<" "<<particle_b.fParamIn[0].GetParameter()[3]<<" "<<particle_b.fParamIn[0].GetParameter()[4]<<"\n";
      


      Double_t pxyz_end[3];
      particle_b.fParamIn[0].GetPxPyPz(pxyz_end);
      //float signb = +1;
      //if ( pxyz_end[2]>0 && particle_b.fParamMC[particle_b.fParamMC.size()-1].GetParameter()[1]<particle_b.fParamMC[0].GetParameter()[1]) signb = -1;
      //if ( pxyz_end[2]<0 && particle_b.fParamMC[particle_b.fParamMC.size()-1].GetParameter()[1]>particle_b.fParamMC[0].GetParameter()[1]) signb = -1;
      //std::cout<<"xstart ystart "<<xstart<<" "<<ystart<<std::endl;
      //std::cout<<"ALICE end/start x: "<<particle_f.fParamMC[particle_f.fParamMC.size()-1].GetParameter()[1]<<" "<<particle_f.fParamMC[0].GetParameter()[1]<<std::endl;
      //Double_t mod_pend = signb*particle_b.fParamIn[0].GetP();
      Double_t cb=TMath::Cos(-particle_b.fParamIn[0].GetAlpha()), sb=TMath::Sin(-particle_b.fParamIn[0].GetAlpha());
      Double_t sfb=particle_b.fParamIn[0].GetParameter()[2];
      Double_t cfb=TMath::Sqrt((1.- sfb)*(1.+sfb));
      Double_t sfrotb = sfb*cb - cfb*sb;
      //std::cout<<"Ending dir ALICE: "<<pxyz_end[2]/mod_pend<<" "<<pxyz_end[1]/mod_pend<<" "<<pxyz_end[0]/mod_pend<<"\n";
      //std::cout<<"Ending par ALICE conv (x,y,z,q/r,inhi,lambda): "<<xyz_end[2]<<" "<<xyz_end[1]+(GArCenter[1]-yend)<<" "<<xyz_end[0]+(GArCenter[2]-xend)<<" ";
      //std::cout<<" "<<particle_b.fParamIn[0].GetParameter()[4]*(5*0.299792458e-3)<<" "<<sfrotb<<" "<<TMath::ATan(particle_b.fParamIn[0].GetParameter()[3])<<"\n";
      
      cb=TMath::Cos(-particle_b.fParamOut[particle_b.fParamOut.size()-1].GetAlpha()), sa=TMath::Sin(-particle_b.fParamOut[particle_b.fParamOut.size()-1].GetAlpha());
      sfb=particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[2];
      cfb=TMath::Sqrt((1.- sf)*(1.+sf));
      sfrotb = sfb*cb - cfb*sb;
      //std::cout<<"Starting dir ALICE: "<<pxyz_beg[2]/mod_p<<" "<<pxyz_beg[1]/mod_p<<" "<<pxyz_beg[0]/mod_p<<"\n";
      if(NTPCClusters>100) std::cout<<"Starting par ALICE Out conv (x,y,z,q/r,sinphi,tanlambda): "<<xyz_start_out[2]<<" "<<xyz_start_out[1]+(GArCenter[1]-yend)<<" "<<xyz_start_out[0]+(GArCenter[2]-xend)<<" ";
      if(NTPCClusters>100) std::cout<<" "<<particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[4]*(5*0.299792458e-3)<<" "<<sfrotb<<" "<<particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[3]<<"\n";

      double dirb[3];
      dirb[0] = TMath::Tan(-TMath::ATan(particle_b.fParamIn[0].GetParameter()[3]));
      dirb[1] = TMath::Sin(TMath::ASin(sfrotb));
      dirb[2] = TMath::Cos(TMath::ASin(sfrotb));
      float signba = +1;
      if ( dirb[0]>0 && particle_b.fParamMC[particle_b.fParamMC.size()-1].GetParameter()[1]<particle_b.fParamMC[0].GetParameter()[1]) signba = -1;
      if ( dirb[0]<0 && particle_b.fParamMC[particle_b.fParamMC.size()-1].GetParameter()[1]>particle_b.fParamMC[0].GetParameter()[1]) signba = -1;
      float normb = signba * TMath::Sqrt( 1.0 + dirb[0]*dirb[0]);
      dirb[0] /= normb;
      dirb[1] /= normb;
      dirb[2] /= normb;
      //std::cout<<"Ending dir ALICE conv (xyz): "<<dirb[0]<<" "<<dirb[1]<<" "<<dirb[2]<<std::endl;



      ////////////////////////////////////////////////////////////////////////////////////////////////////
      
      std::vector<float> tparend(6);
      float covmatend[25];
      float chisqforwards = 0;
      float lengthforwards = 0;
      std::set<int> unused_TPCClusters;
      std::vector<std::pair<float,float>> dSigdXs_FWD;
      std::vector<TVector3> trajpts_FWD;

      int retcode = KalmanFit(TPCClusters,hlf,tparend,chisqforwards,lengthforwards,covmatend,unused_TPCClusters,dSigdXs_FWD,trajpts_FWD);
      if (retcode != 0) return 1;

      // the "backwards" fit is in decreasing x.  Track parameters are at the end of the fit, the other end of the track

      std::vector<float> tparbeg(6);
      float covmatbeg[25];
      float chisqbackwards = 0;
      float lengthbackwards = 0;
      std::vector<std::pair<float,float>> dSigdXs_BAK;
      std::vector<TVector3> trajpts_BAK;

      retcode = KalmanFit(TPCClusters,hlb,tparbeg,chisqbackwards,lengthbackwards,covmatbeg,unused_TPCClusters,dSigdXs_BAK,trajpts_BAK);
      if (retcode != 0) return 1;
      double dir[3];
      dir[0] = TMath::Tan(tparbeg[4]);
      dir[1] = TMath::Sin(tparbeg[3]);
      dir[2] = TMath::Cos(tparbeg[3]);
      float sigh = +1;
      if ( dir[0]>0 && tparend[5]<tparbeg[5]) sigh = -1;
      if ( dir[0]<0 && tparend[5]>tparbeg[5] ) sigh = -1;
      float norm = sigh * TMath::Sqrt( 1.0 + dir[0]*dir[0]);
      dir[0] /= norm;
      dir[1] /= norm;
      dir[2] /= norm;
      //std::cout<<"GAr end/start x: "<<tparend[5]<<" "<<tparbeg[5]<<std::endl;
      //std::cout<<"GAr end/start y: "<<tparend[0]<<" "<<tparbeg[0]<<std::endl;
      if(NTPCClusters>100) std::cout<<"Starting par GAr (x,y,z,q/r,sinphi,tanlambda): "<<tparbeg[5]<<" "<<tparbeg[0]<<" "<<tparbeg[1]<<" "<<tparbeg[2]<<" "<<TMath::Sin(tparbeg[3])<<" "<<TMath::Sin(tparbeg[4])<<"\n";
     

      double dirgb[3];
      dirgb[0] = TMath::Tan(tparend[4]);
      dirgb[1] = TMath::Sin(tparend[3]);
      dirgb[2] = TMath::Cos(tparend[3]);
      float sighb = +1;
      if ( dirgb[0]>0 && tparend[5]<tparend[5]) sighb = -1;
      if ( dirgb[0]<0 && tparend[5]>tparend[5] ) sighb = -1;
      float normgb = sighb * TMath::Sqrt( 1.0 + dirgb[0]*dirgb[0]);
      dirgb[0] /= normgb;
      dirgb[1] /= normgb;
      dirgb[2] /= normgb;
      //std::cout<<"Starting dir GAr (xyz): "<<dir[0]<<" "<<dir[1]<<" "<<dir[2]<<std::endl<<std::endl;
      if(NTPCClusters>100) std::cout<<"Ending par GAr (x,y,z,q/r,sinphi,tanlambda): "<<tparend[5]<<" "<<tparend[0]<<" "<<tparend[1]<<" "<<tparend[2]<<" "<<TMath::Sin(tparend[3])<<" "<<TMath::Tan(tparend[4])<<"\n\n";
      size_t nTPCClusters=0;
      if (TPCClusters.size()>unused_TPCClusters.size())
        { nTPCClusters = TPCClusters.size()-unused_TPCClusters.size(); }
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

      trackions.setData(dSigdXs_FWD,dSigdXs_BAK);
      tracktraj.setData(trajpts_FWD,trajpts_BAK);

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------

    // KalmanFit does a forwards or backwards Kalman fit using the sorted TPCCluster list
    // variables:  x is the independent variable
    // 0: y
    // 1: z
    // 2: curvature
    // 3: phi
    // 4: lambda

    int tpctrackfit2::KalmanFit( std::vector<TPCCluster> &TPCClusters,
                                 std::vector<int> &TPCClusterlist,    // sort ordered list
                                 std::vector<float> &trackparatend,
                                 float &chisquared,
                                 float &length,
                                 float *covmat,                     // 5x5 covariance matrix
                                 std::set<int> &unused_TPCClusters,
                                 std::vector<std::pair<float,float>>& dSigdXs,
                                 std::vector<TVector3>& trajpts)
    {

      // set some default values in case we return early

      size_t nTPCClusters = TPCClusterlist.size();
      chisquared = 0;
      length = 0;
      for (size_t i=0; i<5; ++i) trackparatend[i] = 0;
      for (size_t i=0; i<25; ++i) covmat[i] = 0;

      float roadsq = fRoadYZinFit*fRoadYZinFit;

      // estimate curvature, lambda, phi, xpos from the initial track parameters
      float curvature_init=0.1;
      float phi_init = 0;
      float lambda_init = 0;
      float xpos_init=0;
      float ypos_init=0;
      float zpos_init=0;
      float x_other_end = 0;
      if ( gar::rec::initial_trackpar_estimate(TPCClusters,
                                               TPCClusterlist,
                                               curvature_init,
                                               lambda_init,
                                               phi_init,
                                               xpos_init,
                                               ypos_init,
                                               zpos_init,
                                               x_other_end,
                                               fInitialTPNTPCClusters,
                                               fPrintLevel) != 0)
        {
          //std::cout << "kalman fit failed on initial trackpar estimate" << std::endl;
          return 1;
        }

      // Kalman fitter variables

      float xpos = xpos_init;

      TMatrixF P(5,5);  // covariance matrix of parameters
      // fill in initial guesses -- generous uncertainties on first value.
      P.Zero();
      P[0][0] = TMath::Sq(1);   // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1);   // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5);  // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // lambda uncertainty

      TMatrixF PPred(5,5);

      // per-step additions to the covariance matrix
      TMatrixF Q(5,5);
      Q.Zero();
      Q[2][2] = fKalCurvStepUncSq;     // allow for some curvature uncertainty between points
      Q[3][3] = fKalPhiStepUncSq;      // phi
      Q[4][4] = fKalLambdaStepUncSq;   // lambda

      // Noise covariance on the measured points.
      // 16 cm2 initially, might reasonably be lowered to typicalResidual near line 552-67
      TMatrixF R(2,2);
      R.Zero();
      R[0][0] = TMath::Sq(fKalCovZYMeasure);  // in cm^2
      R[1][1] = TMath::Sq(fKalCovZYMeasure);  // in cm^2

      // add the TPCClusters and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
      // scattering and energy loss can change the track parameters along the way.

      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
      TMatrixF F(5,5);
      TMatrixF FT(5,5);
      TVectorF parvec(5);
      parvec[0] = ypos_init;
      parvec[1] = zpos_init;
      parvec[2] = curvature_init;
      parvec[3] = phi_init;
      parvec[4] = lambda_init;
      TVectorF predstep(5);

      TMatrixF H(2,5);   // partial(obs)/partial(params)
      H.Zero();
      H[0][0] = 1;  // y
      H[1][1] = 1;  // z
      TMatrixF HT(5,2);

      TVectorF z(2);
      TVectorF ytilde(2);
      TVectorF hx(2);
      TMatrixF S(2,2);
      TMatrixF K(5,2);

      TMatrixF I(5,5);
      I.Zero();
      for (int i=0;i<5;++i) I[i][i] = 1;

      for (size_t iTPCCluster=1; iTPCCluster<nTPCClusters; ++iTPCCluster)
        {

          float xh = TPCClusters[TPCClusterlist[iTPCCluster]].Position()[0];
          float yh = TPCClusters[TPCClusterlist[iTPCCluster]].Position()[1];
          float zh = TPCClusters[TPCClusterlist[iTPCCluster]].Position()[2];

          if (fPrintLevel > 0)
            {
              std::cout << std::endl;
              std::cout << "Adding a new TPCCluster: " << xh << " " << yh << " " << zh << std::endl;
            }

          // for readability
          float curvature = parvec[2];
          float phi = parvec[3];
          float lambda = parvec[4];

          // update prediction to the plane containing x.  Maybe we need to find
          // the closest point on the helix to the TPCCluster we are adding,
          // and not necessarily force it to be at this x

          F.Zero();

          // y = yold + slope*dx*Sin(phi).   F[0][i] = dy/dtrackpar[i], where f is the update function slope*dx*Sin(phi)

          float slope = TMath::Tan(lambda);
          if (slope != 0)
            {
              slope = 1.0/slope;
            }
          else
            {
              slope = 1E9;
            }

          // relocate dx to be the location along the helix which minimizes
		  // [ (Xhit -Xhelix)/sigmaX ]^2 + [ (Yhit -Yhelix)/sigmaY ]^2 + [ (Zhit -Zhelix)/sigmaZ ]^2
          // Linearize for now near xpos:
          //        x = xpos + dx
          //        y = parvec[0] + slope * dx * sin(phi)
          //        z = parvec[1] + slope * dx * cos(phi)
          // parvec is updated as the fit progresses so the 'zero point' where y_0, z_0, phi_0
          // are defined is at the end of the fit, not at the place where the fit begins.
          //
          // old calc was just based on TPCCluster position in x:
          // float dx = xh - xpos;


          float dxnum = (slope/(fTPCClusterResolYZ*fTPCClusterResolYZ))*( (yh - parvec[0])*TMath::Sin(phi) + (zh - parvec[1])*TMath::Cos(phi) )
            + (xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX);
          float dxdenom = slope*slope/(fTPCClusterResolYZ*fTPCClusterResolYZ) + 1.0/(fTPCClusterResolX*fTPCClusterResolX);
          float dx = dxnum/dxdenom;
          if (dx == 0) dx = 1E-3;
          //std::cout << "dxdenom, dxnum: " << dxdenom << " " << dxnum << std::endl;
          //std::cout << "Track pos: " << xpos << " " << parvec[0] << " " << parvec[1] << " " << " TPCCluster pos: " << xh << " " << yh << " " << zh << std::endl;
          //std::cout << "dx old and new: " << xh - xpos << " " << dx << std::endl;


          //TODO check this -- are these the derivatives?
          // y = yold + dx*slope*TMath::Sin(phi)
          // slope = cot(lambda), so dslope/dlambda = -csc^2(lambda) = -1 - slope^2
          F[0][0] = 1.;
          F[0][3] = dx*slope*TMath::Cos(phi);
          F[0][4] = dx*TMath::Sin(phi)*(-1.0-slope*slope);

          // z = zold + slope*dx*Cos(phi)
          F[1][1] = 1.;
          F[1][3] = -dx*slope*TMath::Sin(phi);
          F[1][4] = dx*TMath::Cos(phi)*(-1.0-slope*slope);

          // curvature = old curvature -- doesn't change but put in an uncertainty
          F[2][2] = 1.;

          // phi = old phi + curvature*slope*dx
          // need to take the derivative of a product here
          F[3][2] = dx*slope;
          F[3][3] = 1.;
          F[3][4] = dx*curvature*(-1.0-slope*slope);

          // lambda -- same -- but put in an uncertainty in case it changes
          F[4][4] = 1.;

          // predicted step
          if (fPrintLevel > 1)
            {
              std::cout << "F Matrix: " << std::endl;
              F.Print();
              std::cout << "P Matrix: " << std::endl;
              P.Print();
            }
          if (fPrintLevel > 0)
            {
              std::cout << "x: " << xpos << " dx: " << dx <<  std::endl;
              std::cout << " Parvec:   y " << parvec[0] << " z " << parvec[1] << " c " << parvec[2] << " phi " << parvec[3] << " lambda " << parvec[4] << std::endl;
            }

          predstep = parvec;
          predstep[0] += slope*dx*TMath::Sin(phi);  // update y
          predstep[1] += slope*dx*TMath::Cos(phi);  // update z
          predstep[3] += slope*dx*curvature;        // update phi

          if (fPrintLevel > 1)
            {
              std::cout << " Predstep: y " << predstep[0] << " z " << predstep[1] << " c " << predstep[2] << " phi " << predstep[3] << " lambda " << predstep[4] << std::endl;
            }
          // equations from the extended Kalman filter
          FT.Transpose(F);
          PPred = F*P*FT + Q;
          if (fPrintLevel > 1)
            {
              std::cout << "PPred Matrix: " << std::endl;
              PPred.Print();
            }

          ytilde[0] = yh - predstep[0];
          ytilde[1] = zh - predstep[1];
          float ydistsq = ytilde.Norm2Sqr();
          if (ydistsq > roadsq)
            {
              unused_TPCClusters.insert(iTPCCluster);
              continue;
            }

          // Now the chisquared calculation.  Need the angle of track re line from axis of
          // TPC, that is alpha; and to determine which part of the detector we are
          // in.  The residual is determined based on a simple linear parameterization
          float impactAngle;
          TVector3 trajPerp(0.0, predstep[0],predstep[1]);
          float rTrj = trajPerp.Mag();
          TVector3 trajStepPerp(0.0, sin(predstep[3]),cos(predstep[3]));
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
          if (fPrintLevel > 0)
            {
              std::cout << "ytilde (residuals): " << std::endl;
              ytilde.Print();
            }
          if (fPrintLevel > 1)
            {
              std::cout << "H Matrix: " << std::endl;
              H.Print();
            }

          HT.Transpose(H);
          S = H*PPred*HT + R;
          if (fPrintLevel > 1)
            {
              std::cout << "S Matrix: " << std::endl;
              S.Print();
            }

          S.Invert();
          if (fPrintLevel > 1)
            {
              std::cout << "Inverted S Matrix: " << std::endl;
              S.Print();
            }

          K = PPred*HT*S;
          if (fPrintLevel > 1)
            {
              std::cout << "K Matrix: " << std::endl;
              K.Print();
            }

          float yprev = parvec[0];
          float zprev = parvec[1];
          parvec = predstep + K*ytilde;
          P = (I-K*H)*PPred;
          xpos = xpos + dx;
          //std::cout << " Updated xpos: " << xpos << " " << dx << std::endl;
          trajpts.emplace_back(xpos,parvec[0],parvec[1]);

          float d_length = TMath::Sqrt( dx*dx + TMath::Sq(parvec[0]-yprev) + TMath::Sq(parvec[1]-zprev) );
          length += d_length;

          // Save the ionization data - skip large gaps from sector boundaries
          float valSig = TPCClusters[TPCClusterlist[iTPCCluster]].Signal();
          if (d_length < fMinIonizGapCut)
            {
              std::pair pushme = std::make_pair(valSig,d_length);
              dSigdXs.push_back( pushme );
            }
          else
            // Have to remove the fellow before the large gap, too
            {
              if (dSigdXs.size()>0) dSigdXs.pop_back();
            }
        }

      for (size_t i=0; i<5; ++i)
        {
          trackparatend[i] = parvec[i];
        }
      trackparatend[5] = xpos;  // tack this on so we can specify where the track endpoint is
      if (fPrintLevel > 1)
        {
          std::cout << "Track params at end (y, z, curv, phi, lambda) " << trackparatend[0] << " " << trackparatend[1] << " " <<
            trackparatend[2] << " " << trackparatend[3] <<" " << trackparatend[4] << std::endl;
          S.Print();
        }

      // just for visualization of the initial track parameter guesses.
      //  Comment out when fitting tracks.
      //trackparatend[0] = ypos_init;
      //trackparatend[1] = zpos_init;
      //trackparatend[2] = curvature_init;
      //trackparatend[3] = phi_init;
      //trackparatend[4] = lambda_init;
      //trackparatend[5] = xpos_init;


      size_t icov=0;
      for (size_t i=0; i<5; ++i)
        {
          for (size_t j=0; j<5; ++j)
            {
              covmat[icov] = P[i][j];
            }
        }

      return 0;
    }

    //--------------------------------------------------------------------------------------------------------------


    DEFINE_ART_MODULE(tpctrackfit2)

  } // namespace rec
} // namespace gar
