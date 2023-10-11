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
      // variables: z is independent variable with rotating fram in radial coordinates
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
      int PDGcode = 211;

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
      
      particle_f.reconstructParticleFullOut(geom,PDGcode,10000); ///outwards reconstruction

      ///////////////Converting to garsoft convention
      std::vector<float> tparend(6,0);
      float covmatend[25];
      float chisqforwards = 0;
      float lengthforwards = 0;
      std::set<int> unused_TPCClusters;
      std::vector<std::pair<float,float>> dSigdXs_FWD;
      std::vector<TVector3> trajpts_FWD;


      Double_t xyz_end_out[3];
      particle_f.fParamOut[particle_f.fParamOut.size()-1].GetXYZ(xyz_end_out);
      Double_t ca=TMath::Cos(-particle_f.fParamOut[particle_f.fParamOut.size()-1].GetAlpha()), sa=TMath::Sin(-particle_f.fParamOut[particle_f.fParamOut.size()-1].GetAlpha());
      Double_t sf=particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[2];
      Double_t cf=TMath::Sqrt((1.- sf)*(1.+sf));
      Double_t sfrot = sf*ca - cf*sa;
      tparend[0]=xyz_end_out[1]+(GArCenter[1]-ystart);
      tparend[1]=xyz_end_out[0]+(GArCenter[2]-xstart);
      tparend[2]=particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[4]*(5*0.299792458e-3);
      tparend[3]=TMath::ASin(sfrot);
      tparend[4]=TMath::ATan(particle_f.fParamOut[particle_f.fParamOut.size()-1].GetParameter()[3]);
      tparend[5]=xyz_end_out[2];
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
    
      particle_b.reconstructParticleFullOut(geom,PDGcode,10000);

      /////////////Converting to garsoft convention
      std::vector<float> tparbeg(6,0);
      float covmatbeg[25];
      float chisqbackwards = 0;
      float lengthbackwards = 0;
      std::vector<std::pair<float,float>> dSigdXs_BAK;
      std::vector<TVector3> trajpts_BAK;
      
      Double_t xyz_start_out[3];
      particle_b.fParamOut[particle_b.fParamOut.size()-1].GetXYZ(xyz_start_out);
      Double_t cb=TMath::Cos(-particle_b.fParamOut[particle_b.fParamOut.size()-1].GetAlpha()), sb=TMath::Sin(-particle_b.fParamOut[particle_b.fParamOut.size()-1].GetAlpha());
      Double_t sfb=particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[2];
      Double_t cfb=TMath::Sqrt((1.- sf)*(1.+sf));
      Double_t sfrotb = sfb*cb - cfb*sb;

      tparbeg[0]=xyz_start_out[1]+(GArCenter[1]-yend);
      tparbeg[1]=xyz_start_out[0]+(GArCenter[2]-xend);
      tparbeg[2]=particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[4]*(5*0.299792458e-3);
      tparbeg[3]=TMath::ASin(sfrotb);
      tparbeg[4]=TMath::ATan(particle_b.fParamOut[particle_b.fParamOut.size()-1].GetParameter()[3]);
      tparbeg[5]=xyz_start_out[2];
      ////////////////////////////////////////////////////////////////////////////////////////////////////
      


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


    DEFINE_ART_MODULE(tpctrackfit2)

  } // namespace rec
} // namespace gar
