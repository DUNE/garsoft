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
#include "ReconstructionDataProducts/TPCCluster.h"
#include "ReconstructionDataProducts/VecHit.h"
#include "ReconstructionDataProducts/Track.h"
#include "Reco/TrackPar.h"
#include "Reco/tracker2algs.h"
#include "Geometry/Geometry.h"

#include "Geant4/G4ThreeVector.hh"

#include "nutools/MagneticField/MagneticField.h"

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

      std::string fPatRecLabel;     ///< input patrec tracks and associations
      int fPrintLevel;              ///< debug printout:  0: none, 1: just track parameters and residuals, 2: all
      float  fKalCurvStepUncSq;     ///< constant uncertainty term on each step of the Kalman fit -- squared, for curvature
      float  fKalPhiStepUncSq;      ///< constant uncertainty term on each step of the Kalman fit -- squared, for phi
      float  fKalLambdaStepUncSq;   ///< constant uncertainty term on each step of the Kalman fit -- squared, for lambda
      float  fTPCClusterResolYZ;           ///< resolution in cm of a TPCCluster in YZ (pad size)
      float  fTPCClusterResolX;            ///< resolution in cm of a TPCCluster in X (drift direction)
      unsigned int fInitialTPNTPCClusters; ///< number of TPCClusters to use for initial trackpar estimate, if present
      unsigned int fMinNumTPCClusters;     ///< minimum number of TPCClusters to define a track
      int fDumpTracks;              ///< 0: do not print out tracks, 1: print out tracks
      float fTPCClusterResolYZinFit;       ///< TPCCluster resolution parameter to use in fit
      float fRoadYZinFit;           ///< cut in cm for dropping TPCClusters from tracks in fit

      int KalmanFit( std::vector<TPCCluster> &TPCClusters,
		     std::vector<int> &TPCClusterlist,
		     std::vector<float> &trackparatend,
		     float &chisquared,
		     float &length,
		     float *covmat,    // 5x5 covariance matrix
		     std::set<int> &unused_TPCClusters);

      int KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
			    TrackPar &trackpar);

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
      fInitialTPNTPCClusters    = p.get<unsigned int>("InitialTPNTPCClusters",100);
      fDumpTracks        = p.get<int>("DumpTracks",0);
      fTPCClusterResolYZ        = p.get<float>("TPCClusterResolYZ",1.0); // TODO -- think about what this value is
      fTPCClusterResolX         = p.get<float>("TPCClusterResolX",0.5);  // this is probably much better
      fTPCClusterResolYZinFit   = p.get<float>("TPCClusterResolYZinFit",4.0);
      fRoadYZinFit       = p.get<float>("RoadYZinFit",1.0);

      art::InputTag patrecTag(fPatRecLabel);
      consumes< std::vector<gar::rec::Track> >(patrecTag);
      consumes< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >(patrecTag);

      // probably don't need the vector hits at this point if we have the TPCClusters
      //consumes< std::vector<gar::rec::VecHit> >(patrecTag);
      //consumes< art::Assns<gar::rec::VecHit, gar::rec::Track> >(patrecTag);

      produces< std::vector<gar::rec::Track> >();
      produces< art::Assns<gar::rec::TPCCluster, gar::rec::Track> >();
    }

    void tpctrackfit2::produce(art::Event& e)
    {
      // output collections

      std::unique_ptr< std::vector<gar::rec::Track> > trkCol(new std::vector<gar::rec::Track>);
      std::unique_ptr< art::Assns<gar::rec::TPCCluster,gar::rec::Track> > TPCClusterTrkAssns(new ::art::Assns<gar::rec::TPCCluster,gar::rec::Track>);

      // inputs

      auto patrecTrackHandle = e.getValidHandle< std::vector<gar::rec::Track> >(fPatRecLabel);
      auto const& patrecTracks = *patrecTrackHandle;

      auto const trackPtrMaker = art::PtrMaker<gar::rec::Track>(e);
      //auto const TPCClusterPtrMaker = art::PtrMaker<gar::rec::TPCCluster>(e, TPCClusterHandle.id());

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      art::ServiceHandle<geo::Geometry> geo;
      double xtpccent = geo->TPCXCent();
      double ytpccent = geo->TPCYCent();
      double ztpccent = geo->TPCZCent();
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
	  if (KalmanFitBothWays(TPCClusters,trackparams) == 0)   // to think about -- unused TPCClusters?  Or just ignore them in the fit?
	    {
	      trkCol->push_back(trackparams.CreateTrack());
	      auto const trackpointer = trackPtrMaker(trkCol->size()-1);
	      for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
		{
		  TPCClusterTrkAssns->addSingle(TPCClustersFromPatRecTracks.at(itrack).at(iTPCCluster),trackpointer);
		}
	    }
	}

      e.put(std::move(trkCol));
      e.put(std::move(TPCClusterTrkAssns));
    }

    int tpctrackfit2::KalmanFitBothWays(std::vector<gar::rec::TPCCluster> &TPCClusters,
				        TrackPar &trackpar)

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
      gar::rec::sort_TPCClusters_along_track(TPCClusters,hlf,hlb,fPrintLevel,lftmp,lbtmp);

      std::vector<float> tparend(6);
      float covmatend[25];
      float chisqforwards = 0;
      float lengthforwards = 0;
      std::set<int> unused_TPCClusters;

      int retcode = KalmanFit(TPCClusters,hlf,tparend,chisqforwards,lengthforwards,covmatend,unused_TPCClusters);
      if (retcode != 0) return 1;

      // the "backwards" fit is in decreasing x.  Track paramters are at the end of the fit, the other end of the track

      std::vector<float> tparbeg(6);
      float covmatbeg[25];
      float chisqbackwards = 0;
      float lengthbackwards = 0;

      retcode = KalmanFit(TPCClusters,hlb,tparbeg,chisqbackwards,lengthbackwards,covmatbeg,unused_TPCClusters);
      if (retcode != 0) return 1;

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
				 std::set<int> &unused_TPCClusters)
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
      P[0][0] = TMath::Sq(1); // initial position uncertainties -- y
      P[1][1] = TMath::Sq(1); // and z
      P[2][2] = TMath::Sq(.5);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
      P[3][3] = TMath::Sq(.5); // phi uncertainty
      P[4][4] = TMath::Sq(.5);  // lambda uncertainty

      TMatrixF PPred(5,5);

      // per-step additions to the covariance matrix
      TMatrixF Q(5,5);
      Q.Zero();
      Q[2][2] = fKalCurvStepUncSq;     // allow for some curvature uncertainty between points
      Q[3][3] = fKalPhiStepUncSq;      // phi
      Q[4][4] = fKalLambdaStepUncSq;   // lambda

      // uncertainties on the measured points  (big for now)
      TMatrixF R(2,2);
      R.Zero();
      R[0][0] = TMath::Sq(fTPCClusterResolYZinFit);  // in cm^2
      R[1][1] = TMath::Sq(fTPCClusterResolYZinFit);  // in cm^2

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

	  // update prediction to the plane containing x.  Maybe we need to find the closest point on the helix to the TPCCluster we are adding,
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

	  // relocate dx to be the location along the helix of the closest point.  Linearize for now near xpos.
	  // old calc

	  float dx = xh - xpos;

	  float dxdenom = slope*slope/(fTPCClusterResolYZ*fTPCClusterResolYZ) + 1.0/(fTPCClusterResolX*fTPCClusterResolX);
	  float dxnum = (slope/(fTPCClusterResolYZ*fTPCClusterResolYZ))*( (yh - parvec[0])*TMath::Sin(phi) + (zh - parvec[1])*TMath::Cos(phi) )
	    + (xh - xpos)/(fTPCClusterResolX*fTPCClusterResolX);
	  dx = dxnum/dxdenom;
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
	  chisquared += ytilde.Norm2Sqr()/TMath::Sq(fTPCClusterResolYZ);
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

	  length += TMath::Sqrt( dx*dx + TMath::Sq(parvec[0]-yprev) + TMath::Sq(parvec[1]-zprev) );
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

      // just for visualization of the initial track parameter guesses.  Comment out when fitting tracks

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
