////////////////////////////////////////////////////////////////////////
// Class:       tracker1
// Plugin Type: producer (art v2_11_02)
// File:        tracker1_module.cc
//
// Generated at Thu Jun 14 15:47:04 2018 by Thomas Junk using cetskelgen
// from cetlib version v3_03_01.
// Sort hits in X and add nearby hits.
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
#include "cetlib_except/exception.h"

#include <memory>

#include "TVectorF.h"
#include "TMatrix.h"

#include "Geant4/G4ThreeVector.hh"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/MagneticField/MagneticField.h"

// GArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "ReconstructionDataProducts/Hit.h"
#include "ReconstructionDataProducts/Track.h"

namespace gar {
  namespace rec {

    class tracker1 : public art::EDProducer {
    public:
      explicit tracker1(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      tracker1(tracker1 const &) = delete;
      tracker1(tracker1 &&) = delete;
      tracker1 & operator = (tracker1 const &) = delete;
      tracker1 & operator = (tracker1 &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;

    private:

      float fHitResolYZ;           ///< resolution in cm of a hit in YZ (pad size)
      float fHitResolX;            ///< resolution in cm of a hit in X (drift direction)
      float fSigmaRoad;            ///< how many sigma away from a track a hit can be and still add it
      float fXGapToEndTrack;       ///< how big a gap must be before we end a track and start a new one
      unsigned int fMinNumHits;    ///< minimum number of hits to define a track
      std::string fHitLabel;       ///< label of module creating hits
    };


    tracker1::tracker1(fhicl::ParameterSet const & p)
    {

      produces< std::vector<rec::Track> >();
      produces< ::art::Assns<rec::Hit, rec::Track> >();

      fHitResolYZ     = p.get<float>("HitResolYZ",0.7);
      fHitResolX      = p.get<float>("HitResolX",0.5);  // this is probably much better
      fSigmaRoad      = p.get<float>("SigmaRoad",5.0);
      fMinNumHits     = p.get<unsigned int>("MinNumHits",10);
      fHitLabel       = p.get<std::string>("HitLabel","hit");
    }

    void tracker1::produce(art::Event & e)
    {
      std::unique_ptr< std::vector<rec::Track> > trkCol(new std::vector<rec::Track>);
      std::unique_ptr< ::art::Assns<rec::Hit,rec::Track> > hitTrkAssns(new ::art::Assns<rec::Hit,rec::Track>);

      auto hitHandle = e.getValidHandle< std::vector<Hit> >(fHitLabel);
      auto const& hits = *hitHandle;

      art::ServiceHandle<mag::MagneticField> magFieldService;
      G4ThreeVector zerovec(0,0,0);
      G4ThreeVector magfield = magFieldService->FieldAtPoint(zerovec);

      // make an array of hit indices sorted by hit X position
      std::vector<float> hitx;
      for (size_t i=0; i<hits.size(); ++i)
	{
	  hitx.push_back(hits[i].Position()[0]);
	}
      std::vector<int> hsi(hitx.size());
      TMath::Sort((int) hitx.size(),hitx.data(),hsi.data());

      float roadsq = fSigmaRoad*fSigmaRoad;

      // record which hits we have assigned to which tracks
      std::vector< std::vector<int> > hitlist;
      // rough bounding box for the track -- to be used to get a guess of track parameters
      std::vector<float> tminx;
      std::vector<float> tmaxx;
      std::vector<float> tminy;
      std::vector<float> tmaxy;
      std::vector<float> tminz;
      std::vector<float> tmaxz;

      float resolSq = fHitResolYZ*fHitResolYZ;

      // idea -- find all of the track candidates a hit can be added to that have signfs <= roadsq, and
      // pick the one that is closest to the local linear extrapolation of the existing hits.

      for (size_t i=0; i<hits.size(); ++i)
	{
	  const float *hpos = hits[hsi[i]].Position();
	  float bestsignifs = -1;
	  int ibest = -1;
	  for (size_t itcand = 0; itcand < hitlist.size(); ++itcand)
	    {
	      const float *cpos = hits[hsi[hitlist[itcand].back()]].Position();
	      float signifs = TMath::Sq( (hpos[0]-cpos[0])/fHitResolX ) + 
		(TMath::Sq( (hpos[1]-cpos[1]) ) +
		 TMath::Sq( (hpos[2]-cpos[2]) ))/resolSq;
	      if (bestsignifs < 0 || signifs < bestsignifs)
		{
		  bestsignifs = signifs;
		  ibest = itcand;
		}
	    }
	  if (ibest == -1 || bestsignifs > roadsq)  // start a new track if we're not on the road, or if we had no tracks to begin with
	    {
	      std::vector<int> hloc;
	      hloc.push_back(i);
	      hitlist.push_back(hloc);
	      tminx.push_back(hpos[0]);
	      tmaxx.push_back(hpos[0]);
	      tminy.push_back(hpos[1]);
	      tmaxy.push_back(hpos[1]);
	      tminz.push_back(hpos[2]);
	      tmaxz.push_back(hpos[2]);
	    }
	  else  // add the hit to the existing best track
	    {
	      hitlist[ibest].push_back(i);
	      tminx[ibest] = TMath::Min(tminx[ibest],hpos[0]);
	      tmaxx[ibest] = TMath::Max(tmaxx[ibest],hpos[0]);
	      tminy[ibest] = TMath::Min(tminy[ibest],hpos[1]);
	      tmaxy[ibest] = TMath::Max(tmaxy[ibest],hpos[1]);
	      tminz[ibest] = TMath::Min(tminz[ibest],hpos[2]);
	      tmaxz[ibest] = TMath::Max(tmaxz[ibest],hpos[2]);
	    }
	}

      // now that we have the hits assigned, fit the tracks and adjust the hit lists.
      // fit them in both directions.  The Kalman filter gives the most precise measurement
      // of position and momentum at the end.

      size_t ntracks = hitlist.size();
      for (size_t itrack=0; itrack<ntracks; ++itrack)
	{
	  size_t nhits = hitlist[itrack].size();
	  if ( nhits >= fMinNumHits)
	    {

	      // variables:  x is the independent variable
	      // 0: y
	      // 1: z
	      // 2: curvature
	      // 3: phi
	      // 4: slope

	      // initial track position.  TODO -- start the track off with better guesses for the curvature, slope, and phi

	      float trackbeg[3] = {hits[hsi[hitlist[itrack][0]]].Position()[0],hits[hsi[hitlist[itrack][0]]].Position()[1],hits[hsi[hitlist[itrack][0]]].Position()[2]};
	      float curvature_init = 0.1; // initial guess -- maybe use the min and max above to refine this guess
	      float phi_init = 0;        // initial guess at the track beginning
	      float slope_init = 1;   // d(yz distance)/dx

	      float xpos = trackbeg[0];

	      TMatrixF P(5,5);  // covariance matrix of parameters
	      // fill in initial guesses -- generous uncertainties on first value.
	      P.Zero();
	      P[0][0] = TMath::Sq(1); // initial position uncertainties -- y
	      P[1][1] = TMath::Sq(1); // and z
	      P[2][2] = TMath::Sq(1);  // curvature of zero gets us to infinite momentum, and curvature of 2 is curled up tighter than the pads
	      P[3][3] = TMath::Sq(2*TMath::Pi()); // phi uncertainty
	      P[4][4] = TMath::Sq(2.0);  // slope uncertainty
	  
	      TMatrixF PPred(5,5);

	      // per-step additions to the covariance matrix
	      TMatrixF Q(5,5);
	      Q.Zero();
	      Q[2][2] = 0.0001;  // allow for some curvature uncertainty between points  (1%)
	      Q[3][3] = 0.0001;   // phi
	      Q[4][4] = 0.0001;   // slope

	      // uncertainties on the measured points  (big for now)
	      TMatrixF R(2,2);
	      R[0][0] = 0.5;  // in cm^2
	      R[1][1] = 0.5;  // in cm^2

	      // add the hits and update the track parameters and uncertainties.  Put in additional terms to keep uncertainties from shrinking when
              // scattering and energy loss can change the track parameters along the way.

	      // F = partial(updatefunc)/partial(parvec).  Update functions are in the comments below.
	      TMatrixF F(5,5);
	      TMatrixF FT(5,5);
	      TVectorF parvec(5);
	      parvec[0] = trackbeg[1];
	      parvec[1] = trackbeg[2];
	      parvec[2] = curvature_init;
	      parvec[3] = phi_init;
	      parvec[4] = slope_init;
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
	      for (int i=0;i<5;++i) I[i][i] = i;

	      for (size_t ihit=1; ihit<nhits; ++ihit)
		{
	          float xh = hits[hsi[hitlist[itrack][ihit]]].Position()[0];
		  float yh = hits[hsi[hitlist[itrack][ihit]]].Position()[1];
		  float zh = hits[hsi[hitlist[itrack][ihit]]].Position()[2];

		  // for readability

		  float curvature = parvec[2];
		  float phi = parvec[3];
		  float slope = parvec[4];
		  float dx = xh - xpos;

		  // update prediction to the plane containing x

		  F.Zero();

		  // y = yold + slope*dx*Sin(phi)
		  F[0][0] = 1.; 
		  F[0][3] = dx*slope*TMath::Cos(phi);
		  F[0][4] = dx*TMath::Sin(phi);

		  // z = zold + slope*dx*Cos(phi)
		  F[1][1] = 1.;
		  F[1][3] = -dx*slope*TMath::Sin(phi);
		  F[1][4] = dx*TMath::Cos(phi);

		  // curvature = old curvature -- doesn't change but put in an uncertainty
		  F[2][2] = 1.;

		  // phi = old phi + curvature*slope*dx  
		  // need to take the derivative of a product here
		  F[3][2] = dx*slope;
		  F[3][3] = 1.;
		  F[3][4] = dx*curvature;

		  // slope -- same -- but put in an uncertainty in case it changes
		  F[4][4] = 1.;

		  // predicted step

		  predstep = parvec;
		  predstep[0] += slope*dx*TMath::Sin(phi);  // update y
		  predstep[1] += slope*dx*TMath::Cos(phi);  // update z
		  predstep[3] += slope*dx*curvature;        // update phi

		  // equations from the extended Kalman filter
		  FT.Transpose(F);
		  PPred = F*P*FT + Q;
		  ytilde[0] = yh - predstep[0];
		  ytilde[1] = zh - predstep[1];
		  HT.Transpose(H);
		  S = H*PPred*HT + R;
		  S.Invert();
		  K = PPred*HT*S;
		  parvec = predstep + K*ytilde;
		  P = (I-K*H)*PPred;
		  xpos = xh;
		}
	    }
	}

      e.put(std::move(trkCol));
      e.put(std::move(hitTrkAssns));
    }

    DEFINE_ART_MODULE(tracker1)

  } // namespace rec
} // namespace gar