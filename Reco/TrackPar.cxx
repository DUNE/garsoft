// TrackPar class for defining tracks in garsoft and providing algorithms separate from the Track data product
// Tom Junk, July 2018 
// change from slope parameter to lambda December 6, 2018

#include "Reco/TrackPar.h"
#include "TMath.h"

namespace gar
{
  namespace rec
  {

    TrackPar::TrackPar(Track const &t)  // constructor from Track data product
    {
      fNHits = t.NHits();
      fXBeg = t.Vertex()[0];
      fXEnd = t.End()[0];
      for (size_t i=0; i<5; ++i)
	{
	  fTrackParametersBegin[i] = t.TrackParBeg()[i];
	  fTrackParametersEnd[i] = t.TrackParEnd()[i];
	}
      fChisquaredForwards = t.ChisqForward();
      fChisquaredBackwards = t.ChisqBackward();
      fLengthForwards = t.LengthForward();
      fLengthBackwards = t.LengthBackward();
      t.CovMatBegSymmetric(fCovMatBeg);
      t.CovMatEndSymmetric(fCovMatEnd);
      CalcCenter();
    }

    //--------------------------------------------------------------------
    
    TrackPar::TrackPar(const float lengthforwards,  // constructor from parameters
		       const float lengthbackwards,
		       const size_t nhits,
		       const float xbeg,           // x location at beginning of track in cm
		       const float *trackparbeg,   // y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)
		       const float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
		       const float chisqforward,   // chisquared of forwards fit
		       const float xend,           // x location at end of track
		       const float *trackparend,   // y, z, curvature, phi, lambda  -- 5-parameter track (cm, cm, cm-1, radians, radians)
		       const float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
		       const float chisqbackward,  // chisquared of backwards fit
		       const ULong64_t time) // timestamp
    {
      fNHits = nhits;
      fLengthForwards = lengthforwards;
      fLengthBackwards = lengthbackwards;
      fXBeg = xbeg;
      fXEnd = xend;
      for (size_t i=0; i<5; ++i)
	{
	  fTrackParametersBegin[i] = trackparbeg[i];
	  fTrackParametersEnd[i] = trackparend[i];
	}
      for (size_t i=0; i<25; ++i)
	{
	  fCovMatBeg[i] = covmatbeg[i];
	  fCovMatEnd[i] = covmatend[i];
	}
      fChisquaredForwards = chisqforward;
      fChisquaredBackwards = chisqbackward;
      fTime = time;
      CalcCenter();
    }

    //--------------------------------------------------------------------

    void TrackPar::CalcCenter()
    {
      if (fTrackParametersBegin[2] != 0)
	{
	  float phi = fTrackParametersBegin[3];
	  float r = 1.0/fTrackParametersBegin[2];
	  fZCentBeg = fTrackParametersBegin[1] - r*TMath::Sin(phi);
	  fYCentBeg = fTrackParametersBegin[0] + r*TMath::Cos(phi);
	  fBegCentValid = true;
	}
      else
	{
	  fBegCentValid = false;
	}
      if (fTrackParametersEnd[2] != 0)
	{
	  float phi = fTrackParametersEnd[3];
	  float r = 1.0/fTrackParametersEnd[2];
	  fZCentEnd = fTrackParametersEnd[1] - r*TMath::Sin(phi);
	  fYCentEnd = fTrackParametersEnd[0] + r*TMath::Cos(phi);
	  fEndCentValid = true;
	}
      else
	{
	  fEndCentValid = false;
	}
    }

    //--------------------------------------------------------------------

    const float *TrackPar::getTrackParametersBegin()
    {
      return fTrackParametersBegin;
    }

    const float *TrackPar::getTrackParametersEnd()
    {
      return fTrackParametersEnd;
    }

    const float *TrackPar::getCovMatBeg()
    {
      return fCovMatBeg;
    }

    const float *TrackPar::getCovMatEnd()
    {
      return fCovMatEnd;
    }

    size_t TrackPar::getNHits()
    {
      return fNHits;
    }

    float TrackPar::getLengthForwards()
    {
      return fLengthForwards;
    }

    float TrackPar::getLengthBackwards()
    {
      return fLengthBackwards;
    }

    float TrackPar::getChisqForwards()
    {
      return fChisquaredForwards;
    }

    float TrackPar::getChisqBackwards()
    {
      return fChisquaredBackwards;
    }

    float TrackPar::getXBeg()
    {
      return fXBeg;
    }

    float TrackPar::getXEnd()
    {
      return fXEnd;
    }

    ULong64_t TrackPar::getTime()
    {
      return fTime;
    }

    int TrackPar::getChargeBeg()   // just returns +1 or -1 depending on the sign of the curvature at the track beginning point
    {
      int icharge=0;
      if (fTrackParametersBegin[2] > 0)
	{
	  icharge = 1;
	}
      else
	{
	  icharge = -1;
	}
      return icharge;
    }

    int TrackPar::getChargeEnd()  // just returns +1 or -1 depending on the sign of the curvature at the track ending point
    {
      int icharge=0;
      if (fTrackParametersEnd[2] > 0)
	{
	  icharge = 1;
	}
      else
	{
	  icharge = -1;
	}
      return icharge;
    }

    void TrackPar::setNHits(const size_t nhits)
    {
      fNHits = nhits;
    }

    void TrackPar::setTrackParametersBegin(const float *tparbeg)
    {
      for (size_t i=0; i<5; ++i)
	{
	  fTrackParametersBegin[i] = tparbeg[i];
	}
    }

    void TrackPar::setTrackParametersEnd(const float *tparend)
    {
      for (size_t i=0; i<5; ++i)
	{
	  fTrackParametersEnd[i] = tparend[i];
	}
    }

    void TrackPar::setCovMatBeg(const float *covmatbeg)
    {
      for (size_t i=0; i<25; ++i)
	{
	  fCovMatBeg[i] = covmatbeg[i];
	}
    }

    void TrackPar::setCovMatEnd(const float *covmatend)
    {
      for (size_t i=0; i<25; ++i)
	{
	  fCovMatEnd[i] = covmatend[i];
	}
    }

	void TrackPar::setHits(const std::vector<TVector3>* inHits)
	{
		for (size_t ii=0; ii<inHits->size(); ii++)
		{
			fTrackHits.push_back(inHits->at(ii));
		}
	}

    void TrackPar::setLengthForwards(const float lengthforwards)
    {
      fLengthForwards = lengthforwards;
    }

    void TrackPar::setLengthBackwards(const float lengthbackwards)
    {
      fLengthBackwards = lengthbackwards;
    }

    void TrackPar::setChisqForwards(const float chisqforwards)
    {
      fChisquaredForwards = chisqforwards;
    }

    void TrackPar::setChisqBackwards(const float chisqbackwards)
    {
      fChisquaredBackwards = chisqbackwards;
    }

    void TrackPar::setXBeg(const float xbeg)
    {
      fXBeg = xbeg;
    }

    void TrackPar::setXEnd(const float xend)
    {
      fXEnd = xend;
    }

    // distance from a point to a helix

    float TrackPar::DistXYZ(const float *xyz)
    {
      float dmin = TMath::Min(DistXYZ_Aux(xyz,true),
	        	      DistXYZ_Aux(xyz,false));
      return dmin;
    }

    void TrackPar::setTime(const ULong64_t time)
    {
      fTime = time;
    }

    // split off the calc so we can do the two distances, using the track begin parameters and
    // the track end parameters

    float TrackPar::DistXYZ_Aux(const float *xyz, bool useBegin)
    {
      float distance = 0;

      // just to make formulas more readable.  (xt,yt,zt) = test point.

      float xt = xyz[0];
      float yt = xyz[1];
      float zt = xyz[2];
      float x0 = fXBeg;
      float y0 = fTrackParametersBegin[0];
      float z0 = fTrackParametersBegin[1];
      float curv = fTrackParametersBegin[2];
      float phi0 = fTrackParametersBegin[3];
      float s = TMath::Tan(fTrackParametersBegin[4]);
      if (s != 0)
	{
	  s = 1.0/s;
	}
      else
	{
	  s = 1E9;
	}
      float yc = fYCentBeg;
      float zc = fZCentBeg;
      bool valid_center = fBegCentValid;
      if (!useBegin)
	{
	  x0 = fXEnd;
	  y0 = fTrackParametersEnd[0];
	  z0 = fTrackParametersEnd[1];
	  curv = fTrackParametersEnd[2];
	  phi0 = fTrackParametersEnd[3];
	  s = TMath::Tan(fTrackParametersEnd[4]);
	  if (s != 0)
	    {
	      s = 1.0/s;
	    }
	  else
	    {
	      s = 1E9;
	    }


	  yc = fYCentEnd;
	  zc = fZCentEnd;
	  valid_center = fEndCentValid;
	}

      float cosphi0 = TMath::Cos(phi0);
      float sinphi0 = TMath::Sin(phi0);

      if (curv == 0 || !valid_center )
	{
	  // distance from a point to a line -- use the norm of the cross product of the
	  // unit vector along the line and the displacement between the tesst point and 
	  // a point on the line

	  float h = TMath::Sqrt(1.0 + s*s);
	  float xhc = (yt-y0)*(cosphi0/h) - (sinphi0/h)*(zt-z0);
	  float yhc = (s/h)*(zt-z0) - (xt-x0)*(cosphi0/h);
	  float zhc = (xt-x0)*(sinphi0/h) - (s/h)*(yt-y0);

	  distance = TMath::Sqrt( xhc*xhc + yhc*yhc + zhc*zhc );
	}
      else if (s == 0)
	{
	  // zero slope.  The track is another line, this time along the x axis

	  distance = TMath::Sqrt( TMath::Sq(yt-y0) + TMath::Sq(zt-z0) );
	}
      else
	{
	  // general case -- need to compute distance from a point to a helix
	  float r = 1.0/curv;
	  float A = (xt-x0) + (r/s)*phi0;
	  float B = s*(yt-yc);
	  float C = s*(zt-zc);
	  float D = r/s;
	  float E = TMath::Sqrt(B*B + C*C);
	  float a = TMath::ATan2(-B,C);

	  float eps = 0.001;

	  float phicent = (s*curv)*(xt-x0);
	  float philow = phicent - TMath::Pi();
	  float phihi = phicent + TMath::Pi();

	  // solve for phi of the closest point using Newton's method
	  // with a fixed number of iterations.

	  float phi = phicent;
	  for (size_t iter=0; iter<10; ++iter)
	    {
	      float f = A + E*TMath::Cos(phi+a) - D*phi;
	      if (f == 0) 
		{
		  break;
		}
	      float fprime = -E*TMath::Sin(phi+a) -D;
	      if (fprime != 0)
		{
		  phi = phi - f/fprime;
		}
	      else
		{
		  // protect against zero derivative -- just scan a couple of points.
		  float f2 = A + E*TMath::Cos(phi+a+eps) - D*(phi+eps);
		  float f3 = A + E*TMath::Cos(phi+a-eps) - D*(phi-eps);
		  if (TMath::Abs(f2) < TMath::Abs(f3))
		    {
		      phi = phi + eps;
		    }
		  else
		    {
		      phi = phi - eps;
		    }
		}
	      if (phi<philow) 
		{
		  phi=philow;
		}
	      if (phi>phihi) 
		{
		  phi=phihi;
		}
	    }
	  
	  float xp = x0 + r*(phi-phi0)/s;
	  float yp = yc - r*TMath::Cos(phi);
	  float zp = zc + r*TMath::Sin(phi);

	  distance = TMath::Sqrt( TMath::Sq(xt-xp) + TMath::Sq(yt-yp) + TMath::Sq(zt-zp) );
	}

      return distance;

    }

    //--------------------------------------------------------------------

    // We couldn't make a constructor of Track that takes TrackPar, but we can create a track data product from what we have here.

    gar::rec::Track TrackPar::CreateTrack()
    {

      return gar::rec::Track(fLengthForwards,
			     fLengthBackwards,
			     fNHits,
			     fXBeg,
			     fTrackParametersBegin,
			     fCovMatBeg,
			     fChisquaredForwards,
			     fXEnd,
			     fTrackParametersEnd,
			     fCovMatEnd,
			     fChisquaredBackwards,
			     fTime
			     );
    }

	// for when we want hits too
    gar::rec::Track TrackPar::CreateTrackWithHits()
    {

		float hitXs[200];
		float hitYs[200];
		float hitZs[200];
		for(uint ii=0; ii<fTrackHits.size() && ii<200; ii++){
			hitXs[ii]=fTrackHits.at(ii).X();
			hitYs[ii]=fTrackHits.at(ii).Y();
			hitZs[ii]=fTrackHits.at(ii).Z();
		}

      return gar::rec::Track(fLengthForwards,
			     fLengthBackwards,
			     fNHits,
				 hitXs,
				 hitYs,
				 hitZs,
			     fXBeg,
			     fTrackParametersBegin,
			     fCovMatBeg,
			     fChisquaredForwards,
			     fXEnd,
			     fTrackParametersEnd,
			     fCovMatEnd,
			     fChisquaredBackwards,
			     fTime
			     );
    }


    TVector3 TrackPar::getXYZBeg()
    {
      TVector3 result(getXBeg(),getTrackParametersBegin()[0],getTrackParametersBegin()[1]);
      return result;
    }

    TVector3 TrackPar::getXYZEnd()
    {
      TVector3 result(getXEnd(),getTrackParametersEnd()[0],getTrackParametersEnd()[1]);
      return result;
    }

    //-----------------------------------

    // position along a track given X.

    TVector3 TrackPar::getPosAtX(const float x, bool usebegpar)
    {
      float y=0;
      float z=0;
      if (usebegpar)
	{
	  float stmp = TMath::Tan(fTrackParametersBegin[4]);
	  if (stmp != 0)
	    {
	      stmp = 1.0/stmp;
	    }
	  else
	    {
	      stmp = 1E9;
	    }
	  float phi = (x-fXBeg)*stmp*fTrackParametersBegin[2] + fTrackParametersBegin[3];   // (x-x0)*s*c + phi0  -- c=1/r
	  if (fTrackParametersBegin[2] == 0)
	    {
	      y = 0;
	      z = 0;
	    }
	  else
	    {
	      y = fYCentBeg - TMath::Cos(phi)/fTrackParametersBegin[2];
	      z = fZCentBeg + TMath::Sin(phi)/fTrackParametersBegin[2];
	    }	  
	}
      else
	{
	  float stmp = TMath::Tan(fTrackParametersEnd[4]);
	  if (stmp != 0)
	    {
	      stmp = 1.0/stmp;
	    }
	  else
	    {
	      stmp = 1E9;
	    }
	  float phi = (x-fXEnd)*stmp*fTrackParametersEnd[2] + fTrackParametersEnd[3];   // (x-x0)*s*c + phi0  -- c=1/r
	  if (fTrackParametersEnd[2] == 0)
	    {
	      y = 0;
	      z = 0;
	    }
	  else
	    {
	      y = fYCentEnd - TMath::Cos(phi)/fTrackParametersEnd[2];
	      z = fZCentEnd + TMath::Sin(phi)/fTrackParametersEnd[2];
	    }	  
	}
      TVector3 pos(x,y,z);
      return pos;
    }

  }  // namespace rec
} // namespace gar
