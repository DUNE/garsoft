// TrackPar class for defining tracks in garsoft and providing algorithms separate from the Track data product
// Tom Junk, July 2018
// change from slope parameter to lambda December 6, 2018

#include "Reco/TrackPar.h"
#include "TMath.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Geometry/Geometry.h"

namespace gar
{
  namespace rec
  {

    TrackPar::TrackPar(Track const &t, bool reversed)  // constructor from Track data product
    {
      fNTPCClusters = t.NHits();
      if (!reversed)
	{
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
	}
      else
	{
	  fXEnd = t.Vertex()[0];
	  fXBeg = t.End()[0];
	  for (size_t i=0; i<5; ++i)
	    {
	      fTrackParametersEnd[i] = t.TrackParBeg()[i];
	      fTrackParametersBegin[i] = t.TrackParEnd()[i];
	    }
	  fChisquaredBackwards = t.ChisqForward();
	  fChisquaredForwards = t.ChisqBackward();
	  fLengthBackwards = t.LengthForward();
	  fLengthForwards = t.LengthBackward();
	  t.CovMatEndSymmetric(fCovMatBeg);
	  t.CovMatEndSymmetric(fCovMatEnd);
	}
      CalcCenter();
    }

    //--------------------------------------------------------------------

    TrackPar::TrackPar(const float lengthforwards,  // constructor from parameters
		       const float lengthbackwards,
		       const size_t nTPCClusters,
		       const float xbeg,           // x location at beginning of track in cm
		       const float *trackparbeg,   // y, z, curvature, phi, lambda  -- 5-parameter track  (cm, cm, cm-1, radians, radians)
		       const float *covmatbeg,     // covariance matrix at beginning of track -- symmetric 5x5
		       const float chisqforward,   // chisquared of forwards fit
		       const float xend,           // x location at end of track
		       const float *trackparend,   // y, z, curvature, phi, lambda  -- 5-parameter track (cm, cm, cm-1, radians, radians)
		       const float *covmatend,     // covariance matrix at beginning of track -- symmetric 5x5
		       const float chisqbackward,  // chisquared of backwards fit
		       const double time) // timestamp
    {
      fNTPCClusters = nTPCClusters;
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

    const float *TrackPar::getTrackParametersBegin() const
    {
      return fTrackParametersBegin;
    }

    const float *TrackPar::getTrackParametersEnd() const
    {
      return fTrackParametersEnd;
    }

    const float *TrackPar::getCovMatBeg() const
    {
      return fCovMatBeg;
    }

    const float *TrackPar::getCovMatEnd() const
    {
      return fCovMatEnd;
    }


    const float TrackPar::getYCentBeg() const
    {
      return fYCentBeg;
    }

    const float TrackPar::getZCentBeg() const
    {
      return fZCentBeg;
    }

    const float TrackPar::getYCentEnd() const
    {
      return fYCentEnd;
    }

    const float TrackPar::getZCentEnd() const
    {
      return fZCentEnd;
    }

    const bool  TrackPar::getBegCentValid() const
    {
      return fBegCentValid;
    }

    const bool  TrackPar::getEndCentValid() const
    {
      return fEndCentValid;
    }

    size_t TrackPar::getNTPCClusters() const
    {
      return fNTPCClusters;
    }

    float TrackPar::getLengthForwards() const
    {
      return fLengthForwards;
    }

    float TrackPar::getLengthBackwards() const
    {
      return fLengthBackwards;
    }

    float TrackPar::getChisqForwards() const
    {
      return fChisquaredForwards;
    }

    float TrackPar::getChisqBackwards() const
    {
      return fChisquaredBackwards;
    }

    float TrackPar::getXBeg() const
    {
      return fXBeg;
    }

    float TrackPar::getXEnd() const
    {
      return fXEnd;
    }

    double TrackPar::getTime() const
    {
      return fTime;
    }

    //int TrackPar::getChargeBeg() const
    // Not built yet.  Use Track::ChargeBeg() for now
    //int TrackPar::getChargeEnd() const
    // Not built yet.  Use Track::ChargeEnd() for now

    void TrackPar::setNTPCClusters(const size_t nTPCClusters)
    {
      fNTPCClusters = nTPCClusters;
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

    void TrackPar::setTime(const double time)
    {
      fTime = time;
    }


    //--------------------------------------------------------------------

    // We couldn't make a constructor of Track that takes TrackPar, but we can
    // create a track data product from what we have here.

    gar::rec::Track TrackPar::CreateTrack()
    {

      return gar::rec::Track(fLengthForwards,
			     fLengthBackwards,
			     fNTPCClusters,
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

    TVector3 TrackPar::getXYZBeg() const
    {
      TVector3 result(getXBeg(),getTrackParametersBegin()[0],getTrackParametersBegin()[1]);
      return result;
    }

    TVector3 TrackPar::getXYZEnd() const
    {
      TVector3 result(getXEnd(),getTrackParametersEnd()[0],getTrackParametersEnd()[1]);
      return result;
    }


  }  // namespace rec

} // namespace gar
