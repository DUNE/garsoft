////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.cxx
/// \brief Interface to algorithm class for the standard, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov, trj@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapStandardAlg.h"

#include "TMath.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
  namespace geo{

    //----------------------------------------------------------------------------
    ChannelMapStandardAlg::ChannelMapStandardAlg(fhicl::ParameterSet const& p)
    {
    }

    //----------------------------------------------------------------------------

    void ChannelMapStandardAlg::Initialize(GeometryCore &geo)
    {
      // start over:
      Uninitialize();

      fROC = AliTPCROC::Instance();

      //std::string driftvolname = geo.GetGArTPCVolumeName();

      fTPCCenter.x = geo.TPCXCent();
      fTPCCenter.y = geo.TPCYCent();
      fTPCCenter.z = geo.TPCZCent();

      // get these from the geometry?
      // all dimensions are in cm

      fNumSectors = 18;
      fSectorOffsetAngleDeg = 10.0;
      fPhiSectorWidth = 2.0*TMath::Pi()/fNumSectors;
      fNumPadRowsIROC = fROC->GetNRowLow();
      fNumPadRowsOROCI = fROC->GetNRowUp1();
      fNumPadRowsOROCO = fROC->GetNRowUp2();
      fPadHeightIROC = fROC->GetInnerPadLength();
      fPadWidthIROC = fROC->GetInnerPadWidth();
      fPadWidthOROC = fROC->GetOuterPadWidth();
      fPadHeightOROCI = fROC->GetOuter1PadLength();
      fPadHeightOROCO = fROC->GetOuter2PadLength();

      fIROCInnerRadius = fROC->GetInnerRadiusLow() + 1.575;
      fIROCOuterRadius = fROC->GetInnerRadiusUp();
      fOROCInnerRadius = fROC->GetOuterRadiusLow()+1.6;
      fOROCPadHeightChangeRadius = 198.6;
      fOROCOuterRadius = fROC->GetOuterRadiusUp();

      // old hardcoded numbers
      //fOROCOuterRadius = 246.6;
      //fOROCInnerRadius = 134.6;
      //fIROCOuterRadius = 132.1;

      fGapChannelNumber = 10000000;

      fCenterPadWidth = 0.6;
      fXPlaneLoc = geo.TPCLength()/2.0;
      //std::cout << "Plane loc from geometry service " << fXPlaneLoc << std::endl;
      //float TsectorH = TMath::Tan(TMath::DegToRad()*(360/(fNumSectors*2)));

      // count up channels
      // get the xyz center for each pixel and fill pad row information -- nominal geometry, negative x side


      for (UInt_t isector = 0; isector < fNumSectors; ++isector)
	{
	  UInt_t ipadacc=0;

	  // Inner Readout Chamber

	  for (UInt_t irow = 0; irow < fNumPadRowsIROC; ++irow)
	    {
	      UInt_t numpads = fROC->GetNPads(isector,irow);
	      if (isector == 0)
		{
		  fNumPadsPerRow.push_back(numpads);
		  fFirstPadInRow.push_back(ipadacc);
		}
	      for (UInt_t ipad = 0; ipad < numpads; ++ipad)
		{
		  Float_t pos[3];
		  fROC->GetPositionGlobal(isector,irow,ipad,pos);
		  pos[2] = pos[0];
		  pos[0] = -fXPlaneLoc;
		  fPixelCenters.emplace_back(pos[0],pos[1],pos[2]);
		  ipadacc++;
		}
	    }

	  // Outer readout chamber

	  for (UInt_t irow = 0; irow < fNumPadRowsOROCI + fNumPadRowsOROCO; ++irow)
	    {

	      UInt_t numpads = fROC->GetNPads(isector + 36,irow);
	      if (isector == 0)
		{
		  fNumPadsPerRow.push_back(numpads);
		  fFirstPadInRow.push_back(ipadacc);
		}
	      for (UInt_t ipad = 0; ipad < numpads; ++ipad)
		{
		  Float_t pos[3];
		  fROC->GetPositionGlobal(isector + 36,irow,ipad,pos);
		  pos[2] = pos[0];
		  pos[0] = -fXPlaneLoc;
		  fPixelCenters.emplace_back(pos[0],pos[1],pos[2]);
		  ipadacc++;
		}
	      if (isector == 0) fNumChansPerSector = ipadacc;
	    }
	} // end of loop over sectors

	  // fill in the hole with a rectangular grid of pixels.   Make x,y,z=0 a pad corner.  Put pads under the cover electrode for now

      fNumChansCenter = 0;
      float rcenter = fIROCInnerRadius - 2.0;
      UInt_t nrowscenter = 2*TMath::Floor(rcenter/fCenterPadWidth);
      for (UInt_t irow = 0; irow < nrowscenter; ++irow)
	{
	  float yloc =  ( (float) irow - (float) nrowscenter/2 + 0.5 )*fCenterPadWidth;
	  UInt_t numpads = 2*TMath::Floor( TMath::Sqrt(rcenter*rcenter - yloc*yloc)/fCenterPadWidth );
	  fCenterNumPadsPerRow.push_back(numpads);
	  fCenterFirstPadInRow.push_back(fNumChansCenter);
	  float zloc = 0;
	  for (UInt_t ipad = 0; ipad < numpads; ++ipad)
	    {
	      zloc = ( (float) ipad - (float) numpads/2 + 0.5 )*fCenterPadWidth;
	      XYZPos pixpos(-fXPlaneLoc+fTPCCenter.x,yloc+fTPCCenter.y,zloc+fTPCCenter.z);
	      fPixelCenters.push_back(pixpos);
	      fNumChansCenter++;
	    }
	}

      //done with all channels on one side

      UInt_t numpixside = fPixelCenters.size();

      // duplicate the yz geometry for the opposite side -- positive x

      for (UInt_t ipix = 0; ipix < numpixside; ++ipix)
	{
	  XYZPos pixpos(fXPlaneLoc+fTPCCenter.x,fPixelCenters[ipix].y,fPixelCenters[ipix].z);
	  fPixelCenters.push_back(pixpos);
	  //std::cout << "trjpix " << fPixelCenters[ipix].z << " " << fPixelCenters[ipix].y << std::endl;
	}


      LOG_DEBUG("ChannelMapStandardAlg")
	<< "There are "
	<< fPixelCenters.size()
	<< " pixels and each sector has "
	<< fNumPadsPerRow.size()
	<< " pad rows";

      CheckPositions();

      return;
    }
    
    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::Uninitialize()
    {
      fROC = 0;
    }

    void ChannelMapStandardAlg::CheckPositions()
    {
      std::cout << "gar::ChannelMapStandardAlg::CheckPositions -- checking positions" << std::endl;

      UInt_t numchans = Nchannels();
      float xyz[3] = {0, 0, 0};

      for (UInt_t ichan = 0; ichan < numchans; ++ichan)
	{
	  ChannelToPosition(ichan, xyz);
	  UInt_t chancheck = NearestChannel(xyz);
	  if (chancheck != ichan)
	    {
	      std::cout << "gar::ChannelMapStandardAlg::CheckPositions mismatch, input chan, xyz, output chan " << ichan << " "
			<< xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << chancheck << std::endl;
	    }
	}

      std::cout << "gar::ChannelMapStandardAlg::CheckPositions -- done checking positions" << std::endl;
    }

    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::Nchannels() const
    {
      return (fPixelCenters.size());
    }

    //----------------------------------------------------------------------------
    unsigned int ChannelMapStandardAlg::NearestChannel(float const* xyz_input) const
    {

      float xyz[3] = {xyz_input[0] - fTPCCenter.x, xyz_input[1] - fTPCCenter.y, xyz_input[2] - fTPCCenter.z};

      // so far the calculations are formulaic lookups not relying on fPixelCenters.  So
      UInt_t ichan=fGapChannelNumber;  // 10 million is the flag channel for no channel at this point

      float phi = TMath::ATan2(xyz[1],xyz[2]);
      if (phi<0) phi += 2.0*TMath::Pi();
      float phisc = phi/( fPhiSectorWidth );
      UInt_t isector = TMath::Floor(phisc); // assumes the sector boundary is at phi=0  // goes from 0 to 17
	  
      // rotate this back down to a single sector -- to do -- test this!

      float rotang = TMath::DegToRad()*( isector*360/fNumSectors + fSectorOffsetAngleDeg );
      float crot = TMath::Cos(rotang);
      float srot = TMath::Sin(rotang);
      float zrot =    xyz[2]*crot + xyz[1]*srot;
      float yrot =  - xyz[2]*srot + xyz[1]*crot;

      //std::cout << "zrot, yrot, isector: " << zrot << " " << yrot << " " << isector << std::endl;
      if (zrot>fIROCInnerRadius-0.1)
	{

	  //std::cout << "Rotation: " << rotang << " " << xyz[1] << " " << xyz[2] << " " << zrot << " " << yrot << std::endl;
	  // zrot is r, and yrot=0 is centered on the midplane

	  float rtmp=zrot;
	  if (zrot > fOROCOuterRadius) return fGapChannelNumber;

	  //std::cout << "in iroc rtmp: " << rtmp << std::endl;

	  float padwidthloc = fPadWidthOROC;
	  UInt_t irow = 0;
	  if (rtmp <= fIROCOuterRadius)
	    {
	      irow = TMath::Floor( (rtmp-(fIROCInnerRadius-0.1))/fPadHeightIROC );
	      irow = TMath::Min(fNumPadRowsIROC-1,irow);
	      padwidthloc = fPadWidthIROC;
	    }
	  else if (rtmp < fOROCPadHeightChangeRadius)
	    {
	      //std::cout << "in OROCI rtmp: " << rtmp << std::endl;

	      irow =  TMath::Floor( (rtmp-(fOROCInnerRadius-0.1))/fPadHeightOROCI ) + fNumPadRowsIROC;
	      //std::cout << "Inner OROC row calc: " << irow << " " << fNumPadRowsIROC << std::endl;
	      padwidthloc = fPadWidthOROC;
	    }
	  else
	    {
	      //std::cout << "in OROCO rtmp: " << rtmp << std::endl;

	      irow = TMath::Floor( (rtmp-fOROCPadHeightChangeRadius)/fPadHeightOROCO ) + fNumPadRowsIROC + fNumPadRowsOROCI;
	      padwidthloc = fPadWidthOROC;
	    }

	  //std::cout << "irow: " << irow << std::endl;

	  UInt_t totpadrows = fNumPadRowsIROC + fNumPadRowsOROCI + fNumPadRowsOROCO;
	  if (irow >= totpadrows) return fGapChannelNumber;


	  int ichanrowoff = TMath::Floor(yrot/padwidthloc) + fNumPadsPerRow.at(irow)/2;
	  if (ichanrowoff < 0) return fGapChannelNumber;
	  if ( (UInt_t) ichanrowoff >= fNumPadsPerRow.at(irow)) return fGapChannelNumber;
	  UInt_t ichansector = fFirstPadInRow.at(irow) + ichanrowoff;
	  //std::cout << "ichansector calc: " << yrot/padwidthloc << " " << fNumPadsPerRow.at(irow) << " " << fFirstPadInRow.at(irow) << " " << ichansector << std::endl;

	  ichan = ichansector + fNumChansPerSector * isector;

	  //if (ichan>1000000)
	  //{
	  //  std::cout << "Channel ID. in gaps " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << ichan << std::endl;
	  //  std::cout << irow << " " << rtmp << " " << r <<  " " << ichansector << " " << isector << " " << fNumChansPerSector << std::endl;
	  //  std::cout << yrot << " " << zrot << " " << crot << " " << srot << std::endl;
	  //}

	} // end test if we are outside the inner radius of the ALICE chambers
      else  // must be in the hole filler
	{
	  float tvar = xyz[1]/fCenterPadWidth + fCenterNumPadsPerRow.size()/2;
	  if (tvar<0) return fGapChannelNumber;
	  UInt_t irow = TMath::Floor(tvar);
	  if (irow > fCenterFirstPadInRow.size()-1) return fGapChannelNumber;
	  tvar = xyz[2]/fCenterPadWidth + fCenterNumPadsPerRow.at(irow)/2.0;
	  if (tvar<0) return fGapChannelNumber;
	  UInt_t ivar =  TMath::Floor(tvar);
	  if (ivar >=  fCenterNumPadsPerRow.at(irow)) return fGapChannelNumber;
	  ichan = ivar + fCenterFirstPadInRow.at(irow) + fNumSectors*fNumChansPerSector;
	  // 	  ichan = fCenterFirstPadInRow.at(irow) + TMath::Min((UInt_t) fCenterNumPadsPerRow.at(irow)-1,
	  //						     (UInt_t) TMath::Floor(TMath::Max((float)0.0,(float) (xyz[2]/fCenterPadWidth + fCenterNumPadsPerRow.at(irow)/2))))
	  // +  fNumSectors*fNumChansPerSector;

	  //if (ichan>1000000)
	  //{
	  //  std::cout << "Problem Channel ID, inner filler " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << ichan << std::endl;
	  //  std::cout << irow << std::endl;
	  //}

	}

      if (xyz[0] > 0) ichan += fNumSectors*fNumChansPerSector + fNumChansCenter;  // the opposite side of the TPC.

      ichan = TMath::Max(TMath::Min(ichan, (UInt_t) fPixelCenters.size() - 1), (UInt_t) 0);
      return ichan;
    }

    //----------------------------------------------------------------------------
    void ChannelMapStandardAlg::ChannelToPosition(unsigned int chan,
						  float*       xyz)  const
    {
      xyz[0] = fPixelCenters[chan].x;
      xyz[1] = fPixelCenters[chan].y;
      xyz[2] = fPixelCenters[chan].z;

      return;
    }

  } // namespace
} // namespace gar
