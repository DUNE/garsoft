////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapStandardAlg.cxx
/// \brief Interface to algorithm class for the standard, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov, trj@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geometry/ChannelMapAlgs/ChannelMapStandardAlg_toad.h"
#include "TMath.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace gar {
    namespace geo{
        namespace seg {

            //----------------------------------------------------------------------------
            ChannelMapStandardAlgToad::ChannelMapStandardAlgToad(fhicl::ParameterSet const& )
            {
            }

            //----------------------------------------------------------------------------

            void ChannelMapStandardAlgToad::Initialize(GeometryCore &geo)
            {
                // start over:
                Uninitialize();

		fGeo = &geo;

                fROC = AliTPCROC::Instance();

                //std::string driftvolname = geo.GetGArTPCVolumeName();

                // get these from the geometry?
                // all dimensions are in cm

                fNumSectors = 18;
                fSectorOffsetAngleDeg = 10.0;
                fPhiSectorWidth = 2.0*TMath::Pi()/18;
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
                //fXPlaneLoc = geo.TPCLength()/2.0; // in cm
                fXPlaneLoc = 62.5/2 - 28.4 + 1.0; // hardcoded, should probably write some methods to acess these
                // Note: OROC Main is at distance
                // offset_z = self.RailLength/2 - self.TerminatorOffset + self.TerminatorOROCoffset + self.orocThickness_dy1/2
                // ignore thickness, place it on the face of the OROC main
                // Q("625mm") - Q("284mm") +  Q("10mm")
                std::cout << "Channel map location in the drift direction " << fXPlaneLoc << " cm" << std::endl;
                //float TsectorH = TMath::Tan(TMath::DegToRad()*(360/(fNumSectors*2)));

                fTPCCenter.x = geo.TPCXCent();
                fTPCCenter.y = geo.TPCYCent() + fOROCInnerRadius + (fOROCOuterRadius - fOROCInnerRadius)/2 + 7.75 ;
                fTPCCenter.z = geo.TPCZCent();
                // Initialize at 0,0+ offset in y,0
                // Offset in y is fOROCInnerRadius to make top of the OROC at 0, 0, 0
                // To make it centered at 0, 0, 0, must add (fOROCOuterRadius - fOROCInnerRadius)/2
                // To make it overlap with the geometry, must add y offset from gdml file -> OROC frame position in y 77.5 millimeter
    
                std::cout << "Initializing TPC channel standard map alg with TPC Center:  " << fTPCCenter.x << " " << fTPCCenter.y << " " << fTPCCenter.z << std::endl;

                // count up channels
                // get the xyz center for each pixel and fill pad row information -- nominal geometry, negative x side


                //for (UInt_t isector = 0; isector < fNumSectors; ++isector)
                //{
                    UInt_t isector = 13;
                    UInt_t ipadacc=0;

                    // Inner Readout Chamber
                    /*for (UInt_t irow = 0; irow < fNumPadRowsIROC; ++irow)
                    {
                        UInt_t numpads = fROC->GetNPads(isector,irow);
                        //std::cout << "Row " << irow << " numpads " << numpads << std::endl;

                        if (isector == 13)
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
                            fPixelCenters.emplace_back(pos[0]+fTPCCenter.x,pos[1]+fTPCCenter.y,pos[2]+fTPCCenter.z);
                            ipadacc++;
                        }
                    }*/

                    // Outer readout chamber

                    for (UInt_t irow = 0; irow < fNumPadRowsOROCI + fNumPadRowsOROCO; ++irow)
                    {

                        UInt_t numpads = fROC->GetNPads(isector + 36,irow);
                        //std::cout << "Row " << irow << " numpads " << numpads << std::endl;
                        if (isector == 13)
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
                            fPixelCenters.emplace_back(pos[0]+fTPCCenter.x,pos[1]+fTPCCenter.y,pos[2]+fTPCCenter.z);
                            // overlay with the geoemtry here I think
                            // pos[0] is x
                            // pos[1] is y
                            // pos[3] is z
                            ipadacc++;
                        }
                        if (isector == 13) fNumChansPerSector = ipadacc;
                    }
                //} // end of loop over sectors

                // fill in the hole with a rectangular grid of pixels.   Make x,y,z=0 a pad corner.  Put pads under the cover electrode for now

                /*fNumChansCenter = 0;
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
                */

                //done with all channels on one side

                UInt_t numpixside = fPixelCenters.size();

                // duplicate the yz geometry for the opposite side -- positive x

                for (UInt_t ipix = 0; ipix < numpixside; ++ipix)
                {
                    XYZPos pixpos(fXPlaneLoc+fTPCCenter.x,fPixelCenters[ipix].y,fPixelCenters[ipix].z);
                    fPixelCenters.push_back(pixpos);
                    //std::cout << "trjpix " << fPixelCenters[ipix].z << " " << fPixelCenters[ipix].y << std::endl;
                }


                MF_LOG_DEBUG("ChannelMapStandardAlg")
                << "There are "
                << fPixelCenters.size()
                << " pixels and each sector has "
                << fNumPadsPerRow.size()
                << " pad rows";

                CheckPositions();

                return;
            }

            //----------------------------------------------------------------------------
            void ChannelMapStandardAlgToad::Uninitialize()
            {
                fROC = 0;
            }

            void ChannelMapStandardAlgToad::CheckPositions()
            {
                //std::cout << "gar::ChannelMapStandardAlgToad::CheckPositions -- checking positions" << std::endl;

                UInt_t numchans = Nchannels();
                float xyz[3] = {0, 0, 0};

                for (UInt_t ichan = 0; ichan < numchans; ++ichan)
                {

                    ChannelToPosition(ichan, xyz);
		    // if we are in just one drift volume, check only the postiive X channels
                    if (xyz[0] < fGeo->TPCXCent() && fGeo->TPCNumDriftVols() == 1) continue; 

                    UInt_t chancheck = NearestChannel(xyz);
                    if (chancheck != ichan)
                    {
                        std::cout << "gar::ChannelMapStandardAlgToad::CheckPositions mismatch, input chan, xyz, output chan " << ichan << " "
                        << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << chancheck << std::endl;
                    }

                    //if (ichan == 30000 || ichan == 100000 || ichan == 230000 || ichan == 13 || ichan == 3001 || ichan == 6001 || ichan == 309000)
                    //{
                    //std::cout << "trjc " << ichan  << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
                    //std::cout << " ====================================================" << std::endl;
                    //std::cout << "trjc " << ichan << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
                    //float varxyz[3] = {0,0,0};
                    //for (float dy=-2.0; dy<2.0; dy += 0.01)
                    //	{
                    //	  varxyz[1] = xyz[1] + dy;
                    //	  for (float dz=-2.0; dz<2.0; dz += 0.01)
                    //	    {
                    //	      varxyz[2] = xyz[2] + dz;
                    //        UInt_t chancheck2 = NearestChannel(varxyz);
                    //	      if (chancheck2 == ichan)
                    //		{
                    //		  std::cout << "trjc " << ichan << " " << varxyz[1] << " " << varxyz[2] << std::endl;
                    //		}
                    //	    }
                    //	}
                    //}

                }

                std::cout << "gar::ChannelMapStandardAlgToad::CheckPositions -- done checking positions" << std::endl;
            }

            //----------------------------------------------------------------------------
            unsigned int ChannelMapStandardAlgToad::Nchannels() const
            {
                return (fPixelCenters.size());
            }

            //----------------------------------------------------------------------------
            // wrapper for backward compatibility -- versiont that does not return the roctype.

            unsigned int ChannelMapStandardAlgToad::NearestChannel(float const* xyz) const
            {
                gar::geo::ROCType roctype = HFILLER;
                unsigned int nearestchannel;
                NearestChannelWithROCType(xyz, roctype, nearestchannel);
                return nearestchannel;
            }

            //----------------------------------------------------------------------------
            void ChannelMapStandardAlgToad::NearestChannelWithROCType(float const *xyz_input, gar::geo::ROCType &roctype, unsigned int &nearestchannel) const
            {
                roctype = HFILLER;

                float xyz[3] = {xyz_input[0] - fTPCCenter.x, xyz_input[1] - fTPCCenter.y, xyz_input[2] - fTPCCenter.z};

                // so far the calculations are formulaic lookups not relying on fPixelCenters.  So
                UInt_t ichan=fGapChannelNumber;  // 10 million is the flag channel for no channel at this point

                float phi = TMath::ATan2(xyz[1],xyz[2]);
                if (phi<0) phi += 2.0*TMath::Pi();
                float phisc = phi/( fPhiSectorWidth );
                UInt_t isector = TMath::Floor(phisc); // assumes the sector boundary is at phi=0  // goes from 0 to 17
                //std::cout << "Isector " << isector << std::endl;
                
                // rotate this back down to a single sector

                float rotang = TMath::DegToRad()*( isector*360/fNumSectors + fSectorOffsetAngleDeg );
                float crot = TMath::Cos(rotang);
                float srot = TMath::Sin(rotang);
                float zrot =    xyz[2]*crot + xyz[1]*srot;
                float yrot =  - xyz[2]*srot + xyz[1]*crot;


                //std::cout << "zrot, yrot, isector: " << zrot << " " << yrot << " " << isector << std::endl;
                if (isector != 13){
                    nearestchannel = fGapChannelNumber;
                    return;
                }
                else if (zrot>fIROCInnerRadius-0.4)
                {
                    
                    roctype = IROC;

                    //std::cout << "Rotation: " << rotang << " " << xyz[1] << " " << xyz[2] << " " << zrot << " " << yrot << std::endl;
                    // zrot is r, and yrot=0 is centered on the midplane
                
                    float rtmp=zrot;
                    if (zrot > fOROCOuterRadius)
                    {
                        nearestchannel = fGapChannelNumber;
                        return;
                    }

                    //std::cout << "in iroc rtmp: " << rtmp << std::endl;

                    float padwidthloc = fPadWidthOROC;
                    UInt_t irow = 0;
                    if (rtmp <= fIROCOuterRadius)
                    {
                        irow = TMath::Floor( (rtmp-(fIROCInnerRadius-0.4))/fPadHeightIROC );
                        if (irow >= fNumPadRowsIROC)
                        {
                            nearestchannel = fGapChannelNumber;
                            return;
                        }
                        // don't be this forgiving irow = TMath::Min(fNumPadRowsIROC-1,irow);
                        padwidthloc = fPadWidthIROC;
                        
                        nearestchannel = fGapChannelNumber;
                        return;
                    
                    }
                    else if (rtmp < fOROCPadHeightChangeRadius)
                    { //std::cout<< "isector is 13 IOROC " << std::endl;
                        //std::cout << "in IOROC rtmp: " << rtmp << std::endl;

                        roctype = IOROC;

                        //irow =  TMath::Floor( (rtmp-(fOROCInnerRadius-0.5))/fPadHeightOROCI ) + fNumPadRowsIROC;
                        irow =  TMath::Floor( (rtmp-(fOROCInnerRadius-0.5))/fPadHeightOROCI );

                        //std::cout << "in IOROC irow " << irow << std::endl;
                        //if (irow < fNumPadRowsIROC)
                        //{
                        //    nearestchannel = fGapChannelNumber;
                        //    std::cout << "IOROC return " << nearestchannel << std::endl;
                        //    return;
                        //}
                        //std::cout << "Inner OROC row calc: " << irow << " " << fNumPadRowsIROC << std::endl;
                        padwidthloc = fPadWidthOROC;
                    }
                    else
                    { //std::cout<< "isector is 13 OOROC " << std::endl;
                        //std::cout << "in OOROC rtmp: " << rtmp << std::endl;

                        roctype = OOROC;

                        //irow = TMath::Floor( (rtmp-fOROCPadHeightChangeRadius)/fPadHeightOROCO ) + fNumPadRowsIROC + fNumPadRowsOROCI;
                        irow = TMath::Floor( (rtmp-fOROCPadHeightChangeRadius)/fPadHeightOROCO ) + fNumPadRowsOROCI;
                        //std::cout << "in OOROC irow " << irow << std::endl;

                        padwidthloc = fPadWidthOROC;
                    }

                    //std::cout << "irow: " << irow << std::endl;

                    //UInt_t totpadrows = fNumPadRowsIROC + fNumPadRowsOROCI + fNumPadRowsOROCO;
                    UInt_t totpadrows = fNumPadRowsOROCI + fNumPadRowsOROCO;
                    //std::cout << "totpadrows: " << totpadrows << std::endl;

                    if (irow >= totpadrows)
                    {
                        nearestchannel = fGapChannelNumber;
                        //std::cout << " irow >= totpadrows nearestchannel: " << nearestchannel << std::endl;
                        return;
                    }

                    // to do -- put in small corrections to make the pads projective

                    int ichanrowoff = TMath::Floor(yrot/padwidthloc) + fNumPadsPerRow.at(irow)/2;
                    if (ichanrowoff < 0 || (UInt_t) ichanrowoff >= fNumPadsPerRow.at(irow))
                    {
                        nearestchannel = fGapChannelNumber;
                        //std::cout << " (ichanrowoff < 0 || (UInt_t) ichanrowoff >= fNumPadsPerRow.at(irow)) nearestchannel: " << nearestchannel << std::endl;
                        return;
                    }
                    UInt_t ichansector = fFirstPadInRow.at(irow) + ichanrowoff;
                    //std::cout << "ichansector" << ichansector << std::endl;
                    //std::cout << "ichansector calc: " << yrot/padwidthloc << " " << fNumPadsPerRow.at(irow) << " " << fFirstPadInRow.at(irow) << " " << ichansector << std::endl;

                    ichan = ichansector + fNumChansPerSector * 1;
                    //std::cout << "ichan" << ichan<< std::endl;

                    //if (ichan>1000000)
                    //{
                    //  std::cout << "Channel ID. in gaps " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << ichan << std::endl;
                    //  std::cout << irow << " " << rtmp << " " << r <<  " " << ichansector << " " << isector << " " << fNumChansPerSector << std::endl;
                    //  std::cout << yrot << " " << zrot << " " << crot << " " << srot << std::endl;
                    //}

                } // end test if we are outside the inner radius of the ALICE chambers
                else  // must be in the hole filler
                { //std::cout<< "isector is 13 Filler " << std::endl;
                    roctype = HFILLER;
                    nearestchannel = fGapChannelNumber;
                    return;
                    /*
                    float tvar = xyz[1]/fCenterPadWidth + fCenterNumPadsPerRow.size()/2;
                    if (tvar<0)
                    {
                        nearestchannel = fGapChannelNumber;
                        return;
                    }
                    UInt_t irow = TMath::Floor(tvar);
                    if (irow > fCenterFirstPadInRow.size()-1)
                    {
                        nearestchannel = fGapChannelNumber;
                        return;
                    }
                    tvar = xyz[2]/fCenterPadWidth + fCenterNumPadsPerRow.at(irow)/2.0;
                    if (tvar<0)
                    {
                        nearestchannel = fGapChannelNumber;
                        return;
                    }
                    UInt_t ivar =  TMath::Floor(tvar);
                    if (ivar >=  fCenterNumPadsPerRow.at(irow))
                    {
                        nearestchannel = fGapChannelNumber;
                        return;
                    }
                    ichan = ivar + fCenterFirstPadInRow.at(irow) + fNumSectors*fNumChansPerSector;
                    std::cout << "Ichan" << ichan << std::endl;
                    // 	  ichan = fCenterFirstPadInRow.at(irow) + TMath::Min((UInt_t) fCenterNumPadsPerRow.at(irow)-1,
                    //						     (UInt_t) TMath::Floor(TMath::Max((float)0.0,(float) (xyz[2]/fCenterPadWidth + fCenterNumPadsPerRow.at(irow)/2))))
                    // +  fNumSectors*fNumChansPerSector;

                    //if (ichan>1000000)
                    //{
                    //  std::cout << "Problem Channel ID, inner filler " << xyz[0] << " " << xyz[1] << " " << xyz[2] << " " << ichan << std::endl;
                    //  std::cout << irow << std::endl;
                    //}
                    */

                }

		// assume positive side X for the lone readout plane if we are using just one

            //if (xyz[0] > fGeo->TPCXCent() || fGeo->TPCNumDriftVols() == 1) 
		    //{
		    //ichan += fNumSectors*fNumChansPerSector + fNumChansCenter;  // the opposite side of the TPC.
		    //}

                nearestchannel = TMath::Max(TMath::Min(ichan, (UInt_t) fPixelCenters.size() - 1), (UInt_t) 0);
                //std::cout << "Nearest channel " << nearestchannel << std::endl;
                return;
            }

            //----------------------------------------------------------------------------
            // for a particular point xyz, find the nearest channel but also make a list of neighboring channels in the pad rows and in the previous and next pad
            // rows, if they don't go out of the ROC.

            void ChannelMapStandardAlgToad::NearestChannelInfo(float const* xyz, gar::geo::ChanWithNeighbors &cwn) const
            {

                // make a list of nearby pads we want to distribute charge using the pad response function
                TVector3 xvec(xyz[0],xyz[1],xyz[2]);
                float xyztest[3] = {0,0,0};
                TVector3 xvectestm(0,0,0);
                bool havem = false;
                TVector3 xvectestp(0,0,0);
                bool havep = false;

                TVector3 zerovec(0,0,0);  // dummy vector to put in until we have directions worked out

                cwn.clear();
                unsigned int chan;
                gar::geo::ROCType roctype;
                //std::cout<<" NearestChannelWithROCType(xyz,roctype,chan) " << std::endl;
                NearestChannelWithROCType(xyz,roctype,chan);   // get the nearest channel ID and add it to the list
                //std::cout<<" ChannelToPosition(chan,xyztest); " << std::endl;
                ChannelToPosition(chan,xyztest);
                TVector3 xvecchan(xyztest[0],xyztest[1],xyztest[2]);
                gar::geo::ChanWithPos centerchanwithpos = {chan, xvecchan, zerovec, roctype}; // put the direction vector in later
                cwn.push_back(centerchanwithpos);
                if (chan == fGapChannelNumber) return;         // we're in a gap so don't look for neighbors

                if (chan>0)   // look on one side in this pad row -- assume channels are numbered along pad rows
                {
                    ChannelToPosition(chan-1,xyztest);        // if we're in a gap, we get -999's here.
                    xvectestm.SetXYZ(xyztest[0],xyztest[1],xyztest[2]);
                    if ( (xvecchan-xvectestm).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                    {
                        gar::geo::ChanWithPos cwp = {chan-1, xvectestm, zerovec, roctype};  // put the direction vector in later
                        cwn.push_back(cwp);
                        havem = true;
                    }
                }
                ChannelToPosition(chan+1,xyztest);  // now look on the other side, same pad row.
                xvectestp.SetXYZ(xyztest[0],xyztest[1],xyztest[2]);
                if ( (xvecchan-xvectestp).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                {
                    gar::geo::ChanWithPos cwp = {chan+1, xvectestp, zerovec, roctype};  // put the direction vector in later
                    cwn.push_back(cwp);
                    havep = true;
                }

                TVector3 dvec(0,0,0);
                if (havem)
                {
                    dvec = xvectestm - xvecchan;
                }
                else if (havep)
                {
                    dvec = xvectestp - xvecchan;
                }
                else
                {
                    return;
                    // we have a lone pad in this row.  No geometrical info to go to another row, so skip.
                    // shouldn't happen, unless a simulated hole-filler pad row has just one pad in it.
                    // throw cet::exception("ChannelMapStandardAlg")
                    // << "Pad has no neighboring pads in its row, geometry problem " << chan;
                }
                float dvm = dvec.Mag();
                if (dvm == 0)
                {
                    throw cet::exception("ChannelMapStandardAlg")
                    << "Pad neighbor has same coordinates as pad, geometry problem " << chan;
                }
                dvec *= 1.0/dvm;
                // fill in the pad row direction vectors now that we know them.
                cwn.at(0).padrowdir = dvec;
                if (havem || havep)
                {
                    cwn.at(1).padrowdir = dvec;
                }
                if (havem && havep)
                {
                    cwn.at(2).padrowdir = dvec;
                }


                TVector3 ndp(0,-dvec.Z(),dvec.Y());
                float mag = ndp.Mag();
                ndp *= (0.8/mag); // go out 8 mm -- guarantee to get to the next pad row.

                TVector3 nextrowhyp = xvecchan + ndp;
                gar::geo::ROCType rtp = HFILLER;
                unsigned int nextrowchan;
                float nextrowhyparray[3] = {0,0,0};
                nextrowhyp.GetXYZ(nextrowhyparray);
                NearestChannelWithROCType(nextrowhyparray,rtp,nextrowchan);
                if (nextrowchan != fGapChannelNumber)
                {
                    ChannelToPosition(nextrowchan,nextrowhyparray);
                    TVector3 nrv(nextrowhyparray[0],nextrowhyparray[1],nextrowhyparray[2]);
                    gar::geo::ChanWithPos cwp = {nextrowchan, nrv, dvec, rtp};
                    cwn.push_back(cwp);
                    if (nextrowchan>0)
                    {
                        ChannelToPosition(nextrowchan-1,xyztest);
                        TVector3 xvt(xyztest[0],xyztest[1],xyztest[2]);
                        if ( (xvecchan-xvt).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                        {
                            gar::geo::ChanWithPos cwp2 = {nextrowchan-1, xvt, dvec, rtp};
                            cwn.push_back(cwp2);
                        }
                    }
                    ChannelToPosition(nextrowchan+1,xyztest);
                    TVector3 xvtp(xyztest[0],xyztest[1],xyztest[2]);
                    if ( (xvecchan-xvtp).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                    {
                        gar::geo::ChanWithPos cwp2 = {nextrowchan+1, xvtp, dvec, rtp};
                        cwn.push_back(cwp2);
                    }
                }

                // look for another row going the other way.
                nextrowhyp = xvecchan - ndp;
                nextrowhyp.GetXYZ(nextrowhyparray);
                NearestChannelWithROCType(nextrowhyparray,rtp,nextrowchan);
                if (nextrowchan != fGapChannelNumber)
                {
                    ChannelToPosition(nextrowchan,nextrowhyparray);
                    TVector3 nrv(nextrowhyparray[0],nextrowhyparray[1],nextrowhyparray[2]);
                    gar::geo::ChanWithPos cwp = {nextrowchan, nrv, dvec, rtp};
                    cwn.push_back(cwp);

                    if (nextrowchan>0)
                    {
                        ChannelToPosition(nextrowchan-1,xyztest);
                        TVector3 xvt(xyztest[0],xyztest[1],xyztest[2]);
                        if ( (xvecchan-xvt).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                        {
                            gar::geo::ChanWithPos cwp2 = {nextrowchan-1, xvt, dvec, rtp};
                            cwn.push_back(cwp2);
                        }
                    }
                    ChannelToPosition(nextrowchan+1,xyztest);
                    TVector3 xvtp(xyztest[0],xyztest[1],xyztest[2]);
                    if ( (xvecchan-xvtp).Mag() < 5.0 )   // if we run off the end of the row (start another), skip this side.
                    {
                        gar::geo::ChanWithPos cwp2 = {nextrowchan+1, xvtp, dvec, rtp};
                        cwn.push_back(cwp2);
                    }
                }
            }


            //----------------------------------------------------------------------------
            void ChannelMapStandardAlgToad::ChannelToPosition(unsigned int chan,
            float*       xyz)  const
            {
                if (chan < fPixelCenters.size())
                {
                    xyz[0] = fPixelCenters[chan].x;
                    xyz[1] = fPixelCenters[chan].y;
                    xyz[2] = fPixelCenters[chan].z;
                }
                else  // gap channel
                {
                    xyz[0] = -999;
                    xyz[1] = -999;
                    xyz[2] = -999;
                }
                return;
            }

        }
    } // namespace
} // namespace gar
