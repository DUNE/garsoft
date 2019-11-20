//////////////////////////////////////////////////////////////////////////
/// \file MPDMagneticField_service.cc
///
/// \author trj@fnal.gov, from the nutools example by David McKee
//////////////////////////////////////////////////////////////////////////
/// \class MPDMagneticField MPDMagneticField.h 
/// The initial implementation will be trivial: simply supporting a
/// constant field in a named detector volume. In principle we should
/// read a full field map from an external file of some kind.
///
/// FHICL values so far
///
///    - "UseField" a integer. When 0 we don't even instantiate a
///      Magnetic field object. Describes the description to use.
///    - "Constant Field" a vector< double > which should have three
///      elements and is interpreted in Tesla
///    - "MagnetizedVolume" names the G4logical volume to which the
///      field should be attached
//////////////////////////////////////////////////////////////////////////

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

// nutools includes
#include "MPDMagneticField.h"

#include "TGeoManager.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include <vector>
#include <string>
#include <sys/stat.h>

#include "cetlib/search_path.h"

namespace mag {

  MPDMagneticField::MPDMagneticField(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------
  void MPDMagneticField::reconfigure(fhicl::ParameterSet const& pset)
  {
    auto fieldDescriptions = pset.get<std::vector<fhicl::ParameterSet> >("FieldDescriptions");
    
    MPDMagneticFieldDescription fieldDescription;
    for(auto itr : fieldDescriptions){
      fieldDescription.fMode   = (mag::MagFieldMode_t)(itr.get<int>("UseField"));

      // if the mode is turned to off, no point looking for a volume or 
      // trying to put a description into the fFieldDescriptions data member.
      // if all input field descriptions are set to fMode = kNoBFieldMode, then the 
      // methods to return the fields will not go into the loop over fFieldDescriptions
      // and will just return a 0 field.
      if(fieldDescription.fMode == mag::kNoBFieldMode) continue;
      fieldDescription.fScaleFactor = itr.get<float>("ScaleFactor",1.0);

      fieldDescription.fVolume = itr.get<std::string>("MagnetizedVolume");
      fieldDescription.fGeoVol = gGeoManager->FindVolumeFast(fieldDescription.fVolume.c_str());

      // check that we have a good volume
      if( fieldDescription.fGeoVol == nullptr )
        throw cet::exception("MagneticField")
	  << "cannot locate volume "
	  << fieldDescription.fVolume
	  << " in gGeoManager, bail";

      // These need to be read as types that FHICL know about, but they
      // are used by Geant, so I store them in Geant4 types.
      std::vector<double> defaultField = {0,0,0};
      std::vector<double> field = itr.get<std::vector<double> >("ConstantField",defaultField);
      
      // Force the dimension of the field definition
      field.resize(3);
      for(size_t i = 0; i < 3; ++i) fieldDescription.fField[i] = field[i];
      
      if (fieldDescription.fMode == mag::kFieldRZMapMode)
	{
	  fieldDescription.fRZFieldMapFilename = itr.get<std::string>("RZFieldMapFilename");
	  cet::search_path sp("FW_SEARCH_PATH");
          std::string fn2 = "MPD/FieldMap/";
          fn2 += fieldDescription.fRZFieldMapFilename;
          std::string fullname;
          sp.find_file(fn2, fullname);
          struct stat sb;
          if (fullname.empty() || stat(fullname.c_str(), &sb)!=0)
	    throw cet::exception("RadioGen") << "Input spectrum file "
					     << fn2
					     << " not found in FW_SEARCH_PATH!\n";

	  std::cout << "MPDMagnetciField: Opening magnetic field RZ map in: " << fullname << std::endl;
          std::ifstream inFile(fullname, std::ios::in);
          std::string line;

	  int rcounter=-1;
	  struct RZFieldMap rzmap;
	  float rcur=0;
	  float zcur = 0;
	  rzmap.dr = 0;
	  rzmap.dz = 0;

          while (std::getline(inFile,line)) 
	    {
	      float x=0;
	      float y=0;
	      float z=0;
	      float bx=0;
	      float by=0;
	      float bz=0;
	      float b=0;
              std::stringstream linestream(line);
	      linestream >> x >> y >> z >> bx >> by >> bz >> b;
	      // we get this in meters.  convert to cm
	      x *= 100.0;
	      y *= 100.0;
	      z *= 100.0;
	      float r = TMath::Sqrt( x*x + y*y );
	      float br = bx;   // Vladimir's map defines it this way.

	      std::vector<float> emptyvec;
	      if (TMath::Abs(r-rcur)>0.001 || rcounter<0)
		{
		  rzmap.br.push_back(emptyvec);
		  rzmap.bz.push_back(emptyvec);
		  rcounter++;
		  if (r > rcur + 0.001)
		    {
		      rzmap.dr = r-rcur;
		    }
		  rcur = r;
		}
	      rzmap.br[rcounter].push_back(br);
	      rzmap.bz[rcounter].push_back(bz);
	      if (z-zcur > 0.001)
		{
		  rzmap.dz = z-zcur;
		}
	      zcur = z;
	    }

          std::vector<double> ZAxis = itr.get<std::vector<double> >("ZAxis");
	  rzmap.ZAxis.SetXYZ(ZAxis[0],ZAxis[1],ZAxis[2]);
          std::vector<double> CoordOffset = itr.get<std::vector<double> >("CoordOffset");
	  rzmap.CoordOffset.SetXYZ(CoordOffset[0],CoordOffset[1],CoordOffset[2]);

	  std::cout << "MPDMagneticField: Read in RZ map, now setting it: " << rzmap.dr << " " << rzmap.dz << std::endl;
	  std::cout << "array sizes: " << rzmap.br.size() << " " << rzmap.br[0].size() << std::endl;
          fieldDescription.fRZFieldMap = rzmap;
	}
      fFieldDescriptions.push_back(fieldDescription);
    }

    fGlobalScaleFactor = pset.get<float>("GlobalScaleFactor",1.0);

    fNBinsXCheckPlots =  pset.get<std::vector<int> >("NBinsXCheckPlots");
    fNBinsYCheckPlots =  pset.get<std::vector<int> >("NBinsYCheckPlots");
    fNBinsZCheckPlots =  pset.get<std::vector<int> >("NBinsZCheckPlots");
    
    fXLowCheckPlots = pset.get<std::vector<float> > ("XLowCheckPlots");
    fXHighCheckPlots = pset.get<std::vector<float> >("XHighCheckPlots");
    fYLowCheckPlots = pset.get<std::vector<float> > ("YLowCheckPlots");
    fYHighCheckPlots = pset.get<std::vector<float> >("YHighCheckPlots");
    fZLowCheckPlots = pset.get<std::vector<float> > ("ZLowCheckPlots");
    fZHighCheckPlots = pset.get<std::vector<float> >("ZHighCheckPlots");

    fXCentCheckPlots = pset.get<std::vector<float> > ("XCentCheckPlots");
    fYCentCheckPlots = pset.get<std::vector<float> > ("YCentCheckPlots");
    fZCentCheckPlots = pset.get<std::vector<float> > ("ZCentCheckPlots");

    MakeCheckPlots();

    return;
  }

  //------------------------------------------------------------
  G4ThreeVector const MPDMagneticField::FieldAtPoint(G4ThreeVector const& p) const
  {
    // check that the input point is in the magnetized volume
    // Use the gGeoManager to determine what node the point
    // is in
    double point[3] = { p.x(), p.y(), p.z() };
    TVector3 pv(p.x(), p.y(), p.z());
    G4ThreeVector zerofield(0,0,0);
    G4ThreeVector calculatedfield = zerofield;

    // loop over the field descriptions to see if the point is in any of them
    for(auto fd : fFieldDescriptions){
      // we found a node, see if its name is the same as
      // the volume with the field

      if(fd.fGeoVol->Contains(point)) 
	{
	  // "Automatic" for now just means constant field.
	  if (fd.fMode == mag::kConstantBFieldMode || fd.fMode == mag::kAutomaticBFieldMode )
	    { 
	      calculatedfield = fd.fField;
	    }
	  else if (fd.fMode == mag::kNoBFieldMode)
	    {
	      calculatedfield = zerofield;
	    }
	  else if (fd.fMode == mag::kFieldRZMapMode)
	    {
	      // a handy name
	      auto &rzm = fd.fRZFieldMap;
	      // check for validity

	      if (rzm.dr <= 0) 
		{
		  throw cet::exception("MPDMagneticFieldService: RZ map bad R spacing:") << rzm.dr;
		}
	      if (rzm.dz <= 0) 
		{
		  throw cet::exception("MPDMagneticFieldService: RZ map bad Z spacing:") << rzm.dz;
		}

	      TVector3 ploc = pv - rzm.CoordOffset;
	      float zloc = ploc.Dot(rzm.ZAxis);
	      float rloc = (ploc.Cross(rzm.ZAxis)).Mag();
	      float azloc = TMath::Abs(zloc);
	      if (azloc > rzm.dz * rzm.nz() ||
		  rloc  > rzm.dr * rzm.nr())
		{
		  calculatedfield = zerofield;   // outside the field map region
		}
	      else
		{
		  TVector3 fv(0,0,0);
		  size_t ibinr = rloc / rzm.dr;
		  size_t ibinz = azloc / rzm.dz;
		  float bzcomponent = rzm.bz.at(ibinr).at(ibinz);  // closest grid point

		  // cheesy interpolation in just one direction at a time.  A better interpolation
		  // would be to stretch a surface around four points.
		  // also this interpolation naively only interpolates Bz as a function of Z and
		  // br as a function of r.  the four-point surface would be more symmetrical

		  if (ibinz == 0)
		    {
		      float dbz = ((rzm.bz.at(ibinr).at(ibinz+1) - bzcomponent)/rzm.dz)*
			(azloc - (int) azloc);
		      bzcomponent += dbz;
		    }
		  else
		    {
		      float dbz = (bzcomponent - (rzm.bz.at(ibinr).at(ibinz-1))/rzm.dz)*
			(azloc - rzm.dz*((int) azloc/rzm.dz));
		      bzcomponent += dbz;
		    }
		  fv = bzcomponent * rzm.ZAxis;  // z component of field.  

		  float bradial = rzm.br.at(ibinr).at(ibinz);
		  if (ibinr == 0)
		    {
		      float dbr = ((rzm.br.at(ibinr+1).at(ibinz) - bradial)/rzm.dr)*
			(rloc - rzm.dr*((int) rloc/rzm.dr));
		      bradial += dbr;
		    }
		  else
		    {
		      float dbr = (bradial - (rzm.br.at(ibinr-1).at(ibinz))/rzm.dr)*
			(rloc - rzm.dr*((int) rloc/rzm.dr));
		      bradial += dbr;
		    }

		  // calculate radial direction of the field

		  TVector3 rperp = ploc - zloc*rzm.ZAxis;
		  float mrp = rperp.Mag();

		  // we only can have a radial component of the field if we are at r>0
		  if (mrp != 0)
		    {
		      rperp *= (1.0/mrp);
		      if (zloc > 0)           // swap the radial field direction based on the sign of z.
			{
			  fv += rperp * bradial;
			}
		      else
			{
			  fv -= rperp * bradial;
			}
		    }
		  G4ThreeVector fvg(fv.X(), fv.Y(), fv.Z());
		  calculatedfield = fvg;
		}
	    }
	  else
	    {
	      throw cet::exception("MPDMagneticField_service: Ununderstood field mode: ") << fd.fMode;
	    }
	}
      calculatedfield *= fd.fScaleFactor;
    }

    calculatedfield *= fGlobalScaleFactor;

    //std::cout << "field at point: " << p << " field: " << calculatedfield << std::endl;
    // if we get here, we can't find a field
    return calculatedfield;
  }

  //------------------------------------------------------------
  G4ThreeVector const MPDMagneticField::UniformFieldInVolume(std::string const& volName) const
  {
    // if the input volume name is the same as the magnetized volume
    // return the uniform field
    
    for(auto fd : fFieldDescriptions){
      if (fd.fVolume.compare(volName) == 0) return fd.fField;
    }
    
    // if we get here, we can't find a field
    return G4ThreeVector(0);
  }

  // make 2D plots of bx, by, and bz

  void MPDMagneticField::MakeCheckPlots()
  {
    art::ServiceHandle< ::art::TFileService> tfs;
    std::vector<int> emtpyintvec;
    std::vector<float> emtpyflotvec;
    std::vector<TString> axisname = {"X", "Y", "Z"};
    fCheckPlots.clear();

    for (size_t i=0; i< fNBinsXCheckPlots.size(); ++i)
      {
	int naxes = 0;
	std::vector<float> axislow;
	std::vector<float> axishigh;
	std::vector<int >  axisnumber;
	std::vector<int >  axisnbins;

	if (fNBinsXCheckPlots.at(i) > 0)
	  {
	    naxes ++;
	    axisnumber.push_back(0);
	    axislow.push_back(fXLowCheckPlots.at(i));
	    axishigh.push_back(fXHighCheckPlots.at(i));
	    axisnbins.push_back(fNBinsXCheckPlots.at(i));
	  }
	if (fNBinsYCheckPlots.at(i) > 0)
	  {
	    naxes ++;
	    axisnumber.push_back(1);
	    axislow.push_back(fYLowCheckPlots.at(i));
	    axishigh.push_back(fYHighCheckPlots.at(i));
	    axisnbins.push_back(fNBinsYCheckPlots.at(i));
	  }
	if (fNBinsZCheckPlots.at(i) > 0)
	  {
	    naxes ++;
	    axisnumber.push_back(2);
	    axislow.push_back(fZLowCheckPlots.at(i));
	    axishigh.push_back(fZHighCheckPlots.at(i));
	    axisnbins.push_back(fNBinsZCheckPlots.at(i));
	  }
	TString mfcheckidx = "mfcheckx";
	mfcheckidx += i;
	TString mfcheckidy = "mfchecky";
	mfcheckidy += i;
	TString mfcheckidz = "mfcheckz";
	mfcheckidz += i;
	TString htitle = "Magnetic Field Check; ";
	if (naxes == 1)
	  {
	    htitle += axisname.at(axisnumber.at(0));
	    htitle += " (cm);";
	    fCheckPlots.push_back( (TH1*) tfs->make<TH1F>(mfcheckidx,"Bx "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0)));
	    TH1 *bxhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH1F>(mfcheckidy,"By "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0)));
	    TH1 *byhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH1F>(mfcheckidz,"Bz "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0)));
	    TH1 *bzhist = fCheckPlots.back();
	    for (int j=0; j< axisnbins.at(0); ++j)
	      {
		float axispos = axislow.at(0) + (((float) j)+0.5 )*(axishigh.at(0)-axislow.at(0))/axisnbins.at(0);	         
		G4ThreeVector point( fXCentCheckPlots.at(i), fYCentCheckPlots.at(i), fZCentCheckPlots.at(i) );
		point[axisnumber.at(0)] += axispos;
		G4ThreeVector bfield = FieldAtPoint(point);
		bxhist->SetBinContent(j+1,bfield.x());
		byhist->SetBinContent(j+1,bfield.y());
		bzhist->SetBinContent(j+1,bfield.z());
	      }
	  }
	else if (naxes == 2)
	  {
	    htitle += axisname.at(axisnumber.at(0));
	    htitle += " (cm);";
	    htitle += axisname.at(axisnumber.at(1));
	    htitle += " (cm);";
	    fCheckPlots.push_back( (TH1*) tfs->make<TH2F>(mfcheckidx,"Bx "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1)));
	    TH1 *bxhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH2F>(mfcheckidy,"By "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1)));
	    TH1 *byhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH2F>(mfcheckidz,"Bz "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1)));
	    TH1 *bzhist = fCheckPlots.back();
	    for (int j=0; j<axisnbins.at(0); ++j)
	      {
		float axispos0 = axislow.at(0) + (((float) j)+0.5 )*(axishigh.at(0)-axislow.at(0))/axisnbins.at(0); 
		for (int k=0; k<axisnbins.at(1); ++k)
		  {
		    float axispos1 = axislow.at(1) + (((float) k)+0.5 )*(axishigh.at(1)-axislow.at(1))/axisnbins.at(1);
		    G4ThreeVector point( fXCentCheckPlots.at(i), fYCentCheckPlots.at(i), fZCentCheckPlots.at(i) );
		    point[axisnumber.at(0)] += axispos0;
		    point[axisnumber.at(1)] += axispos1;
		    G4ThreeVector bfield = FieldAtPoint(point);
		    bxhist->SetBinContent(j+1,k+1,bfield.x());
		    byhist->SetBinContent(j+1,k+1,bfield.y());
		    bzhist->SetBinContent(j+1,k+1,bfield.z());
		  }
	      }
	  }
	else if (naxes == 3)
	  {
	    htitle += axisname.at(axisnumber.at(0));
	    htitle += " (cm);";
	    htitle += axisname.at(axisnumber.at(1));
	    htitle += " (cm);";
	    htitle += axisname.at(axisnumber.at(2));
	    htitle += " (cm)";
	    fCheckPlots.push_back( (TH1*) tfs->make<TH3F>(mfcheckidx,"Bx "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1),
							  axisnbins.at(2),axislow.at(2),axishigh.at(2)));
	    TH1 *bxhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH3F>(mfcheckidy,"By "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1),
							  axisnbins.at(2),axislow.at(2),axishigh.at(2)));
	    TH1 *byhist = fCheckPlots.back();
	    fCheckPlots.push_back( (TH1*) tfs->make<TH3F>(mfcheckidz,"Bz "+htitle,axisnbins.at(0),axislow.at(0),axishigh.at(0),
							  axisnbins.at(1),axislow.at(1),axishigh.at(1),
							  axisnbins.at(2),axislow.at(2),axishigh.at(2)));
	    TH1 *bzhist = fCheckPlots.back();
	    for (int j=0; j<axisnbins.at(0); ++j)
	      {
		float axispos0 = axislow.at(0) + (((float) j)+0.5 )*(axishigh.at(0)-axislow.at(0))/axisnbins.at(0); 
		for (int k=0; k<axisnbins.at(1); ++k)
		  {
		    float axispos1 = axislow.at(1) + (((float) k)+0.5)*(axishigh.at(1)-axislow.at(1))/axisnbins.at(1);
		    for (int m=0; m<axisnbins.at(2); ++m)
		      {
			float axispos2 = axislow.at(2) + (((float) m)+0.5)*(axishigh.at(2)-axislow.at(2))/axisnbins.at(2);
			G4ThreeVector point( fXCentCheckPlots.at(i), fYCentCheckPlots.at(i), fZCentCheckPlots.at(i) );
			point[axisnumber.at(0)] += axispos0;
			point[axisnumber.at(1)] += axispos1;
			point[axisnumber.at(2)] += axispos2;
			G4ThreeVector bfield = FieldAtPoint(point);
			bxhist->SetBinContent(j+1,k+1,m+1,bfield.x());
			byhist->SetBinContent(j+1,k+1,m+1,bfield.y());
			bzhist->SetBinContent(j+1,k+1,m+1,bfield.z());
		      }
		  }
	      }
	  }
      }
  }

}// namespace

namespace mag {

  DEFINE_ART_SERVICE(MPDMagneticField)

} // namespace mag
