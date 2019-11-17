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

// nutools includes
#include "MPDMagneticField.h"

#include "TGeoManager.h"
#include "TMath.h"
#include "TVector3.h"

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
	  << "cannot locat volume "
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

}// namespace

namespace mag {

  DEFINE_ART_SERVICE(MPDMagneticField)

} // namespace mag
