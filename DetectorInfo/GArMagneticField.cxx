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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "DetectorInfo/GArMagneticField.h"

#include "TGeoManager.h"
#include "TMath.h"
#include "TVector3.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <vector>

#include "cetlib/search_path.h"

namespace mag {

  GArMagneticField::GArMagneticField(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);
  }

  //------------------------------------------------------------
  void GArMagneticField::reconfigure(fhicl::ParameterSet const& pset)
  {
    fGlobalScaleFactor = pset.get<float>("GlobalScaleFactor", 1.0);
    fUnitFactor = pset.get<float>("UnitFactor", 1.);
    auto fieldDescriptions = pset.get<std::vector<fhicl::ParameterSet>>("FieldDescriptions");

    MagneticFieldDescription fieldDescription;
    for (auto itr : fieldDescriptions) {
      fieldDescription.fMode = (mag::MagFieldMode_t)(itr.get<int>("UseField"));

      // if the mode is turned to off, no point looking for a volume or
      // trying to put a description into the fFieldDescriptions data member.
      // if all input field descriptions are set to fMode = kNoBFieldMode, then the
      // methods to return the fields will not go into the loop over fFieldDescriptions
      // and will just return a 0 field.
      if (fieldDescription.fMode == mag::kNoBFieldMode) continue;

      fieldDescription.fScaleFactor = itr.get<float>("ScaleFactor", 1.0);
      fieldDescription.fVolume = itr.get<std::string>("MagnetizedVolume");
      fieldDescription.fGeoVol = gGeoManager->FindVolumeFast(fieldDescription.fVolume.c_str());

      // check that we have a good volume
      if (fieldDescription.fGeoVol == nullptr) {
        throw cet::exception("GArMagneticField")
          << "cannot locate volume " << fieldDescription.fVolume << " in gGeoManager, bail";
      }

      // These need to be read as types that FHICL know about, but they
      // are used by Geant, so I store them in Geant4 types.
      std::vector<double> defaultField = {0, 0, 0};
      std::vector<double> field = itr.get<std::vector<double>>("ConstantField", defaultField);

      // Force the dimension of the field definition
      field.resize(3);
      for (size_t i = 0; i < 3; ++i)
        fieldDescription.fField[i] = field[i];

      if (fieldDescription.fMode == mag::kFieldRZMapMode) {
        fieldDescription.fFieldMapFilename = itr.get<std::string>("FieldMapFilename");
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fn2 = fieldDescription.fFieldMapFilename;
        std::string fullname;
        sp.find_file(fn2, fullname);

        struct stat sb;
        if (fullname.empty() || stat(fullname.c_str(), &sb) != 0) {
          throw cet::exception("RadioGen")
            << "Input magnetic field file " << fn2 << " not found in FW_SEARCH_PATH!\n";
        }

        auto start = std::chrono::steady_clock::now();
        struct RZFieldMap rzmap;
        ReadRZFile(fullname, rzmap);
        auto end = std::chrono::steady_clock::now();

        std::vector<double> ZAxis = itr.get<std::vector<double>>("ZAxis");
        rzmap.ZAxis.SetXYZ(ZAxis[0], ZAxis[1], ZAxis[2]);

        std::vector<double> CoordOffset = itr.get<std::vector<double>>("CoordOffset");
        rzmap.CoordOffset.SetXYZ(CoordOffset[0], CoordOffset[1], CoordOffset[2]);

        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "GArMagneticField: Finished reading RZ map in " << elapsed_seconds.count()
                  << "s, now setting it: " << rzmap.dr << " " << rzmap.dz << std::endl;
        std::cout << "Array sizes (R, Z): " << rzmap.br.size() << " " << rzmap.br[0].size()
                  << std::endl;
        fieldDescription.fRZFieldMap = rzmap;
      }

      if (fieldDescription.fMode == mag::kFieldXYZMapMode) {
        fieldDescription.fFieldMapFilename = itr.get<std::string>("FieldMapFilename");
        cet::search_path sp("FW_SEARCH_PATH");
        std::string fn2 = fieldDescription.fFieldMapFilename;
        std::string fullname;
        sp.find_file(fn2, fullname);

        struct stat sb;
        if (fullname.empty() || stat(fullname.c_str(), &sb) != 0) {
          throw cet::exception("RadioGen")
            << "Input magnetic field file " << fn2 << " not found in FW_SEARCH_PATH!\n";
        }

        auto start = std::chrono::steady_clock::now();
        struct XYZFieldMap xyzmap;
        ReadXYZFile(fullname, xyzmap, fUnitFactor);
        auto end = std::chrono::steady_clock::now();

        std::vector<double> ZAxis =
          itr.get<std::vector<double>>("ZAxis", std::vector<double>{0.0, 0.0, 1.0});
        xyzmap.ZAxis.SetXYZ(ZAxis[0], ZAxis[1], ZAxis[2]);

        xyzmap.CoordOffset.SetXYZ(xyzmap.xo, xyzmap.yo, xyzmap.zo);
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "GArMagneticField: Finished reading XYZ map in " << elapsed_seconds.count()
                  << "s." << std::endl;

        xyzmap.UseSymmetry = itr.get<bool>("UseSymmetry", false);
        std::cout << "GArMagneticField: Is map symmetric - " << std::boolalpha << xyzmap.UseSymmetry
                  << std::endl;
        fieldDescription.fXYZFieldMap = xyzmap;
      }

      fFieldDescriptions.push_back(fieldDescription);
    }

    return;
  }

  //------------------------------------------------------------
  G4ThreeVector const GArMagneticField::FieldAtPoint(G4ThreeVector const& p) const
  {
    // check that the input point is in the magnetized volume
    // Use the gGeoManager to determine what node the point
    // is in
    double point[3] = {p.x(), p.y(), p.z()};
    G4ThreeVector zerofield(0, 0, 0);
    G4ThreeVector calculatedfield = zerofield;

    // loop over the field descriptions to see if the point is in any of them
    for (const auto& fd : fFieldDescriptions) {
      // we found a node, see if its name is the same as
      // the volume with the field

      if (fd.fGeoVol->Contains(point)) {
        // "Automatic" for now just means constant field.
        if (fd.fMode == mag::kConstantBFieldMode || fd.fMode == mag::kAutomaticBFieldMode) {
          calculatedfield = fd.fField;
        }
        else if (fd.fMode == mag::kNoBFieldMode) {
          calculatedfield = zerofield;
        }
        else if (fd.fMode == mag::kFieldRZMapMode) {
          calculatedfield = CalcRZField(p, fd.fRZFieldMap);
        }
        else if (fd.fMode == mag::kFieldXYZMapMode) {
          calculatedfield = CalcXYZField(p, fd.fXYZFieldMap);
        }
        else {
          throw cet::exception("GArMagneticField_service: Ununderstood field mode: ") << fd.fMode;
        }
      }
      calculatedfield *= fd.fScaleFactor;
    }

    calculatedfield *= fGlobalScaleFactor;

    // std::cout << "field at point: " << p << " field: " << calculatedfield << std::endl;
    return calculatedfield;
  }

  //------------------------------------------------------------
  G4ThreeVector const GArMagneticField::UniformFieldInVolume(std::string const& volName) const
  {
    // if the input volume name is the same as the magnetized volume
    // return the uniform field

    for (auto fd : fFieldDescriptions) {
      if (fd.fVolume.compare(volName) == 0) return fd.fField;
    }

    // if we get here, we can't find a field
    return G4ThreeVector(0);
  }

  //------------------------------------------------------------
  void GArMagneticField::ReadRZFile(const std::string& filename, RZFieldMap& rzmap)
  {
    std::cout << "GArMagneticField: Opening magnetic field RZ map in: " << filename << std::endl;
    std::ifstream inFile(filename, std::ios::in);
    std::string line;

    int rcounter = -1;
    float rcur = 0;
    float zcur = 0;
    rzmap.dr = 0;
    rzmap.dz = 0;

    while (std::getline(inFile, line)) {
      float x = 0;
      float y = 0;
      float z = 0;
      float bx = 0;
      float by = 0;
      float bz = 0;
      float b = 0;
      std::stringstream linestream(line);
      linestream >> x >> y >> z >> bx >> by >> bz >> b;
      // we get this in meters.  convert to cm
      x *= 100.0;
      y *= 100.0;
      z *= 100.0;
      float r = TMath::Sqrt(x * x + y * y);
      float br = bx; // Vladimir's map defines it this way.

      std::vector<float> emptyvec;
      if (TMath::Abs(r - rcur) > 0.001 || rcounter < 0) {
        rzmap.br.push_back(emptyvec);
        rzmap.bz.push_back(emptyvec);
        rcounter++;
        if (r > rcur + 0.001) { rzmap.dr = r - rcur; }
        rcur = r;
      }
      rzmap.br[rcounter].push_back(br);
      rzmap.bz[rcounter].push_back(bz);
      if (z - zcur > 0.001) { rzmap.dz = z - zcur; }
      zcur = z;
    }
  }

  //------------------------------------------------------------
  void GArMagneticField::ReadXYZFile(const std::string& filename,
                                     XYZFieldMap& xyzmap,
                                     const float unitFactor)
  {
    std::cout << "GArMagneticField: Opening magnetic field XYZ map in: " << filename << std::endl;
    std::fstream fin(filename, std::fstream::in);
    std::string line;

    int xcount{-1}, ycount{-1}, zcount{-1};
    float xcurr{0}, ycurr{0}, zcurr{0};

    while (std::getline(fin >> std::ws, line)) {
      std::stringstream ss(line);

      if (ss.str().front() == '#') continue;

      ss >> xyzmap.xo >> xyzmap.yo >> xyzmap.zo >> xyzmap.dx >> xyzmap.dy >> xyzmap.dz;

      // Position coordinates need to be in mm
      xyzmap.xo *= unitFactor;
      xyzmap.yo *= unitFactor;
      xyzmap.zo *= unitFactor;
      xyzmap.dx *= unitFactor;
      xyzmap.dy *= unitFactor;
      xyzmap.dz *= unitFactor;

      break;
    }

    while (std::getline(fin >> std::ws, line)) {
      float x{0}, y{0}, z{0}, fx{0}, fy{0}, fz{0}, f{0};
      std::stringstream ss(line);

      if (ss.str().front() == '#') continue;

      ss >> x >> y >> z >> fx >> fy >> fz >> f;

      // Position coordinates need to be in mm
      x *= unitFactor;
      y *= unitFactor;
      z *= unitFactor;

      if (std::abs(x - xcurr) > 0.0 || xcount < 0) {
        xcurr = x;
        xcount += 1;
        ycount = -1;
        xyzmap.fx.emplace_back(std::vector<std::vector<float>>{});
        xyzmap.fy.emplace_back(std::vector<std::vector<float>>{});
        xyzmap.fz.emplace_back(std::vector<std::vector<float>>{});
      }

      if (std::abs(y - ycurr) > 0.0 || ycount < 0) {
        ycurr = y;
        ycount += 1;
        zcount = -1;
        xyzmap.fx[xcount].emplace_back(std::vector<float>{});
        xyzmap.fy[xcount].emplace_back(std::vector<float>{});
        xyzmap.fz[xcount].emplace_back(std::vector<float>{});
      }

      xyzmap.fx[xcount][ycount].push_back(fx);
      xyzmap.fy[xcount][ycount].push_back(fy);
      xyzmap.fz[xcount][ycount].push_back(fz);

      if (std::abs(z - zcurr) > 0.0 || zcount < 0) {
        zcurr = z;
        zcount += 1;
      }
    }
  }

  //------------------------------------------------------------
  G4ThreeVector GArMagneticField::CalcRZField(G4ThreeVector const& p, RZFieldMap const& rzm) const
  {
    // check for validity
    if (rzm.dr <= 0) { throw cet::exception("GArMagneticField: RZ map bad R spacing:") << rzm.dr; }
    if (rzm.dz <= 0) { throw cet::exception("GArMagneticField: RZ map bad Z spacing:") << rzm.dz; }

    TVector3 pv(p.x(), p.y(), p.z());
    TVector3 ploc = pv - rzm.CoordOffset;
    float zloc = ploc.Dot(rzm.ZAxis);
    float rloc = (ploc.Cross(rzm.ZAxis)).Mag();
    float azloc = TMath::Abs(zloc);
    if (azloc > rzm.dz * rzm.nz() || rloc > rzm.dr * rzm.nr()) { return G4ThreeVector(0, 0, 0); }
    else {
      TVector3 fv(0, 0, 0);
      size_t ibinr = rloc / rzm.dr;
      size_t ibinz = azloc / rzm.dz;
      float bzcomponent = rzm.bz.at(ibinr).at(ibinz); // closest grid point

      // cheesy interpolation in just one direction at a time.  A better interpolation
      // would be to stretch a surface around four points.
      // also this interpolation naively only interpolates Bz as a function of Z and
      // br as a function of r.  the four-point surface would be more symmetrical

      if (ibinz == 0) {
        float dbz =
          ((rzm.bz.at(ibinr).at(ibinz + 1) - bzcomponent) / rzm.dz) * (azloc - (int)azloc);
        bzcomponent += dbz;
      }
      else {
        float dbz = (bzcomponent - (rzm.bz.at(ibinr).at(ibinz - 1)) / rzm.dz) *
                    (azloc - rzm.dz * ((int)azloc / rzm.dz));
        bzcomponent += dbz;
      }
      fv = bzcomponent * rzm.ZAxis; // z component of field.

      float bradial = rzm.br.at(ibinr).at(ibinz);
      if (ibinr == 0) {
        float dbr = ((rzm.br.at(ibinr + 1).at(ibinz) - bradial) / rzm.dr) *
                    (rloc - rzm.dr * ((int)rloc / rzm.dr));
        bradial += dbr;
      }
      else {
        float dbr = (bradial - (rzm.br.at(ibinr - 1).at(ibinz)) / rzm.dr) *
                    (rloc - rzm.dr * ((int)rloc / rzm.dr));
        bradial += dbr;
      }

      // calculate radial direction of the field

      TVector3 rperp = ploc - zloc * rzm.ZAxis;
      float mrp = rperp.Mag();

      // we only can have a radial component of the field if we are at r>0
      if (mrp != 0) {
        rperp *= (1.0 / mrp);
        if (zloc > 0) // swap the radial field direction based on the sign of z.
        {
          fv += rperp * bradial;
        }
        else {
          fv -= rperp * bradial;
        }
      }
      G4ThreeVector fvg(fv.X(), fv.Y(), fv.Z());
      return fvg;
    }
  }

  //------------------------------------------------------------
  G4ThreeVector GArMagneticField::CalcXYZField(G4ThreeVector const& p, XYZFieldMap const& fd) const
  {
    const float x = fd.UseSymmetry ? std::abs(p.x() - fd.xo) : p.x() - fd.xo;
    const float y = fd.UseSymmetry ? std::abs(p.y() - fd.yo) : p.y() - fd.yo;
    const float z = fd.UseSymmetry ? std::abs(p.z() - fd.zo) : p.z() - fd.zo;

    Interpolator interp;
    const float bx = interp.interpolate(x, y, z, fd.fx, fd.dx, fd.dy, fd.dz, 0.0, 0.0, 0.0);
    const float by = interp.interpolate(x, y, z, fd.fy, fd.dx, fd.dy, fd.dz, 0.0, 0.0, 0.0);
    const float bz = interp.interpolate(x, y, z, fd.fz, fd.dx, fd.dy, fd.dz, 0.0, 0.0, 0.0);

    return G4ThreeVector(bx, by, bz);
  }

  //------------------------------------------------------------
  Interpolator::Interpolator() {}

  //------------------------------------------------------------
  float Interpolator::interpolate(const float* point,
                                  const std::vector<std::vector<std::vector<float>>>& g,
                                  const float* delta,
                                  const float* offset) const
  {
    return interpolate(point[0],
                       point[1],
                       point[2],
                       g,
                       delta[0],
                       delta[1],
                       delta[2],
                       offset[0],
                       offset[1],
                       offset[2]);
  }

  //------------------------------------------------------------
  float Interpolator::interpolate(float x,
                                  float y,
                                  float z,
                                  const std::vector<std::vector<std::vector<float>>>& g,
                                  float hx,
                                  float hy,
                                  float hz,
                                  float xo,
                                  float yo,
                                  float zo) const
  {
    float v = 0;

    const int xsize = g.size();
    const int ysize = g[0].size();
    const int zsize = g[0][0].size();

    const float xp = (x - xo) / hx;
    const float yp = (y - yo) / hy;
    const float zp = (z - zo) / hz;

    const int x_idx = std::floor(xp);
    const int y_idx = std::floor(yp);
    const int z_idx = std::floor(zp);

    for (auto i : {x_idx - 1, x_idx, x_idx + 1, x_idx + 2}) {
      for (auto j : {y_idx - 1, y_idx, y_idx + 1, y_idx + 2}) {
        for (auto k : {z_idx - 1, z_idx, z_idx + 1, z_idx + 2}) {
          if (i < 0 || j < 0 || k < 0) continue;

          if (i >= xsize || j >= ysize || k >= zsize) continue;

          const float x_c = conv_kernel(xp - i);
          const float y_c = conv_kernel(yp - j);
          const float z_c = conv_kernel(zp - k);
          v += g[i][j][k] * x_c * y_c * z_c;
        }
      }
    }

    return v;
  }

  //------------------------------------------------------------
  float Interpolator::conv_kernel(float s) const
  {
    float v = 0;
    float z = std::abs(s);

    if (0 <= z && z < 1)
      v = 1 + (0.5) * z * z * (3 * z - 5);

    else if (1 < z && z < 2)
      v = 2 - z * (4 + 0.5 * z * (z - 5));

    return v;
  }

} // namespace mag
