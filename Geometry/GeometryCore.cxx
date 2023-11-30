/**
 * @file   GeometryCore.cxx
 * @brief  Access the description of detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    GeometryCore.h
 *
 */

// class header
#include "Geometry/GeometryCore.h"
#include "Geometry/ChannelMapAlgs/ChannelMapAlg.h"
#include "Geometry/ChannelMapAlgs/SegmentationAlg.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoPgon.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"
#include "TGeoCompositeShape.h"

// C/C++ includes
#include <cstddef> // size_t
#include <cctype> // ::tolower()
#include <cmath> // std::abs() ...
#include <vector>
#include <algorithm> // std::for_each(), std::transform()
#include <utility> // std::swap()
#include <limits> // std::numeric_limits<>
#include <memory> // std::default_deleter<>
#include <regex>
#include <csignal>

#include <boost/format.hpp>

namespace gar {
  namespace geo {

    template <typename T>
    inline T sqr(T v) { return v * v; }

    //......................................................................
    // Constructor.  Class variables that can't come from gdml get initialized here
    //               from Geometry.fcl
    GeometryCore::GeometryCore(fhicl::ParameterSet const& pset)
      : fSurfaceY         (pset.get<double     >("SurfaceY"                ))
      , fDetectorName     (pset.get<std::string>("Name"                    ))
      , fPositionWiggle   (pset.get<double     >("PositionEpsilon",   1.e-4))
      , fPointInWarnings  (pset.get<bool       >("PointInWarnings",   false))
      , fECALEndcapOutside(pset.get<bool       >("ECALEndcapOutsidePV",true))
      , fIndsOfRefract    (pset.get<fhicl::ParameterSet>
                                                ("IndsOfRefract",fhicl::ParameterSet()))
   {
      std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);

      InitVariables();

      fIofRECAL    = fIndsOfRefract.get("ECALindex",    1.59);
      fIofRMuID    = fIndsOfRefract.get("MuIDindex",    1.59);
      fIofRMinerva = fIndsOfRefract.get("Minervaindex", 1.59);
    } // GeometryCore::GeometryCore()


    //......................................................................
    GeometryCore::~GeometryCore()
    {
      ClearGeometry();
    } // GeometryCore::~GeometryCore()


    //......................................................................
    void GeometryCore::ApplyChannelMap(std::shared_ptr<geo::seg::ChannelMapAlg> pChannelMap)
    {
      pChannelMap->Initialize(*this);
      fChannelMapAlg = pChannelMap;
    } // GeometryCore::ApplyChannelMap()

    //......................................................................
    void GeometryCore::ApplyECALSegmentationAlg(std::shared_ptr<geo::seg::SegmentationAlg> pECALSegmentationAlg)
    {
      pECALSegmentationAlg->Initialize(*this);
      fECALSegmentationAlg = pECALSegmentationAlg;
    } // GeometryCore::ApplyECALSegmentationAlg()

    //......................................................................
    void GeometryCore::ApplyMinervaSegmentationAlg(std::shared_ptr<geo::seg::SegmentationAlg> pMinervaSegmentationAlg)
    {
      pMinervaSegmentationAlg->Initialize(*this);
      fMinervaSegmentationAlg = pMinervaSegmentationAlg;
    } // GeometryCore::ApplyMinervaSegmentationAlg()

    //......................................................................
    void GeometryCore::ApplyMuIDSegmentationAlg(std::shared_ptr<geo::seg::SegmentationAlg> pMuIDSegmentationAlg)
    {
      pMuIDSegmentationAlg->Initialize(*this);
      fMuIDSegmentationAlg = pMuIDSegmentationAlg;
    } // GeometryCore::ApplyMuIDSegmentationAlg()

    //......................................................................
    void GeometryCore::LoadGeometryFile(std::string const& gdmlfile, std::string const& rootfile, bool bForceReload /* = false*/)
    {
      if (gdmlfile.empty()) {
        throw cet::exception("GeometryCore")
          << "No GDML Geometry file specified!\n";
      }

      if (rootfile.empty()) {
        throw cet::exception("GeometryCore")
          << "No ROOT Geometry file specified!\n";
      }

      ClearGeometry();

      // Open the GDML file, and convert it into ROOT TGeoManager format.
      // Then lock the gGeoManager to prevent future imports, for example
      // in AuxDetGeometry
      if( !gGeoManager || bForceReload ){
        if (gGeoManager) TGeoManager::UnlockGeometry();
        else { // very first time (or so it should)
          // [20210630, petrillo@slac.stanford.edu]
          // ROOT 6.22.08 allows us to choose the representation of lengths
          // in the geometry objects parsed from GDML.
          // In LArSoft we want them to be centimeters (ROOT standard).
          // This was tracked as Redmine issue #25990, and I leave this mark
          // because I feel that we'll be back to it not too far in the future.
          // Despite the documentation (ROOT 6.22/08),
          // it seems the units are locked from the beginning,
          // so we unlock without prejudice.
          TGeoManager::LockDefaultUnits(false);
          TGeoManager::SetDefaultUnits(TGeoManager::kRootUnits);
          TGeoManager::LockDefaultUnits(true);
        }
        TGeoManager::Import(rootfile.c_str());
        gGeoManager->LockGeometry();
      }

      fGDMLfile = gdmlfile;
      fROOTfile = rootfile;

      MF_LOG_INFO("GeometryCore")
        << "New detector geometry loaded from "
        << "\n\t" << fROOTfile
        << "\n\t" << fGDMLfile;

    } // GeometryCore::LoadGeometryFile()

    //......................................................................
    float GeometryCore::GetSensVolumeThickness(const TVector3& point) const
    {
      TGeoNode *node = gGeoManager->FindNode(point.x(), point.y(), point.z());

      if(node) {
    return this->FindShapeSize(node)[2] * 2;
      }
      else{
    return 0.;
      }
    }

    //......................................................................
    const std::array<double, 3> GeometryCore::FindShapeSize(const TGeoNode *node) const
    {
      TGeoVolume *vol = node->GetVolume();
      std::string volname = vol->GetName();

      //Check if it is ECAL endcap -> layer size is not the BBox! It is the apothem
      bool isBarrel = true;
      bool isECAL = false;

      if (volname.find("barrel") == std::string::npos &&
          volname.find("Barrel") == std::string::npos) {
        isBarrel =  false;
      }
      if (volname.find("PV") != std::string::npos) {
        isBarrel =  false;
      }

      if (volname.find("ECal") != std::string::npos || volname.find("ecal") != std::string::npos)
        isECAL = true;

      // Old code commented out here but keep it a bit for reference
      // if(volname.find("endcap") != std::string::npos || volname.find("Endcap") != std::string::npos ) isBarrel = false;

      std::array<double, 3> shape = {0., 0., 0.};

      if (vol) {
        TGeoBBox *box = (TGeoBBox*)(vol->GetShape());

        if(isBarrel) {
          shape[0] = box->GetDX();
          shape[1] = box->GetDY();
        } else {
          if (fECALEndcapOutside && isECAL) {
            shape[0] = GetECALEndcapApothemLength() / 2.;
            shape[1] = GetECALEndcapApothemLength() / 2.;
          } else {
            shape[0] = box->GetDX();
            shape[1] = box->GetDY();
          }
        }

        shape[2] = box->GetDZ();

        return shape; //return half size in cm
      } else {
        throw cet::exception("GeometryCore::FindShapeSize")
          << "Could not find the volume associated to node "
          << node->GetName() <<"\n";
      }
    }

    //......................................................................
    void GeometryCore::ClearGeometry() {

      return;
    } // GeometryCore::ClearGeometry()


    //......................................................................
    TGeoManager* GeometryCore::ROOTGeoManager() const {
      return gGeoManager;
    }

    //......................................................................
    unsigned int GeometryCore::NChannels() const {
      return fChannelMapAlg->Nchannels();
    }

    //......................................................................
    struct NodeNameMatcherClass {
      std::set<std::string> const* vol_names;

      NodeNameMatcherClass(std::set<std::string> const& names)
    : vol_names(&names) {}

      /// Returns whether the specified node matches a set of names
      bool operator() (TGeoNode const& node) const
      {
        if (!vol_names) return true;
        return vol_names->find(node.GetVolume()->GetName()) != vol_names->end();
      }

    }; // NodeNameMatcherClass

    //......................................................................
    struct CollectNodesByName {
      std::vector<TGeoNode const*> nodes;

      CollectNodesByName(std::set<std::string> const& names): matcher(names) {}

      /// If the name of the node matches, records the end node
      void operator() (TGeoNode const& node)
      { if (matcher(node)) nodes.push_back(&node); }

      void operator() (ROOTGeoNodeForwardIterator const& iter)
      { operator() (**iter); }

      protected:
        NodeNameMatcherClass matcher;
    }; // CollectNodesByName

    //......................................................................
    struct CollectPathsByName {
      std::vector<std::vector<TGeoNode const*>> paths;

      CollectPathsByName(std::set<std::string> const& names): matcher(names) {}

      /// If the name of the node matches, records the node full path
      void operator() (ROOTGeoNodeForwardIterator const& iter)
      { if (matcher(**iter)) paths.push_back(iter.get_path()); }

      protected:
        NodeNameMatcherClass matcher;
    }; // CollectPathsByName

    //......................................................................
    std::vector<TGeoNode const*> GeometryCore::FindVolumePath(std::string const& vol_name) const
    {
      std::vector<TGeoNode const*> path = { ROOTGeoManager()->GetTopNode() };
      if (!FindFirstVolume(vol_name, path)) path.clear();

      return path;
    } // GeometryCore::FindVolumePath()

    //......................................................................
    bool GeometryCore::FindFirstVolume(std::string const& name, std::vector<const TGeoNode*>& path) const
    {
      assert(!path.empty());
      auto const* pCurrent = path.back();
      auto const* pCurrentVolume = pCurrent->GetVolume();

      // first check the current layer
      if (strncmp(name.c_str(), pCurrentVolume->GetName(), name.length()) == 0)
        return true;

      //explore the next layer down
      unsigned int nd = pCurrentVolume->GetNdaughters();
      for(unsigned int i = 0; i < nd; ++i) {
        path.push_back(pCurrentVolume->GetNode(i));
        if (FindFirstVolume(name, path)) return true;
        path.pop_back();
      } // for
      return false;
    } // GeometryCore::FindFirstVolume()

    //......................................................................
    void GeometryCore::StoreECALNodes(std::map<std::string, std::vector<const TGeoNode*>> &map) const {

      //Loop over all nodes for the ecal active volume and store it in a map
      //Active volume has the form [Barrel-Endcap]ECal_staveXX_moduleYY_layer_ZZ_slice2_vol
      std::string det[2] = { "Barrel", "Endcap" };
      unsigned int nstave = GetECALInnerSymmetry() + 1;
      unsigned int nmodule = 7;
      unsigned int nlayer = fECALSegmentationAlg->nLayers() + 1;

      for(unsigned int idet = 0; idet < 2; idet++) {
        for(unsigned int istave = 0; istave < nstave; istave++) {
          for(unsigned int imodule = 0; imodule < nmodule; imodule++) {
            for(unsigned int ilayer = 0; ilayer < nlayer; ilayer++) {
              boost::format bname = boost::format("%sECal_stave%02i_module%02i_layer_%02i_slice2_vol") % det[idet].c_str() % istave % imodule % ilayer;
              std::string vol_name = bname.str();
              std::vector<const TGeoNode*> path = FindVolumePath(vol_name);
              if(!path.empty()) map.emplace(vol_name, path); //insert in the map
            }
          }
        }
      }

      return;
    } // GeometryCore::StoreECALNodes()

    //......................................................................
    std::vector<TGeoNode const*> GeometryCore::FindAllVolumes(std::set<std::string> const& vol_names) const
    {
      CollectNodesByName node_collector(vol_names);

      ROOTGeoNodeForwardIterator iNode(ROOTGeoManager()->GetTopNode());
      TGeoNode const* pCurrentNode;

      while ((pCurrentNode = *iNode)) {
        node_collector(*pCurrentNode);
        ++iNode;
      } // while

      return node_collector.nodes;
    } // GeometryCore::FindAllVolumes()

    //......................................................................
    std::vector<std::vector<TGeoNode const*>> GeometryCore::FindAllVolumePaths(std::set<std::string> const& vol_names) const
    {
      CollectPathsByName path_collector(vol_names);

      ROOTGeoNodeForwardIterator iNode(ROOTGeoManager()->GetTopNode());

      while (*iNode) {
        path_collector(iNode);
        ++iNode;
      } // while

      return path_collector.paths;
    } // GeometryCore::FindAllVolumePaths()

    //......................................................................
    template<>
    TGeoNode* GeometryCore::FindNode<float>(float const &x, float const &y, float const &z) const
    {
      return gGeoManager->FindNode(x, y, z);
    }

    //......................................................................
    template<>
    TGeoNode* GeometryCore::FindNode<double>(double const &x, double const &y, double const &z) const
    {
      return gGeoManager->FindNode(x, y, z);
    }

    //......................................................................
    TGeoNode* GeometryCore::FindNode(std::array<double, 3> const& point) const
    {
      return gGeoManager->FindNode(point[0], point[1], point[2]);
    }

    //......................................................................
    TGeoNode* GeometryCore::FindNode(TVector3 const& point) const
    {
      return gGeoManager->FindNode(point.x(), point.y(), point.z());
    }

    //......................................................................
    bool GeometryCore::WorldToLocal(std::array<double, 3> const& world, std::array<double, 3> &local, gar::geo::LocalTransformation<TGeoHMatrix> &trans) const
    {
      TVector3 point(world[0], world[1], world[2]);
      std::string name = VolumeName(point);
      std::vector<const TGeoNode*> path;

      // std::cout << "WorldToLocal -- Finding volume " << name << " ..." << std::endl;
      if( fECALNodePath.find(name) != fECALNodePath.end() ) {
        // std::cout << "Found volume " << name << " in fECALNodePath" << std::endl;
        path = fECALNodePath.at(name);
      } else {
        path = FindVolumePath(name);
      }

      if (path.empty()){
        throw cet::exception("GeometryCore::WorldToLocal") << " can't find volume '" << name << "'\n";
        return false;
      }

      //Change to local frame
      trans.SetPath(path, path.size() - 1);
      std::array<double, 3> wd{ {world[0], world[1], world[2]} }, loc;
      trans.WorldToLocal(wd.data(), loc.data());
      local = {loc.at(0), loc.at(1), loc.at(2)};

      return true;
    }

    //......................................................................
    bool GeometryCore::LocalToWorld(std::array<double, 3> const& local, std::array<double, 3> &world, gar::geo::LocalTransformation<TGeoHMatrix> const &trans) const
    {
      if (trans.GetNodes().empty()){
        throw cet::exception("GeometryCore::LocalToWorld")
          << " LocalTransformation has no nodes! -- Call WorldToLocal first!" << "\n";
        return false;
      }

      //Change to world frame
      std::array<double, 3> loc{ {local[0], local[1], local[2]} }, wd;
      trans.LocalToWorld(loc.data(), wd.data());
      world = {wd.at(0), wd.at(1), wd.at(2)};

      return true;
    }

    //......................................................................
    //
    // Return the ranges of x,y and z for the "world volume" that the
    // entire geometry lives in. If any pointers are 0, then those
    // coordinates are ignored.
    //
    // \param xlo : On return, lower bound on x positions
    // \param xhi : On return, upper bound on x positions
    // \param ylo : On return, lower bound on y positions
    // \param yhi : On return, upper bound on y positions
    // \param zlo : On return, lower bound on z positions
    // \param zhi : On return, upper bound on z positions
    //
    void GeometryCore::WorldBox(float* xlo, float* xhi,
                float* ylo, float* yhi,
                float* zlo, float* zhi) const
    {
      const TGeoShape* s = gGeoManager->GetVolume("volWorld")->GetShape();
      if (!s)
        throw cet::exception("GeometryCore") << "no pointer to world volume TGeoShape\n";

      if (xlo || xhi) {
        double x1, x2;
        s->GetAxisRange(1,x1,x2); if (xlo) *xlo = x1; if (xhi) *xhi = x2;
      }
      if (ylo || yhi) {
        double y1, y2;
        s->GetAxisRange(2,y1,y2); if (ylo) *ylo = y1; if (yhi) *yhi = y2;
      }
      if (zlo || zhi) {
        double z1, z2;
        s->GetAxisRange(3,z1,z2); if (zlo) *zlo = z1; if (zhi) *zhi = z2;
      }
    }

    //......................................................................
    bool GeometryCore::PointInWorld(TVector3 const& point) const
    {
      // check that the given point is in the World volume at least
      TGeoVolume *volWorld = gGeoManager->FindVolumeFast("volWorld");
      float halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
      float halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
      float halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
      if (std::abs(point.x()) > halfwidth  ||
          std::abs(point.y()) > halfheight ||
          std::abs(point.z()) > halflength){

        if (fPointInWarnings) {
          MF_LOG_WARNING("GeometryCoreBadInputPoint")
            << "point ("
            << point.x() << ","
            << point.y() << ","
            << point.z() << ") "
            << "is not inside the world volume "
            << " half width = "  << halfwidth
            << " half height = " << halfheight
            << " half length = " << halflength;
        }
        return false;
      }

      return true;
    }

    //......................................................................
    bool GeometryCore::PointInDetEnclosure(TVector3 const& point) const
    {
      // check that the given point is in the enclosure volume at least
      if (std::abs(point.x()) > fEnclosureHalfWidth  ||
          std::abs(point.y()) > fEnclosureHalfHeight ||
          std::abs(point.z()) > fEnclosureLength){

        if (fPointInWarnings) {
          MF_LOG_WARNING("GeometryCoreBadInputPoint")
            << "point ("
            << point.x() << ","
            << point.y() << ","
            << point.z() << ") "
            << "is not inside the detector enclosure volume "
            << " half width = "  << fEnclosureHalfWidth
            << " half height = " << fEnclosureHalfHeight
            << " half length = " << fEnclosureLength;
        }
        return false;
      }

      return true;
    }

    //......................................................................
    bool GeometryCore::PointInMPD(TVector3 const& point) const
    {
      TVector3 tpc_origin(GetOriginX(), GetOriginY(), GetOriginZ());
      TVector3 new_point = point - tpc_origin;
      // check that the given point is in the enclosure volume at least
      if (std::abs(new_point.x()) > fMPDHalfWidth  ||
          std::abs(new_point.y()) > fMPDHalfHeight ||
          std::abs(new_point.z()) > fMPDLength){

        if (fPointInWarnings) {
          MF_LOG_WARNING("GeometryCoreBadInputPoint")
          << "point ("
          << point.x() << ","
          << point.y() << ","
          << point.z() << ") "
          << "is not inside the MPD volume "
          << " half width = "  << fMPDHalfWidth
          << " half height = " << fMPDHalfHeight
          << " half length = " << fMPDLength;
        }
        return false;
      }

      return true;
    }

    //......................................................................
    bool GeometryCore::PointInGArTPC(TVector3 const& point) const
    {
      // check that the given point is in the enclosure volume at least
      float y = std::abs(point.y() - GetOriginY());
      float z = std::abs(point.z() - GetOriginZ());
      if (std::abs(point.x() - GetOriginX()) > fTPCLength/2.0  ||
          std::hypot(z,y)                    > fTPCRadius) {
        if (fPointInWarnings) {
          MF_LOG_WARNING("GeometryCoreBadInputPoint")
            << "point ("
            << point.x() << ","
            << point.y() << ","
            << point.z() << ") "
            << "is not inside the GArTPC volume "
            << " radius = " << fTPCRadius
            << " length = " << fTPCLength;
        }
        return false;
      }

      return true;
    }

    //......................................................................
    bool GeometryCore::PointInLArTPC(TVector3 const& point) const
    {
      //Name to change depending on the detector..
      TGeoVolume *volLArTPC = gGeoManager->FindVolumeFast("volArgonCubeActive");

      float halflength = ((TGeoBBox*)volLArTPC->GetShape())->GetDZ();
      float halfheight = ((TGeoBBox*)volLArTPC->GetShape())->GetDY();
      float halfwidth  = ((TGeoBBox*)volLArTPC->GetShape())->GetDX();
      if (std::abs(point.x()) > halfwidth  ||
      std::abs(point.y()) > halfheight ||
      std::abs(point.z()) > halflength){
        if (fPointInWarnings) {
           MF_LOG_WARNING("GeometryCoreBadInputPoint")
              << "point ("
              << point.x() << ","
              << point.y() << ","
              << point.z() << ") "
              << "is not inside the LArTPC volume "
              << " half width = "  << halfwidth
              << " half height = " << halfheight
              << " half length = " << halflength;
        }
        return false;
      }

      return true;
    }

    //......................................................................
    bool GeometryCore::PointInECALBarrel(TVector3 const& point) const
    {
      std::string vol_name = this->VolumeName(point);

      if (vol_name.find("barrel") == std::string::npos &&
          vol_name.find("Barrel") == std::string::npos) {
        return false;
      }

      if (vol_name.find("ECal") != std::string::npos || 
          vol_name.find("ECAL") != std::string::npos || 
          vol_name.find("ecal") != std::string::npos) return true;

      return false;
    }

    //......................................................................
    bool GeometryCore::PointInECALEndcap(TVector3 const& point) const
    {
      std::string vol_name = this->VolumeName(point);

      if (vol_name.find("endcap") == std::string::npos &&
          vol_name.find("Endcap") == std::string::npos) {
        return false;
      }

      if (vol_name.find("ECal") != std::string::npos || 
          vol_name.find("ECAL") != std::string::npos || 
          vol_name.find("ecal") != std::string::npos) return true;

      return false;
    }

    //......................................................................
    bool GeometryCore::PointInMuIDBarrel(TVector3 const& point) const
    {
      std::string vol_name = this->VolumeName(point);

      if (vol_name.find("barrel") == std::string::npos &&
          vol_name.find("Barrel") == std::string::npos) {
        return false;
      }

      if (vol_name.find("Yoke") != std::string::npos || 
          vol_name.find("yoke") != std::string::npos) return true;

      return false;
    }

    //......................................................................
    bool GeometryCore::PointInMuIDEndcap(TVector3 const& point) const
    {
      std::string vol_name = this->VolumeName(point);

      if (vol_name.find("endcap") == std::string::npos &&
          vol_name.find("Endcap") == std::string::npos) {
        return false;
      }

      if (vol_name.find("Yoke") != std::string::npos || 
          vol_name.find("yoke") != std::string::npos) return true;

      return false;
    }

    //......................................................................
    const std::string GeometryCore::VolumeName(TVector3 const& point) const
    {
      if( !this->PointInWorld(point) ){
    const std::string unknown("unknownVolume");
    return unknown;
      }

      const std::string name(this->FindNode(point)->GetVolume()->GetName());
      return name;
    }

    //......................................................................
    const std::string GeometryCore::MaterialName(TVector3 const& point)
    {
      if( !this->PointInWorld(point) ){
    const std::string unknown("unknownVolume");
    return unknown;
      }

      const std::string name(this->FindNode(point)->GetMedium()->GetMaterial()->GetName());
      return name;
    }

    //......................................................................
    //
    // Return the total mass of the detector
    //
    //
    double GeometryCore::TotalMass(const char *vol) const
    {
      //the TGeoNode::GetVolume() returns the TGeoVolume of the detector outline
      //and ROOT calculates the mass in kg for you
      TGeoVolume *gvol = gGeoManager->FindVolumeFast(vol);
      if(gvol) return gvol->Weight();

      throw cet::exception("GeometryCore") << "could not find specified volume "
                       << vol
                       << " to determine total mass\n";
    }

    //......................................................................
    //
    // Return the column density between 2 points
    //
    // \param p1  : pointer to array holding xyz of first point in world coordinates
    // \param p2  : pointer to array holding xyz of second point in world coorinates
    //
    double GeometryCore::MassBetweenPoints(double *p1, double *p2) const
    {

      //The purpose of this method is to determine the column density
      //between the two points given.  Do that by starting at p1 and
      //stepping until you get to the node of p2.  calculate the distance
      //between the point just inside that node and p2 to get the last
      //bit of column density
      double columnD = 0.;

      //first initialize a track - get the direction cosines
      double length = std::sqrt( sqr(p2[0]-p1[0])
                 + sqr(p2[1]-p1[1])
                 + sqr(p2[2]-p1[2]));
      double dxyz[3] = {(p2[0]-p1[0])/length, (p2[1]-p1[1])/length, (p2[2]-p1[2])/length};

      gGeoManager->InitTrack(p1,dxyz);

      //might be helpful to have a point to a TGeoNode
      TGeoNode *node = gGeoManager->GetCurrentNode();

      //check that the points are not in the same volume already.
      //if they are in different volumes, keep stepping until you
      //are in the same volume as the second point
      while(!gGeoManager->IsSameLocation(p2[0], p2[1], p2[2])){
    gGeoManager->FindNextBoundary();
    columnD += gGeoManager->GetStep()*node->GetMedium()->GetMaterial()->GetDensity();

    //the act of stepping puts you in the next node and returns that node
    node = gGeoManager->Step();
      }//end loop to get to volume of second point

      //now you are in the same volume as the last point, but not at that point.
      //get the distance between the current point and the last one
      const double *current = gGeoManager->GetCurrentPoint();
      length = std::sqrt(sqr(p2[0]-current[0]) +
             sqr(p2[1]-current[1]) +
             sqr(p2[2]-current[2]));
      columnD += length*node->GetMedium()->GetMaterial()->GetDensity();

      return columnD;
    }

    //--------------------------------------------------------------------
    float GeometryCore::GetIROCInnerRadius() const
    {
      return fChannelMapAlg->GetIROCInnerRadius();
    }

    //--------------------------------------------------------------------
    float GeometryCore::GetIROCOuterRadius() const
    {
      return fChannelMapAlg->GetIROCOuterRadius();
    }

    //--------------------------------------------------------------------
    float GeometryCore::GetOROCInnerRadius() const
    {
      return fChannelMapAlg->GetOROCInnerRadius();
    }

    //--------------------------------------------------------------------
    float GeometryCore::GetOROCPadHeightChangeRadius() const
    {
      return fChannelMapAlg->GetOROCPadHeightChangeRadius();
    }

    //--------------------------------------------------------------------
    float GeometryCore::GetOROCOuterRadius() const
    {
      return fChannelMapAlg->GetOROCOuterRadius();
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::NearestChannel(const float worldLoc[3]) const
    {
      return fChannelMapAlg->NearestChannel(worldLoc);
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::NearestChannel(std::vector<float> const& worldLoc) const
    {
      float loc[3] = {worldLoc[0], worldLoc[1], worldLoc[2]};

      return this->NearestChannel(loc);
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::NearestChannel(TVector3 const& worldLoc) const
    {
      float loc[3] = {(float)worldLoc[0], (float)worldLoc[1], (float)worldLoc[2]};

      return this->NearestChannel(loc);
    }


    //--------------------------------------------------------------------
    void GeometryCore::NearestChannelInfo(float const* xyz, gar::geo::ChanWithNeighbors &cwn) const
    {
      fChannelMapAlg->NearestChannelInfo(xyz, cwn);
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::GapChannelNumber() const
    {
      return fChannelMapAlg->GapChannelNumber();
    }


    //--------------------------------------------------------------------
    void GeometryCore::ChannelToPosition(unsigned int const channel, float* const worldLoc) const
    {
      return fChannelMapAlg->ChannelToPosition(channel, worldLoc);
    }

    //----------------------------------------------------------------------------
    gar::raw::CellID_t GeometryCore::GetCellID(const TGeoNode *node, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
    {
      std::string node_name = node->GetName();
      gar::raw::CellID_t cellID = 0.;

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {

    const std::array<double, 3> shape = this->FindShapeSize(node);
    fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
    cellID = fECALSegmentationAlg->GetCellID(*this, det_id, stave, module, layer, slice, localPosition);

      } else if(node_name.find("TrackerSc") != std::string::npos || node_name.find("trackersc") != std::string::npos || node_name.find("Tracker") != std::string::npos || node_name.find("tracker") != std::string::npos) {

    const std::array<double, 3> shape = this->FindShapeSize(node);
    fMinervaSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
    cellID = fMinervaSegmentationAlg->GetCellID(*this, det_id, 0, 0, layer, slice, localPosition);

      } else if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {

    const std::array<double, 3> shape = this->FindShapeSize(node);
    fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
    cellID = fMuIDSegmentationAlg->GetCellID(*this, det_id, stave, module, layer, slice, localPosition);

      } else {
    MF_LOG_WARNING("GeometryCore::GetCellID")
      << "Detector id " << det_id << " unknown!"
      << " Node name " << node_name;
      }

      return cellID;
    }

    //----------------------------------------------------------------------------
    const std::string GeometryCore::GetECALCellIDEncoding() const
    {
      if(fECALSegmentationAlg)
        return fECALSegmentationAlg->cellEncoding();
      else
        return "";
    }

    //----------------------------------------------------------------------------
    const std::string GeometryCore::GetMinervaCellIDEncoding() const
    {
      if(fMinervaSegmentationAlg)
        return fMinervaSegmentationAlg->cellEncoding();
      else
        return "";
    }

    //----------------------------------------------------------------------------
    const std::string GeometryCore::GetMuIDCellIDEncoding() const
    {
      if(fMuIDSegmentationAlg)
        return fMuIDSegmentationAlg->cellEncoding();
      else
        return "";
    }

    //----------------------------------------------------------------------------
    std::array<double, 3> GeometryCore::GetPosition(const TGeoNode *node, const gar::raw::CellID_t &cID) const
    {
      std::string node_name = node->GetName();
      std::array<double, 3> pos;

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        const std::array<double, 3> shape = this->FindShapeSize(node);
        fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fECALSegmentationAlg->setVariables(GetECALInnerAngle(), GetECALEndcapSideLength());
        pos = fECALSegmentationAlg->GetPosition(*this, cID);
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        const std::array<double, 3> shape = this->FindShapeSize(node);
        fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fMuIDSegmentationAlg->setVariables(GetMuIDInnerAngle(), 0.);
        pos = fMuIDSegmentationAlg->GetPosition(*this, cID);
      }

      return pos;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::isTile(const std::array<double, 3>& point, const gar::raw::CellID_t& cID) const
    {
      bool isTile = false;
      std::string node_name = this->FindNode(point)->GetName();

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        isTile = fECALSegmentationAlg->isTile(cID);
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        isTile = fMuIDSegmentationAlg->isTile(cID);
      }

      return isTile;
    }

    //----------------------------------------------------------------------------
    double GeometryCore::getStripWidth(const std::array<double, 3>& point) const
    {
      double strip_width = 0.;
      std::string node_name = this->FindNode(point)->GetName();

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        strip_width = fECALSegmentationAlg->stripSizeX();
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        strip_width = fMuIDSegmentationAlg->stripSizeX();
      }

      return strip_width;
    }

    //----------------------------------------------------------------------------
    double GeometryCore::getTileSize(const std::array<double, 3>& point) const
    {
      double tilesize = 0.;
      std::string node_name = this->FindNode(point)->GetName();

      if (node_name.find("ECal") != std::string::npos || 
          node_name.find("ECAL") != std::string::npos || 
          node_name.find("ecal") != std::string::npos) {

          tilesize = fECALSegmentationAlg->gridSizeX();
      }

      return tilesize;
    }

    //----------------------------------------------------------------------------
    double GeometryCore::getStripLength(const std::array<double, 3>& point, const gar::raw::CellID_t &cID) const
    {
      double strip_length = 0.;
      const TGeoNode *node = this->FindNode(point);
      std::string node_name = node->GetName();

      std::array<double, 3> localtemp;
      gar::geo::LocalTransformation<TGeoHMatrix> trans;
      this->WorldToLocal(point, localtemp, trans);
      const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

      if (node_name.find("ECal") != std::string::npos || 
          node_name.find("ECAL") != std::string::npos || 
          node_name.find("ecal") != std::string::npos) {

          fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
          fECALSegmentationAlg->setVariables(GetECALInnerAngle(), GetECALEndcapSideLength());
          strip_length = fECALSegmentationAlg->getStripLength(*this, localtemp, cID);
      }

      if (node_name.find("Yoke") != std::string::npos || 
          node_name.find("yoke") != std::string::npos) {

          fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
          fMuIDSegmentationAlg->setVariables(GetMuIDInnerAngle(), 0.);
          strip_length = fMuIDSegmentationAlg->getStripLength(*this, localtemp, cID);
      }

      return strip_length;
    }

    //----------------------------------------------------------------------------
    std::pair<TVector3, TVector3> GeometryCore::GetStripEnds(const std::array<double, 3>& point, const gar::raw::CellID_t &cID) const
    {
      std::string node_name = this->FindNode(point)->GetName();
      std::pair<TVector3, TVector3> localStripEnds;

      //Get the matrix to make the transformation from Local to World
      std::array<double, 3> localtemp;
      gar::geo::LocalTransformation<TGeoHMatrix> trans;
      this->WorldToLocal(point, localtemp, trans);
      const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fECALSegmentationAlg->setVariables(GetECALInnerAngle(), GetECALEndcapSideLength());
        localStripEnds = fECALSegmentationAlg->getStripEnds(*this, localtemp, cID);
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fMuIDSegmentationAlg->setVariables(GetMuIDInnerAngle(), 0.);
        localStripEnds = fMuIDSegmentationAlg->getStripEnds(*this, localtemp, cID);
      }

      //Get the world coordinates from both local coordinates of the strip ends
      std::array<double, 3> stripEnd1local = { localStripEnds.first.X(), localStripEnds.first.Y(), localStripEnds.first.Z() };
      std::array<double, 3> stripEnd2local = { localStripEnds.second.X(), localStripEnds.second.Y(), localStripEnds.second.Z() };
      std::array<double, 3> stripEnd1;
      std::array<double, 3> stripEnd2;
      this->LocalToWorld(stripEnd1local, stripEnd1, trans);
      this->LocalToWorld(stripEnd2local, stripEnd2, trans);

      return std::make_pair( TVector3(stripEnd1[0], stripEnd1[1], stripEnd1[2]), TVector3(stripEnd2[0], stripEnd2[1], stripEnd2[2]) );
    }

    //----------------------------------------------------------------------------
    std::pair<float, float> GeometryCore::CalculateLightPropagation(const std::array<double, 3>& point, const std::array<double, 3> &local, const gar::raw::CellID_t &cID) const
    {
      std::pair<float, float> light_prop;
      std::string node_name = this->FindNode(point)->GetName();
      const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fECALSegmentationAlg->setVariables(GetECALInnerAngle(), GetECALEndcapSideLength());
        light_prop = fECALSegmentationAlg->CalculateLightPropagation(*this, local, cID);
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fMuIDSegmentationAlg->setVariables(GetMuIDInnerAngle(), 0.);
        light_prop = fMuIDSegmentationAlg->CalculateLightPropagation(*this, local, cID);
      }

      return light_prop;
    }

    //----------------------------------------------------------------------------
    std::array<double, 3> GeometryCore::ReconstructStripHitPosition(const std::array<double, 3>& point, const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t &cID) const
    {
      std::array<double, 3> pos;
      std::string node_name = this->FindNode(point)->GetName();
      const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

      if(node_name.find("ECal") != std::string::npos || node_name.find("ECAL") != std::string::npos || node_name.find("ecal") != std::string::npos) {
        fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fECALSegmentationAlg->setVariables(GetECALInnerAngle(), GetECALEndcapSideLength());
        pos = fECALSegmentationAlg->ReconstructStripHitPosition(*this, local, xlocal, cID);
      }

      if(node_name.find("Yoke") != std::string::npos || node_name.find("yoke") != std::string::npos) {
        fMuIDSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);
        fMuIDSegmentationAlg->setVariables(GetMuIDInnerAngle(), 0.);
        pos = fMuIDSegmentationAlg->ReconstructStripHitPosition(*this, local, xlocal, cID);
      }

      return pos;
    }

    //----------------------------------------------------------------------------
    void GeometryCore::GetGeometryParameters()
    {

      // move this ahead of the segmentation call in the service
      GetDetectorsPresent();

      if(fHasGasTPCDetector)
        StoreTPCParameters();
      if(fHasECALDetector)
        StoreECALParameters();
      if(fHasTrackerScDetector) {
        /* no op */
      }
      if(fHasMuonDetector)
        StoreMuIDParameters();

      StoreOtherParameters();

      SetDetectorOrigin(); //Set the origin!
    }

    //----------------------------------------------------------------------------
    void GeometryCore::FinalizeGeometryParameters() 
    {
      if(fHasECALDetector) {
        FindECALnLayers();
        MakeECALLayeredCalorimeterData();
        StoreECALNodes(fECALNodePath);
      }

      if(fHasMuonDetector) {
        FindMuIDnLayers();
      }

      if(fHasTrackerScDetector) {
        FindTrackerScnPlanes();
      }
    }

    //----------------------------------------------------------------------------
    void GeometryCore::GetDetectorsPresent()
    {
      FindWorldVolume();
      fHasRock = FindRockVolume();
      fHasEnclosure = FindEnclosureVolume();
      //MPD Mother volume
      FindMPDVolume();

      fHasLArTPCDetector = FindLArTPCVolume();

      fHasGasTPCDetector = FindGasTPCVolume();
      fHasECALDetector = FindECALVolume();

      fHasTrackerScDetector = FindTrackerScVolume();

      fHasMuonDetector = FindMuIDVolume();
    }

    //----------------------------------------------------------------------------
    void GeometryCore::SetDetectorOrigin()
    {
      //Case ND-GAr full
      if(fHasGasTPCDetector) {
        fOriginX = TPCXCent();
        fOriginY = TPCYCent();
        fOriginZ = TPCZCent();
      }

      //Case ND-GAr Lite
      if(fHasTrackerScDetector) {
        fOriginX = GArLiteXCent();
        fOriginY = GArLiteYCent();
        fOriginZ = GArLiteZCent();
      }

    }

    //----------------------------------------------------------------------------
    void GeometryCore::StoreOtherParameters()
    {
      if(fHasLArTPCDetector) {
        //For the LArTPC
        FindActiveLArTPCVolume();
      }
    }

    //----------------------------------------------------------------------------
    void GeometryCore::StoreTPCParameters()
    {
      if(!fHasGasTPCDetector) return;
      FindActiveTPCVolume();
    }

    //----------------------------------------------------------------------------
    void GeometryCore::StoreECALParameters()
    {
      FindECALInnerBarrelRadius();
      FindECALOuterBarrelRadius();
      FindECALInnerEndcapRadius();
      FindECALOuterEndcapRadius();
      FindPVThickness();
      FindTPCOuterX();
      FindECALInnerSymmetry();
      FindECALEndcapStartX();
      FindECALEndcapOuterX();
      // std::raise(SIGINT);
    }

    //----------------------------------------------------------------------------
    void GeometryCore::StoreMuIDParameters()
    {
      FindMuIDInnerBarrelRadius();
      FindMuIDOuterBarrelRadius();
      FindMuIDInnerSymmetry();
      // std::raise(SIGINT);
    }

    //----------------------------------------------------------------------------
    void GeometryCore::PrintGeometry() const
    {
      std::cout << "Geometry loaded with \n" << std::endl;
      std::cout << "------------------------------" << std::endl;
      std::cout << "ND-GAr detector origin set to " << GetOriginX() << " cm " << GetOriginY() << " cm " << GetOriginZ() << " cm" << std::endl;
        
      //Prints geometry parameters
      std::cout << "------------------------------" << std::endl;
      std::cout << "World Geometry" << std::endl;
      std::cout << "World Origin (x, y, z) " << GetWorldX() << " cm " << GetWorldY() << " cm " << GetWorldZ() << " cm" << std::endl;
      std::cout << "World Size (H, W, L) " << GetWorldHalfWidth() << " cm " << GetWorldHalfHeight() << " cm " << GetWorldLength() << " cm" << std::endl;

      if(this->HasRock()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Rock Geometry" << std::endl;
        std::cout << "Rock Origin (x, y, z) " << GetRockX() << " cm " << GetRockY() << " cm " << GetRockZ() << " cm" << std::endl;
        std::cout << "Rock Size (H, W, L) " << GetRockHalfWidth() << " cm " << GetRockHalfHeight() << " cm " << GetRockLength() << " cm" << std::endl;
      }

      if(this->HasEnclosure()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Enclosure Geometry" << std::endl;
        std::cout << "Enclosure Origin (x, y, z) " << GetEnclosureX() << " cm " << GetEnclosureY() << " cm " << GetEnclosureZ() << " cm" << std::endl;
        std::cout << "Enclosure Size (H, W, L) " << GetEnclosureHalfWidth() << " cm " << GetEnclosureHalfHeight() << " cm " << GetEnclosureLength() << " cm" << std::endl;
      }

      if(this->HasLArTPCDetector()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "LArArgonCube Geometry" << std::endl;
        std::cout << "LArTPC Origin (x, y, z) " << GetLArTPCX() << " cm " << GetLArTPCY() << " cm " << GetLArTPCZ() << " cm" << std::endl;
        std::cout << "LArTPC Size (H, W, L) " << GetLArTPCHalfWidth() << " cm " << GetLArTPCHalfHeight() << " cm " << GetLArTPCLength() << " cm" << std::endl;
        std::cout << "------------------------------" << std::endl;
        std::cout << "LArActiveArgonCube Geometry" << std::endl;
        std::cout << "LArTPCActive Origin (x, y, z) " << GetActiveLArTPCX() << " cm " << GetActiveLArTPCY() << " cm " << GetActiveLArTPCZ() << " cm" << std::endl;
        std::cout << "LArTPCActive Size (H, W, L) " << GetActiveLArTPCHalfWidth() << " cm " << GetActiveLArTPCHalfHeight() << " cm " << GetActiveLArTPCLength() << " cm" << std::endl;

      }

      std::cout << "------------------------------" << std::endl;
      std::cout << "ND-GAr Geometry" << std::endl;
      std::cout << "ND-GAr Origin (x, y, z) " << GetMPDX() << " cm " << GetMPDY() << " cm " << GetMPDZ() << " cm" << std::endl;
      std::cout << "ND-GAr Size (H, W, L) " << GetMPDHalfWidth() << " cm " << GetMPDHalfHeight() << " cm " << GetMPDLength() << " cm" << std::endl;

      if(this->HasGasTPCDetector()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "TPC Geometry" << std::endl;
        std::cout << "TPC Origin (x, y, z) " << TPCXCent() << " cm " << TPCYCent() << " cm " << TPCZCent() << " cm" << std::endl;
        std::cout << "TPC Active Volume Size (R, L) " << TPCRadius() << " cm " << TPCLength() << " cm" << std::endl;
        std::cout << "TPC Number of Drift Volumes " << TPCNumDriftVols() << std::endl;
      }

      if(this->HasECALDetector()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "ECAL Geometry" << std::endl;
        std::cout << "ECAL Barrel inner radius: " << GetECALInnerBarrelRadius() << " cm" << std::endl;
        std::cout << "ECAL Barrel outer radius: " << GetECALOuterBarrelRadius() << " cm" << std::endl;
        std::cout << "ECAL Endcap inside PV: " << !fECALEndcapOutside << std::endl;
        std::cout << "ECAL Endcap inner radius: " << GetECALInnerEndcapRadius() << " cm" << std::endl;
        std::cout << "ECAL Endcap outer radius: " << GetECALOuterEndcapRadius() << " cm" << std::endl;
        std::cout << "ECAL Barrel inner symmetry: " << GetECALInnerSymmetry() << std::endl;
        std::cout << "ECAL Barrel polyhedra angle: " << GetECALInnerAngle()*180/M_PI << " deg" << std::endl;
        std::cout << "ECAL Barrel polyhedra side length: " << GetECALBarrelSideLength() << " cm" << std::endl;
        std::cout << "ECAL Barrel polyhedra apothem length: " << GetECALBarrelApothemLength() << " cm" << std::endl;
        if(fECALEndcapOutside){
          std::cout << "ECAL Endcap polyhedra side length: " << GetECALEndcapSideLength() << " cm" << std::endl;
          std::cout << "ECAL Endcap polyhedra apothem length: " << GetECALEndcapApothemLength() << " cm" << std::endl;
        }
        std::cout << "ECAL Endcap Start X: " << GetECALEndcapStartX() << " cm" << std::endl;
        std::cout << "ECAL Endcap Outer X: " << GetECALEndcapOuterX() << " cm" << std::endl;
        std::cout << "Number of layers: " << GetNLayers("ECAL") << std::endl;
        std::cout << "Pressure Vessel Thickness: " << GetPVThickness() << " cm" << std::endl;
        std::cout << "Distance from TPC to ECAL Endcap: " << GetECALEndcapStartX() -GetTPCOuterX() << " cm" << std::endl;
      }

      if(this->HasTrackerScDetector()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "ND-GAr Lite Geometry" << std::endl;
        std::cout << "ND-GAr Lite Origin (x, y, z) " << GArLiteXCent() << " cm " << GArLiteYCent() << " cm " << GArLiteZCent() << " cm" << std::endl;
        std::cout << "ND-GAr Lite Active Volume Size (R, L) " << GArLiteRadius() << " cm " << GArLiteLength() << " cm" << std::endl;
        std::cout << "Number of planes (from Segmentation alg) " << GetNLayers("TrackerSc") << std::endl;
      }

      if(this->HasMuonDetector()) {
        std::cout << "------------------------------" << std::endl;
        std::cout << "MuID Geometry" << std::endl;
        std::cout << "MuID Barrel inner radius: " << GetMuIDInnerBarrelRadius() << " cm" << std::endl;
        std::cout << "MuID Barrel outer radius: " << GetMuIDOuterBarrelRadius() << " cm" << std::endl;
        std::cout << "MuID Barrel inner symmetry: " << GetMuIDInnerSymmetry() << std::endl;
        std::cout << "MuID Barrel polyhedra angle: " << GetMuIDInnerAngle()*180/M_PI << " deg" << std::endl;
        std::cout << "MuID Barrel polyhedra side length: " << GetMuIDBarrelSideLength() << " cm" << std::endl;
        std::cout << "MuID Barrel polyhedra apothem length: " << GetMuIDBarrelApothemLength() << " cm" << std::endl;
        std::cout << "Number of layers: " << GetNLayers("MuID") << std::endl;
      }
      std::cout << "------------------------------\n" << std::endl;
    }

    //......................................................................
    bool GeometryCore::FindWorldVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("volWorld");
      if(path.size() == 0) return false;

      const TGeoNode *world_node = path.at(path.size()-1);
      if(world_node == nullptr) {
        std::cout << "Cannot find node volWorld_0" << std::endl;
        return false;
      }

      const double *origin = world_node->GetMatrix()->GetTranslation();

      fWorldX = origin[0];
      fWorldY = origin[1];
      fWorldZ = origin[2];

      fWorldHalfWidth = ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDZ();
      fWorldHalfHeight = ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDY();
      fWorldLength = 2.0 * ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDX();

      return true;
    }

    //......................................................................
    bool GeometryCore::FindRockVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("rockBox_lv");
      if(path.size() == 0) return false;

      const TGeoNode *rock_node = path.at(path.size()-1);
      if(rock_node == nullptr) {
        std::cout << "Cannot find node rockBox_lv_0" << std::endl;
        return false;
      }

      const double *origin = rock_node->GetMatrix()->GetTranslation();

      fRockX = origin[0];
      fRockY = origin[1];
      fRockZ = origin[2];

      fRockHalfWidth = ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDZ();
      fRockHalfHeight = ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDY();
      fRockLength = 2.0 * ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDX();

      return true;
    }

    //......................................................................
    bool GeometryCore::FindEnclosureVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("volDetEnclosure");
      if(path.size() == 0) return false;

      const TGeoNode *enc_node = path.at(path.size()-1);
      if(enc_node == nullptr) {
        std::cout << "Cannot find node volDetEnclosure_0" << std::endl;
        return false;
      }

      const double *origin = enc_node->GetMatrix()->GetTranslation();

      fEnclosureX = origin[0];
      fEnclosureY = origin[1];
      fEnclosureZ = origin[2];

      fEnclosureHalfWidth = ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDZ();
      fEnclosureHalfHeight = ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDY();
      fEnclosureLength = 2.0 * ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDX();

      return true;
    }

    //......................................................................
    bool GeometryCore::FindLArTPCVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("volArgonCubeDetector");
      if(path.size() == 0)
        return false;

      const TGeoNode *lar_node = path.at(path.size()-1);
      if(lar_node == nullptr) {
        std::cout << "Cannot find node volArgonCubeDetector_0" << std::endl;
        return false;
      }

      const double *origin = lar_node->GetMatrix()->GetTranslation();

      fLArTPCX = origin[0];
      fLArTPCY = origin[1];
      fLArTPCZ = origin[2];

      fLArTPCHalfWidth = ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDZ();
      fLArTPCHalfHeight = ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDY();
      fLArTPCLength = 2.0 * ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDX();

      return true;
    }

    //......................................................................
    bool GeometryCore::FindMPDVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("volMPD");
      if(path.size() == 0)
        return false;

      const TGeoNode *mpd_node = path.at(path.size()-1);
      if(mpd_node == nullptr) {
        std::cout << "Cannot find node volMPD_0" << std::endl;
        return false;
      }

      const double *origin = mpd_node->GetMatrix()->GetTranslation();

      fMPDX = origin[0];
      fMPDY = origin[1];
      fMPDZ = origin[2];

      fMPDHalfWidth = ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDZ();
      fMPDHalfHeight = ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDY();
      fMPDLength = 2.0 * ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDX();

      return true;
    }


    //......................................................................
    bool GeometryCore::FindActiveLArTPCVolume()
    {
      // check if the current level of the detector is the active TPC volume, if
      // not, then dig a bit deeper
      std::vector<TGeoNode const*> path = this->FindVolumePath("volArgonCubeActive");
      if(path.size() == 0)
        return false;

      const TGeoNode *LArTPC_node = path.at(path.size()-1);
      if(LArTPC_node == nullptr) {
        std::cout << "Cannot find node volArgonCubeActive_0" << std::endl;
        return false;
      }

      TGeoMatrix *mat = LArTPC_node->GetMatrix();
      const double *origin = mat->GetTranslation();

      //Get the origin correctly
      fLArTPCXCent =  fWorldX + fRockX + fEnclosureX + fLArTPCX + origin[0];
      fLArTPCYCent =  fWorldY + fRockY + fEnclosureY + fLArTPCY + origin[1];
      fLArTPCZCent =  fWorldZ + fRockZ + fEnclosureZ + fLArTPCZ + origin[2];

      //Get the dimension of the active volume
      fLArTPCActiveHalfWidth = ((TGeoBBox*)LArTPC_node->GetVolume()->GetShape())->GetDZ();
      fLArTPCActiveHalfHeight = ((TGeoBBox*)LArTPC_node->GetVolume()->GetShape())->GetDY();
      fLArTPCActiveLength = 2.0 * ((TGeoBBox*)LArTPC_node->GetVolume()->GetShape())->GetDX();

      return true;
    }

    //......................................................................
    bool GeometryCore::FindActiveTPCVolume()
    {
      // check if the current level of the detector is the active TPC volume, if
      // not, then dig a bit deeper
      std::vector<TGeoNode const*> path = this->FindVolumePath("TPCGas_vol");
      if(path.size() == 0)
        path = this->FindVolumePath("volTPCGas");
      if(path.size() == 0)
        return false;

      const TGeoNode *GArTPC_node = path.at(path.size()-1);
      if(GArTPC_node == nullptr) {
        std::cout << "Cannot find node TPCGas_vol_0/TPCGas_vol_0" << std::endl;
        return false;
      }

      TGeoMatrix *mat = GArTPC_node->GetMatrix();
      const double *origin = mat->GetTranslation();

      //Get the origin correctly
      fTPCXCent =  fWorldX + fRockX + fEnclosureX + fMPDX + origin[0];
      fTPCYCent =  fWorldY + fRockY + fEnclosureY + fMPDY + origin[1];
      fTPCZCent =  fWorldZ + fRockZ + fEnclosureZ + fMPDZ + origin[2];

      //Get the dimension of the active volume
      TGeoVolume *activeVol = GArTPC_node->GetVolume();
      float rOuter = ((TGeoTube*)activeVol->GetShape())->GetRmax();

      fTPCRadius = rOuter;
      fTPCLength = 2.0 * ((TGeoTube*)activeVol->GetShape())->GetDZ();

      // suggest -- figure out how many TPC drift volumes there are and set fTPCNumDriftVols here.
      fTPCNumDriftVols = 2;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindGasTPCVolume()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("volGArTPC");
      if (!vol)
        return false;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALVolume()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volBarrelECal");
      if (!vol)
        return false;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindMuIDVolume()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("volYokeBarrel");
      if (!vol)
        return false;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALInnerBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volBarrelECal");
      if (!vol)
        return false;

      fECALRinner = ((TGeoPgon*)vol->GetShape())->GetRmin(0);

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALOuterBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volBarrelECal");
      if (!vol)
        return false;

      fECALRouter = ((TGeoPgon*)vol->GetShape())->GetRmax(0);

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALInnerEndcapRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("EndcapECal_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volEndcapECal");
      if (!vol)
        return false;

      fECALECapRinner = 0.;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALOuterEndcapRadius()
    {
      if(fECALEndcapOutside){
        fECALECapRouter = fECALRouter; //should be equal
      } else {
        //Inside the Pressure Vessel
        TGeoVolume *vol = gGeoManager->FindVolumeFast("EndcapECal_vol");
        if (!vol)
          vol = gGeoManager->FindVolumeFast("volEndcapECal");
        if (!vol)
          return false;

        fECALECapRouter = ((TGeoBBox*)vol->GetShape())->GetDX();
      }

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindPVThickness()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("PVBarrel_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volPVBarrel");
      if (!vol)
        { fPVThickness = 0.; return false; }

      float min = ((TGeoTube*)vol->GetShape())->GetRmin();
      float max = ((TGeoTube*)vol->GetShape())->GetRmax();

      fTPCOuterX = std::abs(max - min);

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindTPCOuterX()
    { // We probably never used the name GArTPC_vol in a gdml file, but just in case
      TGeoVolume *vol_tpc = gGeoManager->FindVolumeFast("volGArTPC");
      if (!vol_tpc)
        vol_tpc = gGeoManager->FindVolumeFast("GArTPC_vol");
      if (!vol_tpc)
        { fTPCOuterX = 0.; return false; }

      fTPCOuterX = ((TGeoBBox*)vol_tpc->GetShape())->GetDZ();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALInnerSymmetry()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volBarrelECal");
      if (!vol)
        return false;

      fECALSymmetry = ((TGeoPgon*)vol->GetShape())->GetNedges();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALEndcapStartX()
    {
      //Find the ecal Endcap volume, which is centered at (0,0,0)
      TGeoVolume *vol_e = gGeoManager->FindVolumeFast("EndcapECal_vol");
      if (!vol_e)
        vol_e = gGeoManager->FindVolumeFast("volEndcapECal");
      if (!vol_e)
        return false;

      // There better be a stave 01 in module 00 in an endcap for this geometry!
      TGeoVolume *vol_endStave = gGeoManager->FindVolumeFast("EndcapECal_stave01_module00_vol");
      if (!vol_endStave)
        return false;

      fECALEndcapStartX  =   ((TGeoBBox*)vol_e->GetShape())->GetDZ();
      fECALEndcapStartX -= 2*((TGeoBBox*)vol_endStave->GetShape())->GetDZ();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALEndcapOuterX()
    {
      //Find the ecal Endcap
      TGeoVolume *vol_e = gGeoManager->FindVolumeFast("EndcapECal_vol");
      if (!vol_e)
        vol_e = gGeoManager->FindVolumeFast("volEndcapECal");
      if (!vol_e)
        return false;

      fECALEndcapOuterX = ((TGeoBBox*)vol_e->GetShape())->GetDZ();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindECALnLayers()
    {
      if (fECALSegmentationAlg)
        fECALnLayers = fECALSegmentationAlg->nLayers();

      return true;
    }

    //----------------------------------------------------------------------------
    unsigned int GeometryCore::GetNLayers(std::string det) const
    {
      if (det.compare("ECAL") == 0)
        return fECALnLayers;

      if (det.compare("TrackerSc") == 0)
        return fTrackerScnPlanes;

      if (det.compare("MuID") == 0)
        return fMuIDnLayers;

      return 0;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::MakeECALLayeredCalorimeterData()
    {
      //Barrel
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout] = std::shared_ptr<gar::geo::LayeredCalorimeterData>(new LayeredCalorimeterData());

      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->layoutType = gar::geo::LayeredCalorimeterData::BarrelLayout;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->inner_symmetry = 0;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->outer_symmetry = GetECALInnerSymmetry();
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->phi0 = 0;

      /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ].
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->extent[0] = GetECALInnerBarrelRadius() ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->extent[1] = GetECALOuterBarrelRadius() ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->extent[2] = 0.f ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->extent[3] = GetECALEndcapStartX() ;

      for(unsigned int i = 0; i < GetNLayers("ECAL"); i++) {
        float nRadiationLengths = 0.;
        float nInteractionLengths = 0.;
        float thickness_sum = 0.;
        float layer_thickness = 0.;
        float abs_thickness = 0.;

        gar::geo::LayeredCalorimeterData::Layer caloLayer;
        // caloLayer.cellSize0 = fECALSegmentationAlg->gridSizeX();
        // caloLayer.cellSize1 = fECALSegmentationAlg->gridSizeX();
        caloLayer.cellSize0 = fECALSegmentationAlg->stripSizeX();
        caloLayer.cellSize1 = fECALSegmentationAlg->stripSizeX();

        TGeoVolume *volLayer = gGeoManager->FindVolumeFast(TString::Format("BarrelECal_stave01_module01_layer_%02i_vol", i+1));

        if (volLayer) {
          for(int islice = 0; islice < volLayer->GetNdaughters(); islice++) {

            TGeoVolume *slice = volLayer->GetNode(islice)->GetVolume();
            TGeoMaterial *mat = slice->GetMaterial();
            double rad_len = mat->GetRadLen();
            double int_len = mat->GetIntLen();
            const char *material = mat->GetName();

            double slice_thickness = ((TGeoBBox*)slice->GetShape())->GetDZ(); //half thickness

            nRadiationLengths   += slice_thickness/rad_len;
            nInteractionLengths += slice_thickness/int_len;
            thickness_sum       += slice_thickness;
            layer_thickness     += slice_thickness;

            MF_LOG_DEBUG("GeometryCore::FindECALLayeredCalorimeterData")
              << " Slice " << islice
              << " RadLen " << nRadiationLengths
              << " IntLen " << nInteractionLengths
              << " ThicknessSum " << thickness_sum
              << " Material " << mat->GetName();

            if(strcmp(material, "Copper") == 0 || strcmp(material, "Steel") == 0 || strcmp(material, "Lead") == 0) {
              abs_thickness += slice_thickness * 2;
            }

            if(strcmp(material, "Scintillator") == 0) {
              caloLayer.inner_nRadiationLengths = nRadiationLengths;
              caloLayer.inner_nInteractionLengths = nInteractionLengths;
              caloLayer.inner_thickness = thickness_sum;
              caloLayer.sensitive_thickness = slice_thickness*2.;

              nRadiationLengths = 0.;
              nInteractionLengths = 0.;
              thickness_sum = 0.;
            }

            nRadiationLengths   += slice_thickness/rad_len;
            nInteractionLengths += slice_thickness/int_len;
            thickness_sum       += slice_thickness;
            layer_thickness     += slice_thickness;
          }

          caloLayer.outer_nRadiationLengths = nRadiationLengths;
          caloLayer.outer_nInteractionLengths = nInteractionLengths;
          caloLayer.outer_thickness = thickness_sum;

          caloLayer.distance = GetECALInnerBarrelRadius() + (i+1) * layer_thickness;
          caloLayer.absorberThickness = abs_thickness;
          fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::BarrelLayout]->layers.push_back( caloLayer );
        }
      }

      //Endcap
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout] = std::shared_ptr<gar::geo::LayeredCalorimeterData>(new LayeredCalorimeterData());

      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->layoutType = gar::geo::LayeredCalorimeterData::EndcapLayout;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->inner_symmetry = 0;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->outer_symmetry = GetECALInnerSymmetry();
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->phi0 = 0;

      /// extent of the calorimeter in the r-z-plane [ rmin, rmax, zmin, zmax ].
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->extent[0] = GetECALInnerEndcapRadius() ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->extent[1] = GetECALOuterEndcapRadius() ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->extent[2] = GetECALEndcapStartX() ;
      fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->extent[3] = GetECALEndcapOuterX() ;

      for(unsigned int i = 0; i < GetNLayers("ECAL"); i++) {

        float nRadiationLengths = 0.;
        float nInteractionLengths = 0.;
        float thickness_sum = 0.;
        float layer_thickness = 0.;
        float abs_thickness = 0.;

        gar::geo::LayeredCalorimeterData::Layer caloLayer;
        // caloLayer.cellSize0 = fECALSegmentationAlg->gridSizeX();
        // caloLayer.cellSize1 = fECALSegmentationAlg->gridSizeX();
        caloLayer.cellSize0 = fECALSegmentationAlg->stripSizeX();
        caloLayer.cellSize1 = fECALSegmentationAlg->stripSizeX();

        TGeoVolume *volLayer = gGeoManager->FindVolumeFast(TString::Format("EndcapECal_stave01_module00_layer_%02i_vol", i+1));

        if (volLayer) {
		
          for(int islice = 0; islice < volLayer->GetNdaughters(); islice++) {

            TGeoVolume *slice = volLayer->GetNode(islice)->GetVolume();
            TGeoMaterial *mat = slice->GetMaterial();
            double rad_len = mat->GetRadLen();
            double int_len = mat->GetIntLen();
            const char *material = mat->GetName();

            double slice_thickness = ((TGeoBBox*)slice->GetShape())->GetDZ(); //half thickness

            nRadiationLengths   += slice_thickness/rad_len;
            nInteractionLengths += slice_thickness/int_len;
            thickness_sum       += slice_thickness;
            layer_thickness     += slice_thickness;

            MF_LOG_DEBUG("GeometryCore::FindECALLayeredCalorimeterData")
              << " Slice " << islice
              << " RadLen " << nRadiationLengths
              << " IntLen " << nInteractionLengths
              << " ThicknessSum " << thickness_sum
              << " Material " << mat->GetName();

            if (strcmp(material, "Copper") == 0 || strcmp(material, "Steel") == 0 || strcmp(material, "Lead") == 0) {
              abs_thickness += slice_thickness * 2;
            }

            if(strcmp(material, "Scintillator") == 0) {
              caloLayer.inner_nRadiationLengths = nRadiationLengths;
              caloLayer.inner_nInteractionLengths = nInteractionLengths;
              caloLayer.inner_thickness = thickness_sum;
              caloLayer.sensitive_thickness = slice_thickness*2.;

              nRadiationLengths = 0.;
              nInteractionLengths = 0.;
              thickness_sum = 0.;
            }

            nRadiationLengths   += slice_thickness/rad_len;
            nInteractionLengths += slice_thickness/int_len;
            thickness_sum       += slice_thickness;
            layer_thickness     += slice_thickness;
          }

          caloLayer.outer_nRadiationLengths = nRadiationLengths;
          caloLayer.outer_nInteractionLengths = nInteractionLengths;
          caloLayer.outer_thickness = thickness_sum;

          caloLayer.distance = GetECALEndcapStartX() + i * layer_thickness;
          caloLayer.absorberThickness = abs_thickness;
          fECALLayeredCalorimeterData[gar::geo::LayeredCalorimeterData::EndcapLayout]->layers.push_back( caloLayer );
        }
      }

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindMuIDInnerBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("YokeBarrel_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volYokeBarrel");
      if (!vol)
        return false;

      fMuIDRinner = ((TGeoPgon*)vol->GetShape())->GetRmin(0);

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindMuIDOuterBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("YokeBarrel_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volYokeBarrel");
      if (!vol)
        return false;

      fMuIDRouter = ((TGeoPgon*)vol->GetShape())->GetRmax(0);

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindMuIDInnerSymmetry()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("YokeBarrel_vol");
      if (!vol)
        vol = gGeoManager->FindVolumeFast("volYokeBarrel");
      if (!vol)
        return false;

      fMuIDSymmetry = ((TGeoPgon*)vol->GetShape())->GetNedges();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindMuIDnLayers()
    {
      if (fMuIDSegmentationAlg)
        fMuIDnLayers = fMuIDSegmentationAlg->nLayers();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindTrackerScVolume()
    {
      std::vector<TGeoNode const*> path = this->FindVolumePath("volTracker");
      if (path.size() == 0)
        return false;

      const TGeoNode *Tracker_node = path.at(path.size()-1);
      if (Tracker_node == nullptr) {
        std::cout << "Cannot find node volTracker_0" << std::endl;
        return false;
      }

      TGeoMatrix *mat = Tracker_node->GetMatrix();
      const double *origin = mat->GetTranslation();

      //Get the origin correctly
      fGArLiteXCent =  fWorldX + fRockX + fEnclosureX + fMPDX + origin[0];
      fGArLiteYCent =  fWorldY + fRockY + fEnclosureY + fMPDY + origin[1];
      fGArLiteZCent =  fWorldZ + fRockZ + fEnclosureZ + fMPDZ + origin[2];

      //Get the dimension of the active volume
      TGeoVolume *activeVol = Tracker_node->GetVolume();
      float rOuter = ((TGeoTube*)activeVol->GetShape())->GetRmax();

      fGArLiteRadius = rOuter;
      fGArLiteLength = 2.0 * ((TGeoTube*)activeVol->GetShape())->GetDZ();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::FindTrackerScnPlanes()
    {
      if (fMinervaSegmentationAlg)
        fTrackerScnPlanes = fMinervaSegmentationAlg->nPlanes();

      return true;
    }

    //----------------------------------------------------------------------------
    void GeometryCore::InitVariables()
    {
      fEnclosureX = 0.;
      fEnclosureY = 0.;
      fEnclosureZ = 0.;

      fTPCNumDriftVols = 2;  // default to ALICE-style TPC with cathode in the middle

      fTPCXCent = 0.;
      fTPCYCent = 0.;
      fTPCZCent = 0.;

      fMPDX = 0.;
      fMPDY = 0.;
      fMPDZ = 0.;

      fLArTPCX = 0.;
      fLArTPCY = 0.;
      fLArTPCZ = 0.;

      fLArTPCXCent = 0.;
      fLArTPCYCent = 0.;
      fLArTPCZCent = 0.;

      fTPCRadius = 9999.;
      fTPCLength = 9999.;

      fEnclosureHalfWidth = 9999.;
      fEnclosureHalfHeight = 9999.;
      fEnclosureLength = 9999.;

      fMPDHalfWidth = 9999.;
      fMPDHalfHeight = 9999.;
      fMPDLength = 9999.;

      fLArTPCHalfWidth = 0.;
      fLArTPCHalfHeight = 0.;
      fLArTPCLength = 0.;

      fLArTPCActiveHalfWidth = 0.;
      fLArTPCActiveHalfHeight = 0.;
      fLArTPCActiveLength = 0.;

      fECALRinner = 0.;
      fECALRouter = 0.;
      fECALECapRouter = 0.;
      fECALECapRouter = 0.;
      fPVThickness = 0.;
      fTPCOuterX = 0.;
      fECALSymmetry = -1;
      fECALEndcapStartX = 0.;
      fECALnLayers = 0;
      fECALNodePath.clear();

      fHasRock = false;
      fHasEnclosure = false;
      fHasLArTPCDetector = false;
      fHasGasTPCDetector = false;
      fHasECALDetector = false;
      fHasTrackerScDetector = false;
      fHasMuonDetector = false;

      fMuIDRinner = 0.;
      fMuIDRouter = 0.;
      fMuIDSymmetry = -1;
      fMuIDnLayers = 0;

      fGArLiteXCent = 0.;
      fGArLiteYCent = 0.;
      fGArLiteZCent = 0.;
      fTrackerScnPlanes = 0;

      fGArLiteRadius = 9999.;
      fGArLiteLength = 9999.;
    }

    //--------------------------------------------------------------------
    constexpr details::geometry_iterator_types::BeginPos_t
    details::geometry_iterator_types::begin_pos;
    constexpr details::geometry_iterator_types::EndPos_t
    details::geometry_iterator_types::end_pos;
    constexpr details::geometry_iterator_types::UndefinedPos_t
    details::geometry_iterator_types::undefined_pos;

    //--------------------------------------------------------------------
    //--- ROOTGeoNodeForwardIterator
    //---
    ROOTGeoNodeForwardIterator& ROOTGeoNodeForwardIterator::operator++ () {
      if (current_path.empty()) return *this;
      if (current_path.size() == 1) { current_path.pop_back(); return *this; }

      // I am done; all my descendants were also done already;
      // first look at my younger siblings
      NodeInfo_t& current = current_path.back();
      NodeInfo_t const& parent = current_path[current_path.size() - 2];
      if (++(current.sibling) < parent.self->GetNdaughters()) {
        // my next sibling exists, let's parse his descendents
        current.self = parent.self->GetDaughter(current.sibling);
        reach_deepest_descendant();
      } else
        current_path.pop_back(); // no sibling, it's time for mum
      return *this;
    } // ROOTGeoNodeForwardIterator::operator++


    //--------------------------------------------------------------------
    std::vector<TGeoNode const*> ROOTGeoNodeForwardIterator::get_path() const {

      std::vector<TGeoNode const*> node_path(current_path.size());

      std::transform(current_path.begin(), current_path.end(), node_path.begin(),
             [](NodeInfo_t const& node_info){ return node_info.self; });
      return node_path;

    } // ROOTGeoNodeForwardIterator::path()


    //--------------------------------------------------------------------
    void ROOTGeoNodeForwardIterator::reach_deepest_descendant() {
      Node_t descendent = current_path.back().self;
      while (descendent->GetNdaughters() > 0) {
        descendent = descendent->GetDaughter(0);
        current_path.emplace_back(descendent, 0U);
      } // while
    } // ROOTGeoNodeForwardIterator::reach_deepest_descendant()

    //--------------------------------------------------------------------
    void ROOTGeoNodeForwardIterator::init(TGeoNode const* start_node) {
      current_path.clear();
      if (!start_node) return;
      current_path.emplace_back(start_node, 0U);
      reach_deepest_descendant();
    } // ROOTGeoNodeForwardIterator::init()

  } // namespace geo
} // gar
