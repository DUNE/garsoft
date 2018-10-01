/**
 * @file   GeometryCore.cxx
 * @brief  Access the description of detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    GeometryCore.h
 *
 */

// class header
#include "Geometry/GeometryCore.h"
#include "Geometry/ChannelMapAlg.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"
#include "TGeoTube.h"

// C/C++ includes
#include <cstddef> // size_t
#include <cctype> // ::tolower()
#include <cmath> // std::abs() ...
#include <vector>
#include <algorithm> // std::for_each(), std::transform()
#include <utility> // std::swap()
#include <limits> // std::numeric_limits<>
#include <memory> // std::default_deleter<>

namespace gar {
  namespace geo {

    template <typename T>
    inline T sqr(T v) { return v * v; }


    //......................................................................
    // Constructor.
    GeometryCore::GeometryCore(fhicl::ParameterSet const& pset)
    : fSurfaceY         (pset.get< double            >("SurfaceY"               ))
    , fDetectorName     (pset.get< std::string       >("Name"                   ))
    , fPositionWiggle   (pset.get< double            >("PositionEpsilon",  1.e-4))
    {
      std::transform(fDetectorName.begin(), fDetectorName.end(),
                     fDetectorName.begin(), ::tolower);
    } // GeometryCore::GeometryCore()


    //......................................................................
    GeometryCore::~GeometryCore()
    {
      ClearGeometry();
    } // GeometryCore::~GeometryCore()


    //......................................................................
    void GeometryCore::ApplyChannelMap(std::shared_ptr<geo::ChannelMapAlg> pChannelMap)
    {
      pChannelMap->Initialize(*this);
      fChannelMapAlg = pChannelMap;
    } // GeometryCore::ApplyChannelMap()

    //......................................................................
    void GeometryCore::LoadGeometryFile(std::string const& gdmlfile,
                                        std::string const& rootfile,
                                        bool bForceReload /* = false*/)
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
        TGeoManager::Import(rootfile.c_str());
        gGeoManager->LockGeometry();
      }

      fGDMLfile = gdmlfile;
      fROOTfile = rootfile;

      LOG_INFO("GeometryCore")
      << "New detector geometry loaded from "
      << "\n\t" << fROOTfile
      << "\n\t" << fGDMLfile;

      this->SetEnclosureVolumePosition();
      this->FindActiveTPCVolume();

      //Load in memory ECAL geometry
      this->SetECALLayerThickness();
      this->SetECALEndcapStartPosition();
      this->SetECALInnerBarrelRadius();
      this->SetECALOuterBarrelRadius();

      this->PrintGeometry();

    } // GeometryCore::LoadGeometryFile()

    //......................................................................
    void GeometryCore::SetEnclosureVolumePosition()
    {
      TGeoManager *geo = ROOTGeoManager();
      TGeoNode *world_node = geo->GetTopVolume()->FindNode("volDetEnclosure_0");
      const double *origin = world_node->GetMatrix()->GetTranslation();

      fEnclosureX = origin[0];
      fEnclosureY = origin[1];
      fEnclosureZ = origin[2];

      return;
    }

    //......................................................................
    void GeometryCore::FindActiveTPCVolume()
    {
      // check if the current level of the detector is the active TPC volume, if
      // not, then dig a bit deeper

      TGeoManager *geo = ROOTGeoManager();
      TGeoVolume *det_vol = geo->FindVolumeFast("volDetEnclosure");
      TGeoNode *GAr_node = det_vol->FindNode("volNDHPgTPC_0");
      TGeoMatrix *mat = GAr_node->GetMatrix();
      const double *origin = mat->GetTranslation();

      //Get the origin correctly
      fTPCXCent =            fEnclosureX + origin[0];
      fTPCYCent =            fEnclosureY + origin[1];
      fTPCZCent =            fEnclosureZ + origin[2];

      //Get the dimension of the active volume
      TGeoVolume *activeVol = geo->FindVolumeFast("volGArTPC");
      fDetHalfWidth  =       ((TGeoBBox*)activeVol->GetShape())->GetDZ();
      fDetHalfHeight =       ((TGeoBBox*)activeVol->GetShape())->GetDY();
      fDetLength     = 2.0 * ((TGeoBBox*)activeVol->GetShape())->GetDX();

      return;
    }

    //......................................................................
    void GeometryCore::ClearGeometry() {

      return;
    } // GeometryCore::ClearGeometry()


    //......................................................................
    TGeoManager* GeometryCore::ROOTGeoManager() const
    {
      return gGeoManager;
    }

    //......................................................................
    unsigned int GeometryCore::NChannels() const
    {
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
    const std::string GeometryCore::GetWorldVolumeName() const
    {
        // For now, and possibly forever, this is a constant (given the
        // definition of "nodeNames" above).

      return std::string(ROOTGeoManager()->GetTopNode()->GetName());
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
      if(!s)
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
      TGeoVolume *volWorld = gGeoManager->FindVolumeFast(this->GetWorldVolumeName().c_str());
      float halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
      float halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
      float halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
      if(std::abs(point.x()) > halfwidth  ||
         std::abs(point.y()) > halfheight ||
         std::abs(point.z()) > halflength){
        LOG_WARNING("GeometryCoreBadInputPoint")
        << "point ("
        << point.x() << ","
        << point.y() << ","
        << point.z() << ") "
        << "is not inside the world volume "
        << " half width = "  << halfwidth
        << " half height = " << halfheight
        << " half length = " << halflength;
        return false;
      }

      return true;
    }

    //......................................................................
    const std::string GeometryCore::VolumeName(TVector3 const& point) const
    {
      if( !this->PointInWorld(point) ){
        const std::string unknown("unknownVolume");
        return unknown;
      }

      const std::string name(gGeoManager->FindNode(point.x(), point.y(), point.z())->GetName());
      return name;
    }

    //......................................................................
    const std::string GeometryCore::MaterialName(TVector3 const& point)
    {
      if( !this->PointInWorld(point) ){
        const std::string unknown("unknownVolume");
        return unknown;
      }

      const std::string name(gGeoManager->FindNode(point.x(),
                                                   point.y(),
                                                   point.z())->GetMedium()->GetMaterial()->GetName());
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
    unsigned int GeometryCore::NearestChannel(const float worldLoc[3]) const
    {
      return fChannelMapAlg->NearestChannel(worldLoc);
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::NearestChannel(std::vector<float> const& worldLoc) const
    {
      float loc[3] = {worldLoc[0],
                      worldLoc[1],
                      worldLoc[2]};

      return this->NearestChannel(loc);
    }

    //--------------------------------------------------------------------
    unsigned int GeometryCore::NearestChannel(TVector3 const& worldLoc) const
    {
      float loc[3] = {(float)worldLoc[0],
                      (float)worldLoc[1],
                      (float)worldLoc[2]};

      return this->NearestChannel(loc);
    }

    //--------------------------------------------------------------------
    void GeometryCore::ChannelToPosition(unsigned int const channel,
                                         float*       const worldLoc) const
    {
      return fChannelMapAlg->ChannelToPosition(channel, worldLoc);
    }

    //--------------------------------------------------------------------
    bool GeometryCore::SetECALLayerThickness()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("IECLayer_vol");
      if(!vol) return false;
      fECALLayerThickness = ((TGeoBBox*)vol->GetShape())->GetDZ() * 2;

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::SetECALInnerBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("InnerBarrelECal_vol");
      if(!vol) return false;
      fECALRinner = ((TGeoTube*)vol->GetShape())->GetRmin();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::SetECALOuterBarrelRadius()
    {
      TGeoVolume *vol = gGeoManager->FindVolumeFast("OuterBarrelECal_vol");
      if(!vol) return false;
      fECALRouter = ((TGeoTube*)vol->GetShape())->GetRmin();

      return true;
    }

    //----------------------------------------------------------------------------
    bool GeometryCore::SetECALEndcapStartPosition()
    {
      TGeoVolume *det_vol = gGeoManager->FindVolumeFast("volNDHPgTPC");
      TGeoVolume *endcap_vol = gGeoManager->FindVolumeFast("InnerEndcapECal_vol");

      double thickness = 0.;

      if(!endcap_vol) return false;
      thickness = ((TGeoTube*)endcap_vol->GetShape())->GetDZ();

      if(!det_vol) return false;

      TGeoNode *endcap_node = det_vol->FindNode("InnerEndcapECal_vol_0");
      TGeoMatrix *mat = endcap_node->GetMatrix();
      const double *origin = mat->GetTranslation();
      fEndcapStartXPosition = std::abs(origin[0]) - thickness;

      return true;
    }

    //----------------------------------------------------------------------------
    void GeometryCore::PrintGeometry()
    {
      //Prints geometry parameters
      std::cout << "------------------------------" << std::endl;
      std::cout << "Enclosure Geometry" << std::endl;
      std::cout << "Enclosure Origin (x, y, z) " << GetEnclosureX() << " cm " << GetEnclosureY() << " cm " << GetEnclosureZ() << " cm" << std::endl;
      std::cout << "------------------------------" << std::endl;
      std::cout << "TPC Geometry" << std::endl;
      std::cout << "TPC Origin (x, y, z) " << TPCXCent() << " cm " << TPCYCent() << " cm " << TPCZCent() << " cm" << std::endl;
      std::cout << "TPC Active Volume Size (H, W, L) " << DetHalfHeight() << " cm " << DetHalfWidth() << " cm " << DetLength() << " cm" << std::endl;
      std::cout << "------------------------------" << std::endl;
      std::cout << "ECAL Geometry" << std::endl;
      std::cout << "Layer thickness: " << GetECALLayerThickness() << " cm" << std::endl;
      std::cout << "ECAL Endcap front face position in x: " << GetECALEndcapStartPosition() << " cm" << std::endl;
      std::cout << "ECAL Inner Barrel minimum radius: " << GetECALInnerBarrelRadius() << " cm" << std::endl;
      std::cout << "ECAL Outer Barrel minimum radius: " << GetECALOuterBarrelRadius() << " cm" << std::endl;
      std::cout << "------------------------------" << std::endl;
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
      }
      else current_path.pop_back(); // no sibling, it's time for mum
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
