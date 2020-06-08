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
#include "Geometry/ChannelMapAlgs/ECALSegmentationAlg.h"

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
        // Constructor.
        GeometryCore::GeometryCore(fhicl::ParameterSet const& pset)
        : fSurfaceY         (pset.get< double            >("SurfaceY"               ))
        , fDetectorName     (pset.get< std::string       >("Name"                   ))
        , fPositionWiggle   (pset.get< double            >("PositionEpsilon",  1.e-4))
        , fPointInWarnings  (pset.get< bool              >("PointInWarnings",  false))
        {
            std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);

            InitVariables();
        } // GeometryCore::GeometryCore()


        //......................................................................
        GeometryCore::~GeometryCore()
        {
            ClearGeometry();
        } // GeometryCore::~GeometryCore()


        //......................................................................
        void GeometryCore::ApplyChannelMap(std::shared_ptr<geo::seg::ChannelMapAlg> pChannelMap)
        {
            FindWorldVolume();
            FindRockVolume();
            FindEnclosureVolume();
            FindMPDVolume();
            FindActiveTPCVolume();

            //For the LArTPC
            FindLArTPCVolume();
            FindActiveLArTPCVolume();

            pChannelMap->Initialize(*this);
            fChannelMapAlg = pChannelMap;
        } // GeometryCore::ApplyChannelMap()

        //......................................................................
        void GeometryCore::ApplyECALSegmentationAlg(std::shared_ptr<geo::seg::ECALSegmentationAlg> pECALSegmentationAlg)
        {
            pECALSegmentationAlg->Initialize(*this);
            fECALSegmentationAlg = pECALSegmentationAlg;

            StoreECALParameters();
            PrintGeometry();
        } // GeometryCore::ApplyECALSegmentationAlg()

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
                TGeoManager::Import(rootfile.c_str());
                gGeoManager->LockGeometry();
            }

            fGDMLfile = gdmlfile;
            fROOTfile = rootfile;

            LOG_INFO("GeometryCore")
            << "New detector geometry loaded from "
            << "\n\t" << fROOTfile
            << "\n\t" << fGDMLfile;

        } // GeometryCore::LoadGeometryFile()

        //......................................................................
        void GeometryCore::FindWorldVolume()
        {
            std::vector<TGeoNode const*> path = this->FindVolumePath("volWorld");
            if(path.size() == 0) return;

            const TGeoNode *world_node = path.at(path.size()-1);
            if(world_node == nullptr) {
                std::cout << "Cannot find node volWorld_0" << std::endl;
                return;
            }

            const double *origin = world_node->GetMatrix()->GetTranslation();

            fWorldX = origin[0];
            fWorldY = origin[1];
            fWorldZ = origin[2];

            fWorldHalfWidth = ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDZ();
            fWorldHalfHeight = ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDY();
            fWorldLength = 2.0 * ((TGeoBBox*)world_node->GetVolume()->GetShape())->GetDX();

            return;
        }

        //......................................................................
        void GeometryCore::FindRockVolume()
        {
            std::vector<TGeoNode const*> path = this->FindVolumePath("rockBox_lv");
            if(path.size() == 0) return;

            const TGeoNode *rock_node = path.at(path.size()-1);
            if(rock_node == nullptr) {
                std::cout << "Cannot find node rockBox_lv_0" << std::endl;
                return;
            }

            const double *origin = rock_node->GetMatrix()->GetTranslation();

            fRockX = origin[0];
            fRockY = origin[1];
            fRockZ = origin[2];

            fRockHalfWidth = ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDZ();
            fRockHalfHeight = ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDY();
            fRockLength = 2.0 * ((TGeoBBox*)rock_node->GetVolume()->GetShape())->GetDX();

            return;
        }

        //......................................................................
        void GeometryCore::FindEnclosureVolume()
        {
            std::vector<TGeoNode const*> path = this->FindVolumePath("volDetEnclosure");
            if(path.size() == 0) return;

            const TGeoNode *enc_node = path.at(path.size()-1);
            if(enc_node == nullptr) {
                std::cout << "Cannot find node volDetEnclosure_0" << std::endl;
                return;
            }

            const double *origin = enc_node->GetMatrix()->GetTranslation();

            fEnclosureX = origin[0];
            fEnclosureY = origin[1];
            fEnclosureZ = origin[2];

            fEnclosureHalfWidth = ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDZ();
            fEnclosureHalfHeight = ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDY();
            fEnclosureLength = 2.0 * ((TGeoBBox*)enc_node->GetVolume()->GetShape())->GetDX();

            return;
        }

        //......................................................................
        void GeometryCore::FindLArTPCVolume()
        {
            std::vector<TGeoNode const*> path = this->FindVolumePath("volArgonCubeDetector");
            if(path.size() == 0)
                return;

            const TGeoNode *lar_node = path.at(path.size()-1);
            if(lar_node == nullptr) {
                std::cout << "Cannot find node volArgonCubeDetector_0" << std::endl;
                return;
            }

            const double *origin = lar_node->GetMatrix()->GetTranslation();

            fLArTPCX = origin[0];
            fLArTPCY = origin[1];
            fLArTPCZ = origin[2];

            fLArTPCHalfWidth = ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDZ();
            fLArTPCHalfHeight = ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDY();
            fLArTPCLength = 2.0 * ((TGeoBBox*)lar_node->GetVolume()->GetShape())->GetDX();

            return;
        }

        //......................................................................
        void GeometryCore::FindMPDVolume()
        {
            std::vector<TGeoNode const*> path = this->FindVolumePath("volMPD");
            if(path.size() == 0)
                return;

            const TGeoNode *mpd_node = path.at(path.size()-1);
            if(mpd_node == nullptr) {
                std::cout << "Cannot find node volMPD_0" << std::endl;
                return;
            }

            const double *origin = mpd_node->GetMatrix()->GetTranslation();

            fMPDX = origin[0];
            fMPDY = origin[1];
            fMPDZ = origin[2];

            fMPDHalfWidth = ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDZ();
            fMPDHalfHeight = ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDY();
            fMPDLength = 2.0 * ((TGeoBBox*)mpd_node->GetVolume()->GetShape())->GetDX();

            return;
        }


        //......................................................................
        void GeometryCore::FindActiveLArTPCVolume()
        {
            // check if the current level of the detector is the active TPC volume, if
            // not, then dig a bit deeper
            std::vector<TGeoNode const*> path = this->FindVolumePath("volArgonCubeActive");
            if(path.size() == 0)
                return;

            const TGeoNode *LArTPC_node = path.at(path.size()-1);
            if(LArTPC_node == nullptr) {
                std::cout << "Cannot find node volArgonCubeActive_0" << std::endl;
                return;
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

            return;
        }

        //......................................................................
        void GeometryCore::FindActiveTPCVolume()
        {
            // check if the current level of the detector is the active TPC volume, if
            // not, then dig a bit deeper
            std::vector<TGeoNode const*> path = this->FindVolumePath("TPCGas_vol");
            if(path.size() == 0)
                path = this->FindVolumePath("volTPCGas");
            if(path.size() == 0)
                return;

            const TGeoNode *GArTPC_node = path.at(path.size()-1);
            if(GArTPC_node == nullptr) {
                std::cout << "Cannot find node TPCGas_vol_0/TPCGas_vol_0" << std::endl;
                return;
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

            return;
        }

        //......................................................................
        const float GeometryCore::GetSensVolumeThickness(const TVector3& point) const
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
            //Check if it is ECAL endcap -> layer size is not the BBox! It is the apothem
            std::string volname = vol->GetName();
            bool isBarrel = true;
            if(volname.find("endcap") != std::string::npos || volname.find("Endcap") != std::string::npos ) isBarrel = false;

            std::array<double, 3> shape;

            if(vol)
            {
                TGeoBBox *box = (TGeoBBox*)(vol->GetShape());

                if(isBarrel) {
                    shape[0] = box->GetDX();
                    shape[1] = box->GetDY();
                } else {
                    shape[0] = GetECALEndcapApothemLength() / 2.;
                    shape[1] = GetECALEndcapApothemLength() / 2.;
                }

                shape[2] = box->GetDZ();

                return shape; //return half size in cm
            }
            else{
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
            TGeoVolume *volWorld = gGeoManager->FindVolumeFast("volWorld");
            float halflength = ((TGeoBBox*)volWorld->GetShape())->GetDZ();
            float halfheight = ((TGeoBBox*)volWorld->GetShape())->GetDY();
            float halfwidth  = ((TGeoBBox*)volWorld->GetShape())->GetDX();
            if (std::abs(point.x()) > halfwidth  ||
                std::abs(point.y()) > halfheight ||
                std::abs(point.z()) > halflength){
                if (fPointInWarnings) {
                    LOG_WARNING("GeometryCoreBadInputPoint")
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
                    LOG_WARNING("GeometryCoreBadInputPoint")
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
            TVector3 tpc_origin(TPCXCent(), TPCYCent(), TPCZCent());
            TVector3 new_point = point - tpc_origin;
            // check that the given point is in the enclosure volume at least
            if (std::abs(new_point.x()) > fMPDHalfWidth  ||
                std::abs(new_point.y()) > fMPDHalfHeight ||
                std::abs(new_point.z()) > fMPDLength){
                if (fPointInWarnings) {
                    LOG_WARNING("GeometryCoreBadInputPoint")
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
            float y = std::abs(point.y() - fTPCYCent);
            float z = std::abs(point.z() - fTPCZCent);
            if (std::abs(point.x() - fTPCXCent) > fTPCLength/2.0  ||
                std::hypot(z,y)                 > fTPCRadius) {
                if (fPointInWarnings) {
                    LOG_WARNING("GeometryCoreBadInputPoint")
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
                    LOG_WARNING("GeometryCoreBadInputPoint")
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

            if( vol_name.find("barrel") == std::string::npos || vol_name.find("Barrel") == std::string::npos)
                return false;

            return true;
        }

        //......................................................................
        bool GeometryCore::PointInECALEndcap(TVector3 const& point) const
        {
            std::string vol_name = this->VolumeName(point);

            if( vol_name.find("endcap") == std::string::npos || vol_name.find("Endcap") == std::string::npos)
                return false;

            return true;
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
        void GeometryCore::StoreECALParameters()
        {
            this->FindECALInnerBarrelRadius();
            this->FindECALOuterBarrelRadius();
            this->FindECALInnerEndcapRadius();
            this->FindECALOuterEndcapRadius();
            this->FindPVThickness();
            this->FindECALInnerSymmetry();
            this->FindECALEndcapStartX();
            this->FindECALEndcapOuterX();
                        
            StoreECALNodes(fECALNodePath);
            std::cout << "Stored " << fECALNodePath.size() << " ECAL nodes in memory" << std::endl;

            // std::raise(SIGINT);
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALInnerBarrelRadius()
        {
            TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
            if(!vol)
                vol = gGeoManager->FindVolumeFast("volBarrelECal");
            if(!vol)
                return false;

            fECALRinner = ((TGeoPgon*)vol->GetShape())->GetRmin(0);

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALOuterBarrelRadius()
        {
            TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
            if(!vol)
                vol = gGeoManager->FindVolumeFast("volBarrelECal");
            if(!vol)
                return false;

            fECALRouter = ((TGeoPgon*)vol->GetShape())->GetRmax(0);

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALInnerEndcapRadius()
        {
            TGeoVolume *vol = gGeoManager->FindVolumeFast("EndcapECal_vol");
            if(!vol)
                vol = gGeoManager->FindVolumeFast("volEndcapECal");
            if(!vol)
                return false;

            fECALECapRinner = 0.;

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALOuterEndcapRadius()
        {
            // TGeoVolume *vol = gGeoManager->FindVolumeFast("EndcapECal_vol");
            // if(!vol)
            // vol = gGeoManager->FindVolumeFast("volEndcapECal");
            // if(!vol)
            // return false;

            // fECALECapRouter = ((TGeoBBox*)vol->GetShape())->GetDX();//Not great as the shape is not a tube...
            fECALECapRouter = fECALRouter; //should be equal

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindPVThickness()
        {
            TGeoVolume *vol = gGeoManager->FindVolumeFast("PVBarrel_vol");
            if(!vol)
                vol = gGeoManager->FindVolumeFast("volPVBarrel");
            if(!vol)
                { fPVThickness = 0.; return false; }

            float min = ((TGeoTube*)vol->GetShape())->GetRmin();
            float max = ((TGeoTube*)vol->GetShape())->GetRmax();

            fPVThickness = std::abs(max - min);

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALInnerSymmetry()
        {
            TGeoVolume *vol = gGeoManager->FindVolumeFast("BarrelECal_vol");
            if(!vol)
                vol = gGeoManager->FindVolumeFast("volBarrelECal");
            if(!vol)
                return false;

            fECALSymmetry = ((TGeoPgon*)vol->GetShape())->GetNedges();

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALEndcapStartX()
        {
            //Find the PV Endcap
            TGeoVolume *vol_pv = gGeoManager->FindVolumeFast("PVEndcap_vol");
            if(!vol_pv)
                vol_pv = gGeoManager->FindVolumeFast("volPVEndcap");
            if(!vol_pv)
                return false;

            TGeoVolume *vol_tpc_chamber = gGeoManager->FindVolumeFast("volGArTPC");
            if(!vol_tpc_chamber) return false;

            //The start of the endcap is after the pv endcap -> sum of tpc chamber length and pressure vessel bulge
            fECALEndcapStartX = ((TGeoBBox*)vol_pv->GetShape())->GetDZ()*2 + ((TGeoBBox*)vol_tpc_chamber->GetShape())->GetDZ();

            return true;
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::FindECALEndcapOuterX()
        {
            //Find the ecal Endcap
            TGeoVolume *vol_e = gGeoManager->FindVolumeFast("EndcapECal_vol");
            if(!vol_e)
                vol_e = gGeoManager->FindVolumeFast("volEndcapECal");
            if(!vol_e)
                return false;

            //The start of the endcap is after the pv endcap -> sum of tpc chamber length and pressure vessel bulge
            fECALEndcapOuterX = ((TGeoBBox*)vol_e->GetShape())->GetDZ()*2 - fECALEndcapStartX;

            return true;
        }

        //----------------------------------------------------------------------------
        gar::raw::CellID_t GeometryCore::GetCellID(const TGeoNode *node, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const
        {
            const std::array<double, 3> shape = this->FindShapeSize(node);

            fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);

            return fECALSegmentationAlg->GetCellID(*this, det_id, stave, module, layer, slice, localPosition);
        }

        //----------------------------------------------------------------------------
        std::string GeometryCore::GetCellIDEncoding() const
        {
            return fECALSegmentationAlg->cellEncoding();
        }

        //----------------------------------------------------------------------------
        std::array<double, 3> GeometryCore::GetPosition(const TGeoNode *node, const gar::raw::CellID_t &cID) const
        {
            const std::array<double, 3> shape = this->FindShapeSize(node);

            fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);

            return fECALSegmentationAlg->GetPosition(*this, cID);
        }

        //----------------------------------------------------------------------------
        bool GeometryCore::isTile(const gar::raw::CellID_t& cID) const
        {
            return fECALSegmentationAlg->isTile(cID);
        }

        //----------------------------------------------------------------------------
        double GeometryCore::getStripWidth() const { return fECALSegmentationAlg->stripSizeX(); }

        //----------------------------------------------------------------------------
        double GeometryCore::getTileSize() const { return fECALSegmentationAlg->gridSizeX(); }

        //----------------------------------------------------------------------------
        double GeometryCore::getStripLength(std::array<double, 3> const& point, const gar::raw::CellID_t &cID) const
        {
            std::array<double, 3> localtemp;
            gar::geo::LocalTransformation<TGeoHMatrix> trans;
            this->WorldToLocal(point, localtemp, trans);

            const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

            fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);

            return fECALSegmentationAlg->getStripLength(*this, localtemp, cID);
        }

        //----------------------------------------------------------------------------
        std::pair<TVector3, TVector3> GeometryCore::GetStripEnds(std::array<double, 3> const& point, const gar::raw::CellID_t &cID) const
        {
            //Get the matrix to make the transformation from Local to World
            std::array<double, 3> localtemp;
            gar::geo::LocalTransformation<TGeoHMatrix> trans;
            this->WorldToLocal(point, localtemp, trans);

            const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

            fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);

            std::pair<TVector3, TVector3> localStripEnds = fECALSegmentationAlg->getStripEnds(*this, localtemp, cID);

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
        std::pair<float, float> GeometryCore::CalculateLightPropagation(std::array<double, 3> const& point, const std::array<double, 3> &local, const gar::raw::CellID_t &cID) const
        {
            const std::array<double, 3> shape = this->FindShapeSize(this->FindNode(point));

            fECALSegmentationAlg->setLayerDimXY(shape[0] * 2, shape[1] * 2);

            return fECALSegmentationAlg->CalculateLightPropagation(*this, local, cID);
        }

        //----------------------------------------------------------------------------
        std::array<double, 3> GeometryCore::ReconstructStripHitPosition(const std::array<double, 3> &local, const float &xlocal, const gar::raw::CellID_t &cID) const
        {
            return fECALSegmentationAlg->ReconstructStripHitPosition(*this, local, xlocal, cID);
        }

        //----------------------------------------------------------------------------
        void GeometryCore::PrintGeometry() const
        {
            //Prints geometry parameters
            std::cout << "------------------------------" << std::endl;
            std::cout << "World Geometry" << std::endl;
            std::cout << "World Origin (x, y, z) " << GetWorldX() << " cm " << GetWorldY() << " cm " << GetWorldZ() << " cm" << std::endl;
            std::cout << "World Size (H, W, L) " << GetWorldHalfWidth() << " cm " << GetWorldHalfHeight() << " cm " << GetWorldLength() << " cm" << std::endl;
            std::cout << "------------------------------" << std::endl;
            std::cout << "Rock Geometry" << std::endl;
            std::cout << "Rock Origin (x, y, z) " << GetRockX() << " cm " << GetRockY() << " cm " << GetRockZ() << " cm" << std::endl;
            std::cout << "Rock Size (H, W, L) " << GetRockHalfWidth() << " cm " << GetRockHalfHeight() << " cm " << GetRockLength() << " cm" << std::endl;
            std::cout << "------------------------------" << std::endl;
            std::cout << "Enclosure Geometry" << std::endl;
            std::cout << "Enclosure Origin (x, y, z) " << GetEnclosureX() << " cm " << GetEnclosureY() << " cm " << GetEnclosureZ() << " cm" << std::endl;
            std::cout << "Enclosure Size (H, W, L) " << GetEnclosureHalfWidth() << " cm " << GetEnclosureHalfHeight() << " cm " << GetEnclosureLength() << " cm" << std::endl;
            std::cout << "------------------------------" << std::endl;

            if(GetLArTPCX() != 0 && GetLArTPCY() != 0 && GetLArTPCZ() != 0) {
                std::cout << "LArArgonCube Geometry" << std::endl;
                std::cout << "LArTPC Origin (x, y, z) " << GetLArTPCX() << " cm " << GetLArTPCY() << " cm " << GetLArTPCZ() << " cm" << std::endl;
                std::cout << "LArTPC Size (H, W, L) " << GetLArTPCHalfWidth() << " cm " << GetLArTPCHalfHeight() << " cm " << GetLArTPCLength() << " cm" << std::endl;
                std::cout << "------------------------------" << std::endl;
                std::cout << "LArActiveArgonCube Geometry" << std::endl;
                std::cout << "LArTPCActive Origin (x, y, z) " << GetActiveLArTPCX() << " cm " << GetActiveLArTPCY() << " cm " << GetActiveLArTPCZ() << " cm" << std::endl;
                std::cout << "LArTPCActive Size (H, W, L) " << GetActiveLArTPCHalfWidth() << " cm " << GetActiveLArTPCHalfHeight() << " cm " << GetActiveLArTPCLength() << " cm" << std::endl;
                std::cout << "------------------------------" << std::endl;
            }


            std::cout << "MPD Geometry" << std::endl;
            std::cout << "MPD Origin (x, y, z) " << GetMPDX() << " cm " << GetMPDY() << " cm " << GetMPDZ() << " cm" << std::endl;
            std::cout << "MPD Size (H, W, L) " << GetMPDHalfWidth() << " cm " << GetMPDHalfHeight() << " cm " << GetMPDLength() << " cm" << std::endl;
            std::cout << "------------------------------" << std::endl;
            std::cout << "TPC Geometry" << std::endl;
            std::cout << "TPC Origin (x, y, z) " << TPCXCent() << " cm " << TPCYCent() << " cm " << TPCZCent() << " cm" << std::endl;
            std::cout << "TPC Active Volume Size (R, L) " << TPCRadius() << " cm " << TPCLength() << " cm" << std::endl;
            std::cout << "------------------------------\n" << std::endl;


            std::cout << "ECAL Geometry" << std::endl;
            std::cout << "ECAL Barrel inner radius (Barrel): " << GetECALInnerBarrelRadius() << " cm" << std::endl;
            std::cout << "ECAL Barrel outer radius (Barrel): " << GetECALOuterBarrelRadius() << " cm" << std::endl;
            std::cout << "ECAL Barrel inner radius (Endcap): " << GetECALInnerEndcapRadius() << " cm" << std::endl;
            std::cout << "ECAL Barrel outer radius (Endcap): " << GetECALOuterEndcapRadius() << " cm" << std::endl;
            std::cout << "ECAL inner symmetry: " << GetECALInnerSymmetry() << std::endl;
            std::cout << "ECAL polyhedra angle: " << GetECALInnerAngle()*180/M_PI << " deg" << std::endl;
            std::cout << "ECAL polyhedra side length (Barrel): " << GetECALBarrelSideLength() << " cm" << std::endl;
            std::cout << "ECAL polyhedra apothem length (Barrel): " << GetECALBarrelApothemLength() << " cm" << std::endl;
            std::cout << "ECAL polyhedra side length (Endcap): " << GetECALEndcapSideLength() << " cm" << std::endl;
            std::cout << "ECAL polyhedra apothem length (Endcap): " << GetECALEndcapApothemLength() << " cm" << std::endl;
            std::cout << "ECAL Endcap Start X: " << GetECALEndcapStartX() << " cm" << std::endl;
            std::cout << "ECAL Endcap Outer X: " << GetECALEndcapOuterX() << " cm" << std::endl;
            std::cout << "Pressure Vessel Thickness: " << GetPVThickness() << " cm" << std::endl;
            std::cout << "------------------------------\n" << std::endl;
        }

        //----------------------------------------------------------------------------
        void GeometryCore::InitVariables()
        {
            fEnclosureX = 0.;
            fEnclosureY = 0.;
            fEnclosureZ = 0.;

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
            fECALSymmetry = -1;
            fECALEndcapStartX = 0.;

            fECALNodePath.clear();
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
