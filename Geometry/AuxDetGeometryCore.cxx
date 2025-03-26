/**
 * @file   AuxDetGeometryCore.cxx
 * @brief  Access the description of auxiliary detector geometry - implementation file
 * @author brebel@fnal.gov
 * @see    AuxDetGeometryCore.h
 *
 */

// class header
#include "Geometry/AuxDetGeometryCore.h"

#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetSensitiveGeo.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
// #include <Rtypes.h>

// C/C++ includes
#include <algorithm> // std::for_each(), std::transform()
#include <cctype>    // ::tolower()
#include <cmath>     // std::abs() ...
#include <cstddef>   // size_t
#include <limits>    // std::numeric_limits<>
#include <memory>    // std::default_deleter<>
#include <utility>   // std::swap()
#include <vector>

namespace gar {
  namespace geo {

    //......................................................................
    // Constructor.
    AuxDetGeometryCore::AuxDetGeometryCore(fhicl::ParameterSet const& pset)
      : fDetectorName(pset.get<std::string>("Name"))
    {
      std::transform(fDetectorName.begin(), fDetectorName.end(), fDetectorName.begin(), ::tolower);
    } // AuxDetGeometryCore::AuxDetGeometryCore()

    //......................................................................
    AuxDetGeometryCore::~AuxDetGeometryCore()
    {
      ClearGeometry();
    } // AuxDetGeometryCore::~AuxDetGeometryCore()

    //......................................................................
    void AuxDetGeometryCore::ApplyChannelMap(
      std::shared_ptr<geo::seg::AuxDetChannelMapAlg> pChannelMap)
    {
      pChannelMap->Initialize(fGeoData);
      fChannelMapAlg = pChannelMap;
    } // AuxDetGeometryCore::ApplyChannelMap()

    //......................................................................
    void AuxDetGeometryCore::LoadGeometryFile(std::string gdmlfile, std::string rootfile)
    {

      if (gdmlfile.empty()) {
        throw cet::exception("AuxDetGeometryCore") << "No GDML Geometry file specified!\n";
      }

      if (rootfile.empty()) {
        throw cet::exception("AuxDetGeometryCore") << "No ROOT Geometry file specified!\n";
      }

      ClearGeometry();

      // Open the GDML file, and convert it into ROOT TGeoManager format.
      // try to be efficient - if the GeometryCore object already imported
      // the file, then the gGeoManager will be non-null.  If not, import it.
      // Then lock the gGeoManager to prevent future imports.
      if (!gGeoManager) {
        TGeoManager::Import(rootfile.c_str());
        gGeoManager->LockGeometry();
      }

      std::vector<const TGeoNode*> path(8);
      path[0] = gGeoManager->GetTopNode();
      FindAuxDet(path, 0);

      fGDMLfile = gdmlfile;
      fROOTfile = rootfile;

      MF_LOG_INFO("AuxDetGeometryCore") << "New detector geometry loaded from "
                                        << "\n\t" << fROOTfile << "\n\t" << fGDMLfile << "\n";

    } // AuxDetGeometryCore::LoadGeometryFile()

    //......................................................................
    void AuxDetGeometryCore::ClearGeometry()
    {
      // auxiliary detectors
      std::for_each(AuxDets().begin(), AuxDets().end(), std::default_delete<AuxDetGeo>());
      AuxDets().clear();

    } // AuxDetGeometryCore::ClearGeometry()

    //......................................................................
    unsigned int AuxDetGeometryCore::NAuxDetSensitive(size_t const& aid) const
    {
      if (aid > NAuxDets() - 1)
        throw cet::exception("Geometry")
          << "Requested AuxDet index " << aid << " is out of range: " << NAuxDets();

      return AuxDets()[aid]->NSensitiveVolume();
    }

    //......................................................................
    //
    // Return the geometry description of the ith AuxDet.
    //
    // \param ad : input AuxDet number, starting from 0
    // \returns AuxDet geometry for ith AuxDet
    //
    // \throws geo::Exception if "ad" is outside allowed range
    //
    const AuxDetGeo& AuxDetGeometryCore::AuxDet(unsigned int const ad) const
    {
      if (ad >= NAuxDets())
        throw cet::exception("AuxDetGeometryCore") << "AuxDet " << ad << " does not exist\n";

      return *(AuxDets()[ad]);
    }

    //......................................................................
    unsigned int AuxDetGeometryCore::FindAuxDetAtPosition(double const worldPos[3]) const
    {
      return fChannelMapAlg->NearestAuxDet(worldPos, AuxDets());
    } // AuxDetGeometryCore::FindAuxDetAtPosition()

    //......................................................................
    const AuxDetGeo& AuxDetGeometryCore::PositionToAuxDet(double const worldLoc[3],
                                                          unsigned int& ad) const
    {
      // locate the desired Auxiliary Detector
      ad = this->FindAuxDetAtPosition(worldLoc);

      return this->AuxDet(ad);
    }

    //......................................................................
    void AuxDetGeometryCore::FindAuxDetSensitiveAtPosition(double const worldPos[3],
                                                           size_t& adg,
                                                           size_t& sv) const
    {
      adg = this->FindAuxDetAtPosition(worldPos);
      sv = fChannelMapAlg->NearestSensitiveAuxDet(worldPos, AuxDets(), adg);

      return;
    } // AuxDetGeometryCore::FindAuxDetAtPosition()

    //......................................................................
    const AuxDetSensitiveGeo& AuxDetGeometryCore::PositionToAuxDetSensitive(
      double const worldLoc[3],
      size_t& ad,
      size_t& sv) const
    {
      // locate the desired Auxiliary Detector
      this->FindAuxDetSensitiveAtPosition(worldLoc, ad, sv);
      return this->AuxDet(ad).SensitiveVolume(sv);
    }

    //......................................................................
    uint32_t AuxDetGeometryCore::PositionToAuxDetChannel(double const worldLoc[3],
                                                         size_t& ad,
                                                         size_t& sv) const
    {
      return fChannelMapAlg->PositionToAuxDetChannel(worldLoc, AuxDets(), ad, sv);
    }

    //......................................................................
    const TVector3 AuxDetGeometryCore::AuxDetChannelToPosition(uint32_t const& channel,
                                                               std::string const& auxDetName) const
    {
      return fChannelMapAlg->AuxDetChannelToPosition(channel, auxDetName, AuxDets());
    }

    //......................................................................
    const AuxDetGeo& AuxDetGeometryCore::ChannelToAuxDet(std::string const& auxDetName,
                                                         uint32_t const& channel) const
    {
      size_t adIdx = fChannelMapAlg->ChannelToAuxDet(AuxDets(), auxDetName, channel);
      return this->AuxDet(adIdx);
    }

    //......................................................................
    const AuxDetSensitiveGeo& AuxDetGeometryCore::ChannelToAuxDetSensitive(
      std::string const& auxDetName,
      uint32_t const& channel) const
    {
      auto idx = fChannelMapAlg->ChannelToSensitiveAuxDet(AuxDets(), auxDetName, channel);
      return this->AuxDet(idx.first).SensitiveVolume(idx.second);
    }

    //......................................................................
    void AuxDetGeometryCore::FindAuxDet(std::vector<const TGeoNode*>& path, unsigned int depth)
    {
      const char* nm = path[depth]->GetName();
      if ((strncmp(nm, "volAuxDet", 9) == 0)) {
        this->MakeAuxDet(path, depth);
        return;
      }

      //explore the next layer down
      unsigned int deeper = depth + 1;
      if (deeper >= path.size()) {
        throw cet::exception("AuxDetGeometryCore") << "exceeded maximum TGeoNode depth\n";
      }

      const TGeoVolume* v = path[depth]->GetVolume();
      int nd = v->GetNdaughters();
      for (int i = 0; i < nd; ++i) {
        path[deeper] = v->GetNode(i);
        this->FindAuxDet(path, deeper);
      }
    }

    //......................................................................
    void AuxDetGeometryCore::MakeAuxDet(std::vector<const TGeoNode*>& path, int depth)
    {
      AuxDets().push_back(new AuxDetGeo(path, depth));
    }

  } // namespace geo
} // gar
