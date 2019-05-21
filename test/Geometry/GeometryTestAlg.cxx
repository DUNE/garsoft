/**
 * @file   GeometryTestAlg.cxx
 * @brief  Unit test for geometry functionalities: implementation file
 * @date   2011/02/17
 * @author brebel@fnal.gov
 * @see    GeometryTestAlg.h
 */

// our header
#include "test/Geometry/GeometryTestAlg.h"

// GArSoft includes
#include "Geometry/GeometryCore.h"
#include "Geometry/AuxDetGeo.h"

// Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TGeoManager.h"
#include "TStopwatch.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TClass.h"

// C/C++ standard libraries
#include <cmath>
#include <vector>
#include <iterator> // std::inserter()
#include <algorithm> // std::copy()
#include <set>
#include <array>
#include <string>
#include <sstream>
#include <iostream>
#include <cassert>
#include <limits> // std::numeric_limits<>
#include <initializer_list>


namespace {
  template <typename T>
  inline T sqr(T v) { return v*v; }

  template <typename T>
  std::string to_string(const T& v) {
    std::ostringstream sstr;
    sstr << v;
    return sstr.str();
  } // ::to_string()

} // local namespace


namespace simple_geo {

  struct Point2D {
    double y = 0.;
    double z = 0.;

    Point2D() = default;
    Point2D(double new_y, double new_z): y(new_y), z(new_z) {}
  }; // struct Point2D

  Point2D operator+ (Point2D const& a, Point2D const& b)
    { return { a.y + b.y, a.z + b.z }; }
  Point2D operator* (Point2D const& p, double f)
    { return { p.y * f, p.z * f }; }
  Point2D operator/ (Point2D const& p, double f)
    { return { p.y / f, p.z / f }; }
  template <typename Stream>
  Stream& operator<< (Stream& out, Point2D const& p)
    { out << "( " << p.y << " ; " << p.z << " )"; return out; }

  class Area {
      public:
    Area() = default;

    Area(Point2D const& a, Point2D const& b)
      {
        set_sorted(min.y, max.y, a.y, b.y);
        set_sorted(min.z, max.z, a.z, b.z);
      } // Area(Point2D x2)

    Point2D const& Min() const { return min; }
    Point2D const& Max() const { return max; }
    Point2D Center() const { return (min + max) / 2; }
    double DeltaY() const { return Max().y - Min().y; }
    double DeltaZ() const { return Max().z - Min().z; }
    bool isEmpty() const { return (DeltaY() == 0) || (DeltaZ() == 0); }

    void IncludePoint(Point2D const& point)
      {
        set_min_max(min.y, max.y, point.y);
        set_min_max(min.z, max.z, point.z);
      } // Include()

    void Include(Area const& area)
      { IncludePoint(area.min); IncludePoint(area.max); }

    void Intersect(Area const& area)
      {
        set_max(min.y, area.min.y);
        set_min(max.y, area.max.y);
        set_max(min.z, area.min.z);
        set_min(max.z, area.max.z);
      }

      protected:
    Point2D min, max;

    void set_min(double& var, double val) { if (val < var) var = val; }
    void set_max(double& var, double val) { if (val > var) var = val; }
    void set_min_max(double& min_var, double& max_var, double val)
      { set_min(min_var, val); set_max(max_var, val); }
    void set_sorted(double& min_var, double& max_var, double a, double b)
      {
        if (a > b) { min_var = b; max_var = a; }
        else       { min_var = a; max_var = b; }
      }
  }; // class Area

} // namespace simple_geo

namespace gar{
  namespace geo{


      //......................................................................
    GeometryTestAlg::GeometryTestAlg(fhicl::ParameterSet const& pset)
    : geom(nullptr)
    {
      // initialize the list of non-fatal exceptions
      std::vector<std::string> NonFatalErrors(pset.get<std::vector<std::string>>
                                              ("ForgiveExceptions", std::vector<std::string>()));
      std::copy(NonFatalErrors.begin(), NonFatalErrors.end(),
                std::inserter(fNonFatalExceptions, fNonFatalExceptions.end()));

      // initialize the list of tests to be run
      //
      // our name selector accepts everything by default;
      // the default set skips the following:
      fRunTests.AddToDefinition("default",
                                "-CheckOverlaps",
                                "-ThoroughCheck");

      fRunTests.ParseNames("@default"); // let's start from default

      // read and parse the test list from the configuration parameters (if any)
      fRunTests.Parse(pset.get<std::vector<std::string>>("RunTests", {}));

      std::ostringstream sstr;
      fRunTests.PrintConfiguration(sstr);

      LOG_INFO("GeometryTestAlg") << "Test selection:" << sstr.str();

    } // GeometryTestAlg::GeometryTestAlg()

      //......................................................................
    unsigned int GeometryTestAlg::Run()
    {

      if (!geom) {
        throw cet::exception("GeometryTestAlg")
        << "GeometryTestAlg not configured: no valid geometry provided.\n";
      }

      unsigned int nErrors = 0; // currently unused

      // change the printed version number when changing the "GeometryTest" output
      LOG_VERBATIM("GeometryTest")
      << "GeometryTest version 1.0";

      LOG_VERBATIM("GeometryTest")
      << "\tRunning on detector: '"
      << geom->DetectorName()
      << "'"
      << "\n\tGeometry file: "
      << geom->ROOTFile();

      try{
        if (shouldRunTests("CheckOverlaps")) {
          LOG_INFO("GeometryTest") << "test for overlaps ...";
          gGeoManager->CheckOverlaps(1e-5);
          gGeoManager->PrintOverlaps();
          if (!gGeoManager->GetListOfOverlaps()->IsEmpty()) {
            mf::LogError("GeometryTest")
            << gGeoManager->GetListOfOverlaps()->GetSize()
            << " overlaps found in geometry during overlap test!";
            ++nErrors;
          }
          LOG_INFO("GeometryTest") << "complete.";
        }

        if (shouldRunTests("ThoroughCheck")) {
          LOG_INFO("GeometryTest") << "thorough geometry test ...";
          gGeoManager->CheckGeometryFull();
          if (!gGeoManager->GetListOfOverlaps()->IsEmpty()) {
            mf::LogError("GeometryTest")
            << gGeoManager->GetListOfOverlaps()->GetSize()
            << " overlaps found in geometry during thorough test!";
            ++nErrors;
          }
          LOG_INFO("GeometryTest") << "complete.";
        }

        if (shouldRunTests("FindVolumes")) {
          LOG_INFO("GeometryTest") << "test FindAllVolumes method ...";
          testFindVolumes();
          LOG_INFO("GeometryTest") << "complete.";
        }

        if (shouldRunTests("NearestWire")) {
          LOG_INFO("GeometryTest") << "testNearestChannel...";
          testNearestChannel();
          LOG_INFO("GeometryTest") << "complete.";
        }

        if (shouldRunTests("Stepping")) {
          LOG_INFO("GeometryTest") << "testStepping...";
          testStepping();
          LOG_INFO("GeometryTest") << "complete.";
        }

      }
      catch (cet::exception &e) {
        mf::LogWarning("GeometryTest") << "exception caught: \n" << e;
        if (fNonFatalExceptions.count(e.category()) == 0) throw;
      }

      std::ostringstream out;
      if (!fRunTests.CheckQueryRegistry(out)) {
        throw cet::exception("GeometryTest")
        << "(postumous) configuration error detected!\n"
        << out.rdbuf();
      }

      mf::LogInfo log("GeometryTest");
      log << "Tests completed:";
      auto tests_run = fRunTests.AcceptedNames();
      if (tests_run.empty()) {
        log << "\n  no test run";
      }
      else {
        log << "\n  " << tests_run.size() << " tests run:\t ";
        for (std::string const& test_name: tests_run) log << " " << test_name;
      }
      auto tests_skipped = fRunTests.RejectedNames();
      if (!tests_skipped.empty()) {
        log << "\n  " << tests_skipped.size() << " tests skipped:\t ";
        for (std::string const& test_name: tests_skipped) log << " " << test_name;
      }

      return nErrors;
    } // GeometryTestAlg::Run()



    //......................................................................
    void GeometryTestAlg::printChannelSummary()
    {
      uint32_t channels = geom->NChannels();
      LOG_VERBATIM("GeometryTest")
      << "there are "
      << channels
      << " channels";

      return;
    }

    //......................................................................
    // great sanity check for geometry, only call in analyze when debugging
    void GeometryTestAlg::printDetDim()
    {
      LOG_VERBATIM("GeometryTest")
      << "  TPC:    width: "
      << geom->TPCHalfWidth()
      << "    height: "
      << geom->TPCHalfHeight()
      << "    length: "
      << geom->TPCLength();

      TVector3 origin(0., 0., 0.);
      TVector3 outside(1.e6, 1.e6, 1.e6);

      LOG_VERBATIM("GeometryTest")
      << "test check of point in world: \n(0, 0, 0) is in volume "
      << geom->VolumeName(origin)
      << "\n (1e6, 1e6, 1e6) is in volume "
      << geom->VolumeName(outside);

      return;
    }

    //--------------------------------------------------------------------------
    void GeometryTestAlg::printAllGeometry() const {
      LOG_VERBATIM("GeometryTest")
      << "Detector " << geom->DetectorName();
    } // GeometryTestAlg::printAllGeometry()


    //......................................................................
    unsigned int GeometryTestAlg::testFindWorldVolumes() {

      unsigned int nErrors = 0;

      std::set<std::string> volume_names;
      std::vector<TGeoNode const*> nodes;

      // world
      volume_names.insert(geom->GetWorldVolumeName());
      nodes = geom->FindAllVolumes(volume_names);
      {
        mf::LogVerbatim log("GeometryTest");
        log
        << "Found "
        << nodes.size()
        << " world volumes '"
        << geom->GetWorldVolumeName() << "':";
        for (TGeoNode const* node: nodes) {
          TGeoVolume const* pVolume = node->GetVolume();
          log << "\n - '" << pVolume->GetName() << "' (a "
          << pVolume->GetShape()->GetName() << ")";
        } // for
      } // anonymous block
      if (nodes.size() != 1) {
        ++nErrors;
        mf::LogError("GeometryTest")
        << "Found " << nodes.size() << " world volumes '"
        << geom->GetWorldVolumeName() << "! [expecting: one!!]";
      } // if nodes

      return nErrors;
    } // GeometryTestAlg::testFindWorldVolumes()

    //--------------------------------------------------------------------------
    unsigned int GeometryTestAlg::testFindTPCvolumePaths() {

      unsigned int nErrors = 0;

      // search the full path of all TPCs
      std::set<std::string> volume_names({"volTPC"});
      std::vector<std::vector<TGeoNode const*>> node_paths = geom->FindAllVolumePaths(volume_names);

      mf::LogVerbatim log("GeometryTest");
      log
      << "Found "
      << node_paths.size()
      << " TPC volumes:";
      for (auto const& path: node_paths) {
        TGeoNode const* node = path.back();
        TGeoVolume const* pVolume = node->GetVolume();
        log << "\n - '" << pVolume->GetName() << "' (a "
        << pVolume->GetShape()->GetName() << ") with " << (path.size()-1)
        << " ancestors";
        for (TGeoNode const* pNode: path) {
          TGeoVolume const* pVolume = pNode->GetVolume();
          log << "\n      * '" << pVolume->GetName() << "' (a "
          << pVolume->GetShape()->GetName() << ") with a "
          << pNode->GetMatrix()->IsA()->GetName() << " that "
          << (pNode->GetMatrix()->IsTranslation()? "is": "is not")
          << " a simple translation";
        } // for node
      } // for path

      return nErrors;
    } // GeometryTestAlg::testFindTPCvolumePaths()

    //--------------------------------------------------------------------------
    void GeometryTestAlg::testFindVolumes() {
      /*
       * Finds and checks a selected number of volumes by name:
       * - world volume
       * - cryostat volumes
       * - TPCs (returns full paths)
       */

      unsigned int nErrors = 0;

      nErrors += testFindWorldVolumes();
      nErrors += testFindTPCvolumePaths();

      if (nErrors != 0) {
        throw cet::exception("FindVolumes")
        << "Collected " << nErrors << " errors during FindAllVolumes() test!\n";
      }

    } // GeometryTestAlg::testFindVolumes()


    //......................................................................
    void GeometryTestAlg::testNearestChannel()
    {
      // Even if you comment it out, please leave the TStopWatch code
      // in this code for additional testing. The NearestChannel routine
      // is the most frequently called in the simulation, so its execution time
      // is an important component of GArSoft's speed.
      TStopwatch stopWatch;
      stopWatch.Start();

      float posWorld[3] = {0.};
      posWorld[1] = 0.5 * geom->TPCHalfHeight();
      posWorld[2] = 0.5 * geom->TPCLength();

      try{

        LOG_VERBATIM("GeometryTestAlg")
        << " find nearest channel to point ("
        << posWorld[0]
        << ", "
        << posWorld[1]
        << ", "
        << posWorld[2]
        << ")";

        // The float[] version tested here is used by the TVector3 version, so this test both.
        unsigned int nearest = geom->NearestChannel(posWorld);

        LOG_VERBATIM("GeometryTestAlg")
        << "\t ...nearest channel to point: "
        << nearest;

        // We also want to test the std::vector<float> version
        std::array<float, 3> posWorldV;
        for (int i = 0; i < 3; ++i) posWorldV[i] = posWorld[i];

        nearest = geom->NearestChannel(posWorldV.data());

        LOG_VERBATIM("GeometryTestAlg")
        << "\t ...nearest channel to point: "
        << nearest;

      }
      catch(cet::exception &e){
        mf::LogWarning("GeoTestCaughtException") << e;
        if (fNonFatalExceptions.count(e.category()) == 0) throw;
      }

      stopWatch.Stop();
      LOG_DEBUG("GeometryTest") << "\tdone testing nearest channel";
      stopWatch.Print();

      // trigger an exception with NearestChannel
      LOG_VERBATIM("GeometryTest")
      << "\tattempt to cause an exception to be caught "
      << "when looking for a nearest channel";

      // pick a position out of the world
      geom->WorldBox(nullptr, posWorld + 0,
                     nullptr, posWorld + 1,
                     nullptr, posWorld + 2);
      for (int i = 0; i < 3; ++i) posWorld[i] *= 2.;

      bool hasThrown = false;
      unsigned int nearest_to_what = 0;
      try{
        nearest_to_what = geom->NearestChannel(posWorld);
      }
      catch(const cet::exception& e){
        LOG_WARNING("GeoTestCaughtException")
        << "caught execpetion: "
        << e;
        hasThrown = true;
      }

      if (!hasThrown) {
        LOG_WARNING("GeoTestErrorNearestChannel")
        << "GeometryCore::NearestChannel() did not raise an exception on out-of-world position ("
        << posWorld[0]
        << "; "
        << posWorld[1]
        << "; "
        << posWorld[2]
        << "), and returned "
        << nearest_to_what
        << " instead.\n This is normally considered a failure.";
      }

      return;
    }

    //......................................................................
    void GeometryTestAlg::testStepping()
    {
      //
      // Test stepping.
      //
      double xyz[3]  = {0.};
      double dxyz[3] = {0., 0., 1.};


      LOG_VERBATIM("GeometryTest")
      << "initial"
      << "\n\tposition:"  << xyz[0]  << "\t" << xyz[1]  << "\t" << xyz[2]
      << "\n\tdirection:" << dxyz[0] << "\t" << dxyz[1] << "\t" << dxyz[2];

      gGeoManager->InitTrack(xyz, dxyz);
      for (int i=0; i<10; ++i) {
        const double* pos = gGeoManager->GetCurrentPoint();
        const double* dir = gGeoManager->GetCurrentDirection();
        LOG_VERBATIM("GeometryTest") << "\tnode = "
        << gGeoManager->GetCurrentNode()->GetName()
        << "\n\t\tpos=" << "\t"
        << pos[0] << "\t"
        << pos[1] << "\t"
        << pos[2]
        << "\n\t\tdir=" << "\t"
        << dir[0] << "\t"
        << dir[1] << "\t"
        << dir[2]
        << "\n\t\tmat = "
        << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();

        gGeoManager->FindNextBoundary();
        gGeoManager->FindNormal();
        gGeoManager->Step(kTRUE,kTRUE);
      }

      xyz[0] = 306.108; xyz[1] = -7.23775; xyz[2] = 856.757;
      gGeoManager->InitTrack(xyz, dxyz);
      LOG_VERBATIM("GeometryTest") << "\tnode = "
      << gGeoManager->GetCurrentNode()->GetName()
      << "\n\tmat = "
      << gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->GetName();

      gGeoManager->GetCurrentNode()->GetVolume()->GetMaterial()->Print();

    }

    //......................................................................
    inline bool GeometryTestAlg::shouldRunTests(std::string test_name) const {
      return fRunTests(test_name);
    } // GeometryTestAlg::shouldRunTests()

  }//end namespace
} // gar
