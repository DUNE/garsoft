/**
* @file   GeometryCore.h
* @brief  Access the description of detector geometry
* @author brebel@fnal.gov
* @see    GeometryCore.cxx
*
* Structure of the header:
*
*     namespace geo {
*
*       // forward class declarations
*
*       namespace details {
*
*         // geometry iterator base class
*
*       }
*
*       // geometry iterators declaration
*       //  - cryostat_id_iterator
*       //  - TPC_id_iterator
*       //  - plane_id_iterator
*       //  - wire_id_iterator
*
*       // GeometryData_t definition (part of GeometryCore)
*
*       // GeometryCore declaration
*
*     }
*
*
*
* Revised <seligman@nevis.columbia.edu> 29-Jan-2009
*         Revise the class to make it into more of a general detector interface
* Revised <petrillo@fnal.gov> 27-Apr-2015
*         Factorization into a framework-independent GeometryCore.h and a
*         art framework interface
* Revised <petrillo@fnal.gov> 30-Apr-2015
*         Redesign of the iterators
*/
#ifndef GEO_GEOMETRYCORE_H
#define GEO_GEOMETRYCORE_H

// GArSoft libraries
#include "Geometry/LocalTransformation.h"
#include "RawDataProducts/CaloRawDigit.h"

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include <TVector3.h>
// #include <Rtypes.h>

// C/C++ standard libraries
#include <cstddef> // size_t
#include <string>
#include <vector>
#include <set>
#include <memory> // std::shared_ptr<>
#include <iterator> // std::forward_iterator_tag
#include <type_traits> // std::is_base_of<>
#include <map>
#include <utility>

// ROOT class prototypes
class TGeoManager;
class TGeoNode;
class TGeoMaterial;

/// Namespace collecting geometry-related classes utilities
namespace gar {
    namespace geo {

        class GeometryCore;
        namespace seg {
            class ChannelMapAlg;
            class SegmentationAlg;
        }

        typedef enum ROCType_ {HFILLER, IROC, IOROC, OOROC} ROCType;

        typedef struct ChanWithPos_
        {
            unsigned int id;
            TVector3 pos;
            TVector3 padrowdir;  // unit vector along pad row.  Sign not guaranteed.
            gar::geo::ROCType roctype;
        } ChanWithPos;

        typedef std::vector<gar::geo::ChanWithPos> ChanWithNeighbors;

        //Calorimeter Struct Data
        struct LayeredCalorimeterStruct
        {
            enum LayoutType {
                BarrelLayout = 0,
                EndcapLayout,
                ConicalLayout
            };

            LayoutType layoutType;
            double extent[6];
            int outer_symmetry;
            int inner_symmetry;
            double outer_phi0;
            double inner_phi0;
            double phi0;
            double gap0;
            double gap1;
            double gap2;

            struct Layer {
                Layer()
                : distance(0),
                phi0(0),
                absorberThickness(0),
                inner_nRadiationLengths(0),
                inner_nInteractionLengths(0),
                outer_nRadiationLengths(0),
                outer_nInteractionLengths(0),
                inner_thickness(0),
                outer_thickness(0),
                sensitive_thickness(0),
                cellSize0(0),
                cellSize1(0)
                {}

                double distance;
                double phi0;
                double absorberThickness;
                double inner_nRadiationLengths;
                double inner_nInteractionLengths;
                double outer_nRadiationLengths;
                double outer_nInteractionLengths;
                double inner_thickness;
                double outer_thickness;
                double sensitive_thickness;
                double cellSize0;
                double cellSize1;
            };

            std::vector<Layer> layers ;
        };

        typedef LayeredCalorimeterStruct LayeredCalorimeterData;

        //
        // iterators
        //
        namespace details {

            /// Base class for geometry iterators, containing some type definitions
            class geometry_iterator_types {
            public:

                //@{
                /// Structures to distinguish the constructors
                struct BeginPos_t {};
                struct EndPos_t {};
                struct UndefinedPos_t {};

                static constexpr BeginPos_t begin_pos = {};
                static constexpr EndPos_t end_pos = {};
                static constexpr UndefinedPos_t undefined_pos = {};
                //@}

            }; // class geometry_iterator_types

            /// Base class for geometry iterators (note: this is not an iterator)
            class geometry_iterator_base: public geometry_iterator_types {
            public:

                /// Constructor: associates with the specified geometry
                geometry_iterator_base(const gar::geo::GeometryCore *geom): pGeo(geom) {}

            protected:
                /// Returns a pointer to the geometry
                geo::GeometryCore const* geometry() const { return pGeo; }

                /// Default constructor; do not use a default-constructed iterator as-is!
                geometry_iterator_base() {}

            private:
                const gar::geo::GeometryCore *pGeo = nullptr; ///< pointer to the geometry

            }; // class geometry_iterator_base


            // forward declarations:
            template <typename GEOIDITER>
            class geometry_element_iterator;

            //@{
            /// Comparison operator: geometry ID and element point to the same ID
            template <typename GEOIDITER>
            bool operator== (
            geometry_element_iterator<GEOIDITER> const& iter,
            GEOIDITER const& id_iter
            );
            template <typename GEOIDITER>
            inline bool operator== (GEOIDITER const& id_iter,
            geometry_element_iterator<GEOIDITER> const& iter)
            { return iter == id_iter; }
            //@}

            //@{
            /// Comparison operator: geometry ID and element point to different IDs
            template <typename GEOIDITER>
            bool operator!= (
            geometry_element_iterator<GEOIDITER> const& iter,
            GEOIDITER const& id_iter
            );
            template <typename GEOIDITER>
            inline bool operator!= (
            GEOIDITER const& id_iter,
            geometry_element_iterator<GEOIDITER> const& iter
            )
            { return iter != id_iter; }
            //@}

            /**
            * @brief Forward iterator browsing all geometry elements in the detector
            * @tparam GEOITER type of geometry ID iterator
            *
            * This iterator works as the corresponding ID iterator in the template
            * argument. The difference is the dereferenciation operator: this one
            * obtains the geometry element directly, or throws on failure.
            * The boolean conversion operator checks that it can obtain a pointer to
            * the geometry element.
            *
            * In particular, get() and ID() methods still return the pointer to the
            * geometry element and its ID, respectively.
            *
            * It can also be initialized and compare with the corresponding ID
            * iterator.
            */
            template <typename GEOIDITER>
            class geometry_element_iterator:
                public std::forward_iterator_tag, public geometry_iterator_types
                {
                public:
                    using id_iterator_t = GEOIDITER;

                    static_assert(
                    std::is_base_of<geometry_iterator_base, id_iterator_t>::value,
                    "template class for geometry_element_iterator"
                    " must be a geometry iterator"
                    );

                    using iterator = geometry_element_iterator<id_iterator_t>; ///< this type

                    //@{
                    /// types mirrored from the ID iterator
                    using LocalID_t = typename id_iterator_t::LocalID_t;
                    using GeoID_t = typename id_iterator_t::GeoID_t;
                    using UndefinedPos_t = typename id_iterator_t::UndefinedPos_t;
                    using BeginPos_t = typename id_iterator_t::BeginPos_t;
                    using EndPos_t = typename id_iterator_t::EndPos_t;
                    using ElementPtr_t = typename id_iterator_t::ElementPtr_t;
                    //@}

                    //@{
                    /// Expose inherited constants
                    using geometry_iterator_types::undefined_pos;
                    using geometry_iterator_types::begin_pos;
                    using geometry_iterator_types::end_pos;
                    //@}

                    /// Geometry class pointed by the iterator
                    using Element_t = typename std::remove_pointer<ElementPtr_t>::type;

                    /// Default constructor; effect not defined: assign to it before using!
                    geometry_element_iterator() = default;

                    /// Constructor: points to begin
                    geometry_element_iterator(gar::geo::GeometryCore const* geom)
                    : id_iter(geom) {}

                    //@{
                    /// Constructor: points to the same element as the specified ID iterator
                    geometry_element_iterator(id_iterator_t const& iter): id_iter(iter) {}
                    geometry_element_iterator(id_iterator_t&& iter): id_iter(iter) {}
                    //@}

                    /// Constructor: points to the specified geometry element
                    geometry_element_iterator
                    (gar::geo::GeometryCore const* geom, GeoID_t const& start_from)
                    : id_iter(geom, start_from)
                    {}

                    /// Constructor: points to beginning
                    geometry_element_iterator
                    (gar::geo::GeometryCore const* geom, BeginPos_t const pos)
                    : id_iter(geom, pos)
                    {}

                    /// Constructor: points to end
                    geometry_element_iterator
                    (gar::geo::GeometryCore const* geom, EndPos_t const pos)
                    : id_iter(geom, pos)
                    {}

                    /// Returns true if the two iterators point to the same object
                    bool operator== (iterator const& as) const
                    { return id_iterator() == as.id_iterator(); }

                    /// Returns true if the two iterators point to different objects
                    bool operator!= (iterator const& as) const
                    { return id_iterator() != as.id_iterator(); }

                    /**
                    * @brief Returns the geometry element the iterator points to
                    * @return a constant reference to the element the iterator points to
                    * @throw cet::exception (category "geometry_iterator") if no valid
                    *   geometry element is currently pointed by the iterator
                    */
                    Element_t const& operator* () const
                    {
                        ElementPtr_t ptr = get();
                        if (ptr) return *ptr;
                        throw cet::exception("geometry_iterator")
                        << "iterator attempted to obtain geometry element "
                        << std::string(ID());
                    } // operator*()

                    /// Returns a pointer to the element the iterator points to (or nullptr)
                    Element_t const* operator-> () const { return get(); }

                    /// Prefix increment: returns this iterator pointing to the next element
                    iterator& operator++ () { ++id_iterator(); return *this; }

                    /// Postfix increment: returns the current iterator, then increments it
                    iterator operator++ (int)
                    { iterator old(*this); ++id_iterator(); return old; }

                    /// Returns whether the iterator is pointing to a valid geometry element
                    operator bool() const
                    { return bool(id_iterator()) && (id_iterator().get() != nullptr); }

                    /// Returns a pointer to the geometry element, or nullptr if invalid
                    ElementPtr_t get() const { return id_iterator().get(); }

                    /// Returns the ID of the pointed geometry element
                    LocalID_t const& ID() const { return *(id_iterator()); }

                protected:
                    friend bool geo::details::operator== <id_iterator_t>
                    (iterator const& iter, id_iterator_t const& id_iter);
                    friend bool geo::details::operator== <id_iterator_t>
                    (id_iterator_t const& id_iter, iterator const& iter);
                    friend bool geo::details::operator!= <id_iterator_t>
                    (iterator const& iter, id_iterator_t const& id_iter);
                    friend bool geo::details::operator!= <id_iterator_t>
                    (id_iterator_t const& id_iter, iterator const& iter);

                    //@{
                    /// Access to the base ID iterator
                    id_iterator_t const& id_iterator() const { return id_iter; }
                    id_iterator_t& id_iterator() { return id_iter; }
                    //@}

                private:
                    id_iterator_t id_iter; ///< iterator performing the job

                }; // class geometry_element_iterator<>

            } // namespace details

            template <typename Iter, Iter (GeometryCore::*BeginFunc)() const, Iter (GeometryCore::*EndFunc)() const>
            class IteratorBox {
            public:
                using iterator = Iter;

                IteratorBox(GeometryCore const* geom):
                    b((geom->*BeginFunc)()), e((geom->*EndFunc)()) {}

                    iterator begin() const { return b; }
                    iterator end() const { return e; }

                    iterator cbegin() const { return b; }
                    iterator cend() const { return e; }

                protected:
                    iterator b, e;
                }; // IteratorBox<>

                //
                // GeometryCore
                //

                /** **************************************************************************
                * @brief Description of geometry of one entire detector
                *
                * @note All lengths are specified in centimetres
                *
                *
                * How to correctly instantiate a GeometryCore object
                * ---------------------------------------------------
                *
                * Instantiation is a multi-step procedure:
                * 1. construct a GeometryCore object (the "service provider"),
                *    with the full configuration; at this step, configuration is just stored
                * 2. load a geometry with GeometryCore::LoadGeometryFile();
                *    this loads the detector geometry information
                * 3. prepare a channel map algorithm object (might use for example
                *    GeometryCore::DetectorName() or the detector geometry from the
                *    newly created object, but any use of channel mapping related functions
                *    is forbidden and it would yield undefined behaviour (expected to be
                *    catastrophic)
                * 4. acquire the channel mapping algorithm with
                *    GeometryCore::ApplyChannelMap(); at this point, the ChannelMapAlg object
                *    is asked to initialize itself and to perform whatever modifications to
                *    the geometry provider is needed.
                *
                * Step 3 (creation of the channel mapping algorithm object) can be performed
                * at any time before step 4, provided that no GeometryCore instance is needed
                * for it.
                *
                *
                * Configuration parameters
                * -------------------------
                *
                * - *Name* (string; mandatory): string identifying the detector; it can be
                *   different from the base name of the file used to initialize the geometry;
                *   standard names are recommended by each experiment.
                *   This name can be used, for example, to select which channel mapping
                *   algorithm to use.
                * - *SurfaceY* (real; mandatory): depth of the detector, in centimetrs;
                *   see SurfaceY() for details
            * - *MinWireZDist* (real; default: 3)
        * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance
        *   (see DefaultWiggle())
        *
        */
        class GeometryCore {
        public:
            // import iterators
            /**
            * @brief Initialize geometry from a given configuration
            * @param pset configuration parameters
            *
            * This constructor does not load any geometry description.
            * The next step is to do exactly that, by GeometryCore::LoadGeometryFile().
            */
            GeometryCore(fhicl::ParameterSet const& pset);

            /// Destructor
            ~GeometryCore();

            // You shall not copy or move or assign me!
            GeometryCore(GeometryCore const&) = delete;
            GeometryCore(GeometryCore&&) = delete;
            GeometryCore& operator= (GeometryCore const&) = delete;
            GeometryCore& operator= (GeometryCore&&) = delete;


            /**
            * @brief Returns the tolerance used in looking for positions
            * @return the tolerance value
            *
            * This parameter is used as tolerance ("wiggle") for methods that require
            * Typically, it's a additional fraction of tolerance: 0 means no tolerance,
            * 0.1 means 10% tolerance.
            *
            * @todo Confirm the definition of wiggle: this one is taken from other doc
            */
            double DefaultWiggle() const { return fPositionWiggle; }

            /**
            * @brief Returns the full directory path to the geometry file source
            * @return the full directory path to the geometry file source
            *
            * This is the full path of the source of the detector geometry GeometryCore
            * relies on.
            */
            std::string ROOTFile() const { return fROOTfile; }

            /**
            * @brief Returns the full directory path to the GDML file source
            * @return the full directory path to the GDML file source
            *
            * This is the full path of the source of the detector geometry handed to
            * the detector simulation (GEANT).
            */
            std::string GDMLFile() const { return fGDMLfile; }

            /// @name Detector information
            /// @{

            //
            // global features
            //
            /// Returns a string with the name of the detector, as configured
            std::string DetectorName() const { return fDetectorName; }


            //
            // position
            //

            /**
            * @brief Fills the arguments with the boundaries of the world
            * @param xlo (output) pointer to the lower x coordinate
            * @param xlo (output) pointer to the upper x coordinate
            * @param ylo (output) pointer to the lower y coordinate
            * @param ylo (output) pointer to the upper y coordinate
            * @param zlo (output) pointer to the lower z coordinate
            * @param zlo (output) pointer to the upper z coordinate
            * @throw cet::exception (`"GeometryCore"` category) if no world found
            *
            * This method fills the boundaries of the world volume, that is the one
            * known as `"volWorld"` in the geometry.
            *
            * If a pointer is null, its coordinate is skipped.
            *
            * @todo Replace it with a TPC boundaries style thing?
            * @todo Unify the coordinates type
            */
            void WorldBox(float* xlo, float* xhi,
            float* ylo, float* yhi,
            float* zlo, float* zhi) const;

            /**
            * @brief The position of the detector respect to earth surface
            * @return typical y position at surface in units of cm
            *
            * This is the depth (y) of the surface (where earth meets air) for this
            * detector site.
            * The number is expressed in world coordinates and in centimetres,
            * and it represents the y coordinate of earth surface.
            * A negative value means that the origin of coordinates, typically matching
            * the detector centre, is above surface.
            *
            * @todo check that this is actually how it is used
            */
            //
            double SurfaceY() const { return fSurfaceY; }


            //
            // object description and information
            //

            /// Access to the ROOT geometry description manager
            TGeoManager* ROOTGeoManager() const;

            float GetOriginX() const { return fOriginX; }

            float GetOriginY() const { return fOriginY; }

            float GetOriginZ() const { return fOriginZ; }

            float GetWorldX() const { return fWorldX; }

            float GetWorldY() const { return fWorldY; }

            float GetWorldZ() const { return fWorldZ; }

            float GetWorldHalfWidth() const  { return fWorldHalfWidth; }

            float GetWorldHalfHeight() const  { return fWorldHalfHeight; }

            float GetWorldLength() const  { return fWorldLength; }

            float GetRockX() const { return fRockX; }

            float GetRockY() const { return fRockY; }

            float GetRockZ() const { return fRockZ; }

            float GetRockHalfWidth() const  { return fRockHalfWidth; }

            float GetRockHalfHeight() const  { return fRockHalfHeight; }

            float GetRockLength() const  { return fRockLength; }

            float GetEnclosureX() const { return fEnclosureX; }

            float GetEnclosureY() const { return fEnclosureY; }

            float GetEnclosureZ() const { return fEnclosureZ; }

            float GetEnclosureHalfWidth() const  { return fEnclosureHalfWidth; }

            float GetEnclosureHalfHeight() const  { return fEnclosureHalfHeight; }

            float GetEnclosureLength() const  { return fEnclosureLength; }

            float GetMPDX() const { return fMPDX; }

            float GetMPDY() const { return fMPDY; }

            float GetMPDZ() const { return fMPDZ; }

            float GetMPDHalfWidth() const  { return fMPDHalfWidth; }

            float GetMPDHalfHeight() const  { return fMPDHalfHeight; }

            float GetMPDLength() const  { return fMPDLength; }

            float GetLArTPCX() const { return fLArTPCX; }

            float GetLArTPCY() const { return fLArTPCY; }

            float GetLArTPCZ() const { return fLArTPCZ; }

            float GetActiveLArTPCX() const { return fLArTPCXCent; }

            float GetActiveLArTPCY() const { return fLArTPCYCent; }

            float GetActiveLArTPCZ() const { return fLArTPCZCent; }

            float GetLArTPCHalfWidth() const { return fLArTPCHalfWidth; }

            float GetLArTPCHalfHeight() const { return fLArTPCHalfHeight; }

            float GetLArTPCLength() const { return fLArTPCLength; }

            float GetActiveLArTPCHalfWidth() const { return fLArTPCActiveHalfWidth; }

            float GetActiveLArTPCHalfHeight() const { return fLArTPCActiveHalfHeight; }

            float GetActiveLArTPCLength() const { return fLArTPCActiveLength; }

            unsigned int GetNLayers(std::string det) const;

            std::map< gar::geo::LayeredCalorimeterData::LayoutType, std::shared_ptr<gar::geo::LayeredCalorimeterData> > GetECALLayeredCalorimeterData() const { return fECALLayeredCalorimeterData; }

            float GArLiteXCent() const { return fGArLiteXCent; }

            float GArLiteYCent() const { return fGArLiteYCent; }

            float GArLiteZCent() const { return fGArLiteZCent; }

            float GArLiteRadius() const { return fGArLiteRadius; }

            float GArLiteLength() const { return fGArLiteLength; }

            /**
            * @brief Returns the name of the deepest volume containing specified point
            * @param point the location to query, in world coordinates
            * @return name of the volume containing the point
            *
            * @todo Use a reference to TVector3
            * @todo Use a double[3] instead?
            * @todo declare it const
            * @todo what happens if none?
            * @todo Unify the coordinates type
            */
            const std::string VolumeName(TVector3 const& point) const;


            /**
            * @brief Returns all the nodes with volumes with any of the specified names
            * @param vol_names list of names of volumes
            * @return list of nodes found
            *
            * All the nodes in the geometry are checked, and all the ones that contain
            * a volume with a name among the ones specified in vol_names are saved
            * in the collection and returned.
            */
            std::vector<TGeoNode const*> FindVolumePath(std::string const& vol_name) const;

            bool FindFirstVolume(std::string const& name, std::vector<const TGeoNode*>& path) const;

            void StoreECALNodes(std::map<std::string, std::vector<const TGeoNode*>> &map) const;

            /**
            * @brief Returns all the nodes with volumes with any of the specified names
            * @param vol_names list of names of volumes
            * @return list of nodes found
            *
            * All the nodes in the geometry are checked, and all the ones that contain
            * a volume with a name among the ones specified in vol_names are saved
            * in the collection and returned.
            */
            std::vector<TGeoNode const*> FindAllVolumes(std::set<std::string> const& vol_names) const;

            /**
            * @brief Returns paths of all nodes with volumes with the specified names
            * @param vol_names list of names of volumes
            * @return list paths of the found nodes
            *
            * All the nodes in the geometry are checked, and the path of all the ones
            * that contain a volume with a name among the ones specified in vol_names
            * is saved in the collection and returned.
            * A node path is a ordered list of all nodes leading to the final one,
            * starting from thetop level (root) down. The node at the `back()` of the
            * path is the one with name in vol_names.
            * No empty paths are returned.
            */
            std::vector<std::vector<TGeoNode const*>> FindAllVolumePaths(std::set<std::string> const& vol_names) const;

            //Return the node of the point
            template <typename T>
            TGeoNode* FindNode(T const &x, T const &y, T const &z) const;

            //Return the node of the point
            TGeoNode* FindNode(std::array<double, 3> const& point) const;

            //Return the node of the point
            TGeoNode* FindNode(TVector3 const& point) const;

            //Transform world coordinates in local coordinates
            bool WorldToLocal(std::array<double, 3> const& world, std::array<double, 3> &local, gar::geo::LocalTransformation<TGeoHMatrix> &trans) const;

            //Transform local coordinates in world coordinates
            bool LocalToWorld(std::array<double, 3> const& local, std::array<double, 3> &world, gar::geo::LocalTransformation<TGeoHMatrix> const &trans) const;

            /**
            * @brief Name of the deepest material containing the point xyz
            * @return material of the origin by default
            *
            * @todo make this constant
            * @todo remove return value constantness (or make it a reference)
            * @todo Unify the coordinates type
            */
            const std::string MaterialName(TVector3 const& point);


            /// Returns the material at the specified position
            /// @todo Unify the coordinates type
            TGeoMaterial const* Material(double x, double y, double z) const;

        /// Returns the total mass [kg] of the specified volume (default: world)
        /// @todo Use GetWorldVolumeName() as default instead
        double TotalMass(const char* vol = "volWorld") const;

        // this requires a bit more explanation about what this mass density is...
        /**
        * @brief Return the column density between two points
        * @param p1 pointer to array holding (x, y, z) of the first point
        * @param p2 pointer to array holding (x, y, z) of the second point
        * @return the mass
        *
        * Both points are specified in world coordinates.
        *
        * @todo Unify the coordinates type
        */
        /// Returns the mass between two coordinates
        double MassBetweenPoints(double *p1, double *p2) const;

        /// @}

        //
        // single object features
        //

        //@{
        /**
        * @brief Returns the radius of the TPC (y or z direction)
        * @return the radius of the TPC (y or z direction)
        */
        float TPCRadius() const { return fTPCRadius; }
        //@}

        //@{
        /**
        * @brief Returns the length of the TPC (x direction)
        * @return the value of the length of the specified TPC
        */
        float TPCLength() const { return fTPCLength; }
        //@}


        //@{
        /**
        * @brief Returns number of TPC drift volumes
        */
        int TPCNumDriftVols() const { return fTPCNumDriftVols; }
        //@}

        //@{
        /**
        * @brief Returns the X location of the center of the TPC in cm
        */
        float TPCXCent() const { return fTPCXCent; }
        //@}

        //@{
        /**
        * @brief Returns the Y location of the center of the TPC in cm
        */
        float TPCYCent() const { return fTPCYCent; }
        //@}

        //@{
        /**
        * @brief Returns the Z location of the center of the TPC in cm
        */
        float TPCZCent() const { return fTPCZCent; }
        //@}

        unsigned int NChannels() const;

        //
        // object description
        //

        //@{
        /**
        * @brief Return the name of GAr TPC volume
        * @return the name of the GArTPC volume
        *
        * This information is used by Geant4 simulation
        *
        */
        std::string GetGArTPCVolumeName() const {return "TPC_Drift"; }
        //@}

        //
        // geometry queries
        //

        //@{
        /**
        * @brief Returns the ID of the channel nearest to the specified position
        * @param worldLoc 3D coordinates of the point (world reference frame)
        * @return the ID of the channel, or raw::InvalidChannelID if invalid wire
        *
        * The different versions allow different way to provide the position.
        *
        * @todo remove the integers version
        * @todo Verify the raw::InvalidChannelID part
        */
        unsigned int NearestChannel(float              const  worldLoc[3]) const;
        unsigned int NearestChannel(std::vector<float> const& worldLoc)    const;
        unsigned int NearestChannel(TVector3           const& worldLoc)    const;
        void NearestChannelInfo(float const* xyz, gar::geo::ChanWithNeighbors &cwn) const;
        //@}

        //@{
        /**
        * @brief radii query methods passing through to the channel map algorithm
        */
        float GetIROCInnerRadius() const;
        float GetIROCOuterRadius() const;
        float GetOROCInnerRadius() const;
        float GetOROCOuterRadius() const;
        float GetOROCPadHeightChangeRadius() const;
        //@}


        //@{
        /**
        * @brief Returns the ID of the channel representing a gap
        * if you call NearestChannel and get this channel number, then charge is lost
        */

        unsigned int GapChannelNumber() const;

        void ChannelToPosition(unsigned int const channel, float* const worldLoc) const;

        //
        // unsorted methods
        //

        /**
        * @brief Returns whether a value is within the specified range
        * @param value the value to be tested
        * @param min the lower boundary
        * @param max the upper boundary
        * @return whether the value is within range
        *
        * If min is larger than max, they are swapped.
        * A tolerance of 10^-6 (absolute) is used.
        *
        * @todo Use wiggle instead of 10^-6
        * @todo resort source code for a bit of speed up
        */
        bool ValueInRange(double value, double min, double max) const;

        /// @name Geometry initialization
        /// @{

        /**
        * @brief Loads the geometry information from the specified files
        * @param gdmlfile path to file to be used for Geant4 simulation
        * @param rootfile path to file for internal geometry representation
        * @param bForceReload reload even if there is already a valid geometry
        * @see ApplyChannelMap()
        *
        * Both paths must directly resolve to an available file, as no search
        * is performed for them.
        *
        * The gdmlfile parameter does not have to necessarily be in GDML format,
        * as long as it's something supported by Geant4. This file is not used by
        * the geometry, but its path is provided on request by the simulation
        * modules (see GArSoft `LArG4` module).
        * The rootfile also does not need to be a ROOT file, but just anything
        * that TGeoManager::Import() supports. This file is parsed immediately
        * and the internal geometry representation is built out of it.
        *
        * @note After calling this method, the detector geometry information can
        * be considered complete, but the geometry service provider is not fully
        * initialized yet, since it's still necessary to provide or update the
        * channel mapping.
        */
        void LoadGeometryFile(std::string const& gdmlfile, std::string const& rootfile, bool bForceReload = false);

        /**
        * @brief Initializes the geometry to work with this channel map
        * @param pChannelMap a pointer to the channel mapping algorithm to be used
        * @see LoadGeometryFile()
        *
        * The specified channel mapping is used with this geometry.
        * The algorithm object is asked and allowed to make the necessary
        * modifications to the geometry description.
        * These modifications typically involve some resorting of the objects.
        *
        * The ownership of the algorithm object is shared, usually with a calling
        * framework: we maintain it alive as long as we need it (and no other code
        * can delete it), and we delete it only if no other code is sharing the
        * ownership.
        *
        * This method needs to be called after LoadGeometryFile() to complete the
        * geometry initialization.
        */
        void ApplyChannelMap(std::shared_ptr<gar::geo::seg::ChannelMapAlg> pChannelMap);
        /// @}

        void ApplyECALSegmentationAlg(std::shared_ptr<gar::geo::seg::SegmentationAlg> pECALSegmentationAlg);

        void ApplyMinervaSegmentationAlg(std::shared_ptr<gar::geo::seg::SegmentationAlg> pMinervaSegmentationAlg);

        void ApplyMuIDSegmentationAlg(std::shared_ptr<gar::geo::seg::SegmentationAlg> pMuIDSegmentationAlg);

        /// Returns the object handling the channel map
        gar::geo::seg::ChannelMapAlg const* ChannelMap() const { return fChannelMapAlg.get(); }

        /// Returns the object handling the ECAL segmentation
        gar::geo::seg::SegmentationAlg const* ECALSegmentationAlg() const { return fECALSegmentationAlg.get(); }

        ///Returns the object handling the Sc Tracker segmentation
        gar::geo::seg::SegmentationAlg const* MinervaSegmentationAlg() const { return fMinervaSegmentationAlg.get(); }

        ///Returns the object handling the MuID segmentation
        gar::geo::seg::SegmentationAlg const* MuIDSegmentationAlg() const { return fMuIDSegmentationAlg.get(); }

        //Returns the ECAL minimum radius of the Inner Barrel
        float GetECALInnerBarrelRadius() const { return fECALRinner; }

        //Returns the ECAL minimum radius of the Outer Barrel
        float GetECALOuterBarrelRadius() const { return fECALRouter; }

        //Returns the ECAL minimum radius of the Inner Endcap
        float GetECALInnerEndcapRadius() const { return fECALECapRinner; }

        //Returns the ECAL minimum radius of the Outer Endcap
        float GetECALOuterEndcapRadius() const { return fECALECapRouter; }

        //Returns the PV thickness
        float GetPVThickness() const { return fPVThickness; }
		
		//Returns the TPC end in drift direction (x)
		float GetTPCOuterX() const { return fTPCOuterX; }

        //Returns the number of sides of the barrel
        int GetECALInnerSymmetry() const { return fECALSymmetry; }

        //Returns the inner angle of the polyhedra
        float GetECALInnerAngle() const { return (( fECALSymmetry - 2 ) * M_PI / fECALSymmetry); }

        //Returns the side length of the polyhedra
        float GetECALBarrelSideLength() const { return ( 2 * std::sin( M_PI / fECALSymmetry ) * fECALRouter ); }

        //Returns the apothem length of the polyhedra
        float GetECALBarrelApothemLength() const { return ( std::cos( M_PI / fECALSymmetry ) * fECALRouter ); }

        //Returns the side length of the polyhedra
        float GetECALEndcapSideLength() const { return ( 2 * std::sin( M_PI / fECALSymmetry ) * fECALECapRouter ); }

        //Returns the apothem length of the polyhedra
        float GetECALEndcapApothemLength() const { return ( std::cos( M_PI / fECALSymmetry ) * fECALECapRouter ); }

        //Returns the ECAL start of the endcap
        float GetECALEndcapStartX() const { return fECALEndcapStartX; }

        //Returns the ECAL outer x of the endcap
        float GetECALEndcapOuterX() const { return fECALEndcapOuterX; }

        //Has Rock in Geometry
        bool HasRock() const { return fHasRock; }

        //Has Enclosure in Geometry
        bool HasEnclosure() const { return fHasEnclosure; }

        //Has LArTPC in Geometry
        bool HasLArTPCDetector() const { return fHasLArTPCDetector; }

        //Has Gas TPC in Geometry
        bool HasGasTPCDetector() const { return fHasGasTPCDetector; }

        //Has ECAL in Geometry
        bool HasECALDetector() const { return fHasECALDetector; }

        //ND-GAr lite
        bool HasTrackerScDetector() const { return fHasTrackerScDetector; }

        //Muon ID detector (only Barrel so far)
        bool HasMuonDetector() const { return fHasMuonDetector; }

        //Returns the MuID minimum radius of the Inner Barrel
        float GetMuIDInnerBarrelRadius() const { return fMuIDRinner; }

        //Returns the MuID minimum radius of the Outer Barrel
        float GetMuIDOuterBarrelRadius() const { return fMuIDRouter; }

        //Returns the number of sides of the barrel MuID
        int GetMuIDInnerSymmetry() const { return fMuIDSymmetry; }

        //Returns the inner angle of the polyhedra
        float GetMuIDInnerAngle() const { return (( fMuIDSymmetry - 2 ) * M_PI / fMuIDSymmetry); }

        //Returns the side length of the polyhedra MuID
        float GetMuIDBarrelSideLength() const { return ( 2 * std::sin( M_PI / fMuIDSymmetry ) * fMuIDRouter ); }

        //Returns the apothem length of the polyhedra MuID
        float GetMuIDBarrelApothemLength() const { return ( std::cos( M_PI / fMuIDSymmetry ) * fMuIDRouter ); }

        std::string GetWorldVolumeName() const { return "volWorld"; }

        bool PointInWorld(TVector3 const& point) const;

        bool PointInDetEnclosure(TVector3 const& point) const;

        bool PointInMPD(TVector3 const& point) const;

        bool PointInGArTPC(TVector3 const& point) const;

        bool PointInLArTPC(TVector3 const& point) const;

        bool PointInECALBarrel(TVector3 const& point) const;

        bool PointInECALEndcap(TVector3 const& point) const;

        bool PointInMuIDBarrel(TVector3 const& point) const;

        bool PointInMuIDEndcap(TVector3 const& point) const;

        float GetSensVolumeThickness(const TVector3& point) const;

        const std::array<double, 3> FindShapeSize(const TGeoNode *node) const;

        gar::raw::CellID_t GetCellID(const TGeoNode *node, const unsigned int& det_id, const unsigned int& stave, const unsigned int& module, const unsigned int& layer, const unsigned int& slice, const std::array<double, 3>& localPosition) const;

        const std::string GetECALCellIDEncoding() const;

        const std::string GetMinervaCellIDEncoding() const;

        const std::string GetMuIDCellIDEncoding() const;

        std::array<double, 3> GetPosition(const TGeoNode *node, const gar::raw::CellID_t &cID) const;

        bool isTile(const std::array<double, 3>& point, const gar::raw::CellID_t& cID) const;

        double getStripWidth(const std::array<double, 3>& point) const;

        double getTileSize(const std::array<double, 3>& point) const;

        double getStripLength(const std::array<double, 3>& point, const gar::raw::CellID_t &cID) const;

        double getIofRECAL() const { return fIofRECAL; }

        double getIofRMuID() const { return fIofRMuID; }

        double getIofRMinerva() const { return fIofRMinerva; }

        std::pair<TVector3, TVector3> GetStripEnds(const std::array<double, 3>& point, const gar::raw::CellID_t &cID) const;

        std::pair<float, float> CalculateLightPropagation(const std::array<double, 3>& point, const std::array<double, 3> &local, const gar::raw::CellID_t &cID) const;

        std::array<double, 3> ReconstructStripHitPosition(const std::array<double, 3>& point, const std::array<double, 3>& local, const float &xlocal, const gar::raw::CellID_t &cID) const;

    protected:

        /// Sets the detector name
        void SetDetectorName(std::string new_name) { fDetectorName = new_name; }

        void GetGeometryParameters();
        void FinalizeGeometryParameters();
        //Prints information on the detector geometry
        void PrintGeometry() const;

        void StoreTPCParameters();
        void GetDetectorsPresent();

    private:

        void InitVariables();

        /// Deletes the detector geometry structures
        void ClearGeometry();

        void StoreECALParameters();
        //Muon ID detector
        void StoreMuIDParameters();
        void StoreOtherParameters();

        void SetDetectorOrigin();

        bool FindWorldVolume();

        bool FindRockVolume();

        bool FindEnclosureVolume();

        bool FindMPDVolume();

        bool FindLArTPCVolume();

        bool FindActiveTPCVolume();

        bool FindActiveLArTPCVolume();

        bool FindGasTPCVolume();

        bool FindECALVolume();

        bool FindMuIDVolume();

        //Sets the ECAL inner barrel minimum radius
        bool FindECALInnerBarrelRadius();

        //Sets the ECAL outer barrel minimum radius
        bool FindECALOuterBarrelRadius();

        //Sets the ECAL inner endcap minimum radius
        bool FindECALInnerEndcapRadius();

        //Sets the ECAL outer endcap minimum radius
        bool FindECALOuterEndcapRadius();

        //Sets the PV thickness
        bool FindPVThickness();

        //Sets the TPC limit in x
        bool FindTPCOuterX();

        //Sets the number of sides of the ecal barrel
        bool FindECALInnerSymmetry();

        //Sets the position of the xplane of the ECAL endcap
        bool FindECALEndcapStartX();

        //Sets the position of the xplan of the end of the ECAL endcap
        bool FindECALEndcapOuterX();

        //Sets the number of layers in the ECAL
        bool FindECALnLayers();

        //Sets the ECAL Layered struct
        bool MakeECALLayeredCalorimeterData();

        //Sets the MuID inner barrel minimum radius
        bool FindMuIDInnerBarrelRadius();

        //Sets the MuID outer barrel minimum radius
        bool FindMuIDOuterBarrelRadius();

        //Sets the number of sides of the MuID barrel
        bool FindMuIDInnerSymmetry();

        //Sets the number of layers in the MuID
        bool FindMuIDnLayers();

        bool FindTrackerScVolume();

        bool FindTrackerScnPlanes();

        double         fSurfaceY;       ///< The point where air meets earth for this detector.
        std::string    fDetectorName;   ///< Name of the detector.
        std::string    fGDMLfile;       ///< path to geometry file used for Geant4 simulation
        std::string    fROOTfile;       ///< path to geometry file for geometry in GeometryCore
        double         fMinWireZDist;   ///< Minimum distance in Z from a point in which
        ///< to look for the closest wire
        double         fPositionWiggle; ///< accounting for rounding errors when testing positions
        bool           fPointInWarnings; ///< Generate warnings from failed inputs to PointIn* methods
        bool           fECALEndcapOutside; ///< Is the ECAL Endcap outside the PV

        float          fOriginX = 0.; //Origin of the ND-GAr
        float          fOriginY = 0.; //Origin of the ND-GAr
        float          fOriginZ = 0.; //Origin of the ND-GAr

        float          fTPCRadius = 0.;      ///< Radius of the TPC
        float          fTPCLength = 0.;      ///< length of the TPC

	    int            fTPCNumDriftVols = 2; ///< 2 if standard ALICE detector, 1 if single drift vol
        // we may need to add here a flag saying which side of the TPC the ROCs
        // are on if we are to support both of them.  And provide initialization
        // and a getter

        float          fTPCXCent = 0.;       ///< center of TPC: X
        float          fTPCYCent = 0.;       ///< center of TPC: Y
        float          fTPCZCent = 0.;       ///< center of TPC: Z

        float          fWorldX = 0.;
        float          fWorldY = 0.;
        float          fWorldZ = 0.;

        float          fRockX = 0.;
        float          fRockY = 0.;
        float          fRockZ = 0.;

        float          fEnclosureX = 0.;
        float          fEnclosureY = 0.;
        float          fEnclosureZ = 0.;

        float          fMPDX = 0.;
        float          fMPDY = 0.;
        float          fMPDZ = 0.;

        float          fLArTPCX = 0.;
        float          fLArTPCY = 0.;
        float          fLArTPCZ = 0.;

        float          fLArTPCXCent = 0.;
        float          fLArTPCYCent = 0.;
        float          fLArTPCZCent = 0.;

        float          fWorldHalfWidth = 0.;
        float          fWorldHalfHeight = 0.;
        float          fWorldLength = 0.;

        float          fRockHalfWidth = 0.;
        float          fRockHalfHeight = 0.;
        float          fRockLength = 0.;

        float          fEnclosureHalfWidth = 0.;
        float          fEnclosureHalfHeight = 0.;
        float          fEnclosureLength = 0.;

        float          fMPDHalfWidth = 0.;
        float          fMPDHalfHeight = 0.;
        float          fMPDLength = 0.;

        float          fLArTPCHalfWidth = 0.;
        float          fLArTPCHalfHeight = 0.;
        float          fLArTPCLength = 0.;

        float          fLArTPCActiveHalfWidth = 0.;
        float          fLArTPCActiveHalfHeight = 0.;
        float          fLArTPCActiveLength = 0.;

        float          fGArLiteXCent = 0;
        float          fGArLiteYCent = 0;
        float          fGArLiteZCent = 0;

        float          fGArLiteRadius = 0;
        float          fGArLiteLength = 0;

        std::map<std::string, std::vector<const TGeoNode*> > fECALNodePath; ///< Stored map of vectors of nodes for the ecal to speedup node searching

        bool fHasRock;
        bool fHasEnclosure;
        bool fHasLArTPCDetector;
        bool fHasGasTPCDetector;
        bool fHasECALDetector;
        bool fHasTrackerScDetector;

        // Indicies of refraction in scintillators.  Placed here because gdml won't have it.
        // Costructed as a separate parameter set in Geometry.fcl
        fhicl::ParameterSet fIndsOfRefract;
        float               fIofRECAL;
        float               fIofRMuID;
        float               fIofRMinerva;

        //Related to the ECAL
        float fECALRinner;              ///< Minimum radius of the ECAL inner barrel
        float fECALRouter;              ///< Minimum radius of the ECAL outer barrel
        float fECALECapRinner;
        float fECALECapRouter;
        float fPVThickness;             ///< Pressure Vessel thickness
        float fTPCOuterX;               ///< Position of the end xplane of the TPC;
        int fECALSymmetry;              ///< Number of sides of the Barrel
        float fECALEndcapStartX;        ///< Position of the start xplane of the ECAL endcap
        float fECALEndcapOuterX;        ///< Position of the end xplane of the ECAL endcap
        unsigned int fECALnLayers;      ///< number of ECAL layers from the seg algorithm

        //Related to GArLite
        unsigned int fTrackerScnPlanes;

        //Related to the MuID
        bool fHasMuonDetector;
        float fMuIDRinner;              ///< Minimum radius of the MuID inner barrel
        float fMuIDRouter;              ///< Minimum radius of the MuID outer barrel
        int fMuIDSymmetry;              ///< Number of sides of the MuID Barrel
        unsigned int fMuIDnLayers;      ///< number of MuID layers from the seg algorithm

        typedef std::shared_ptr<const gar::geo::seg::ChannelMapAlg> ChannelMapPtr;
        ChannelMapPtr  fChannelMapAlg;  ///< Object containing the channel to wire mapping

        typedef std::shared_ptr<const gar::geo::seg::SegmentationAlg> ECALSegmentationAlgPtr;
        ECALSegmentationAlgPtr fECALSegmentationAlg;  ///< Object containing the segmentation for the ECAL

        typedef std::shared_ptr<const gar::geo::seg::SegmentationAlg> MinervaSegmentationAlgPtr;
        MinervaSegmentationAlgPtr fMinervaSegmentationAlg;  ///< Object containing the segmentation for the Sc Tracker

        typedef std::shared_ptr<const gar::geo::seg::SegmentationAlg> MuIDSegmentationAlgPtr;
        MuIDSegmentationAlgPtr fMuIDSegmentationAlg;  ///< Object containing the segmentation for the Sc Tracker

        //ECAL layered struct containing all informations about the ECAL layers
        typedef std::map< gar::geo::LayeredCalorimeterData::LayoutType, std::shared_ptr<gar::geo::LayeredCalorimeterData> > ECALLayeredCalorimeterData;
        ECALLayeredCalorimeterData fECALLayeredCalorimeterData;

    }; // class GeometryCore

    template<> TGeoNode* GeometryCore::FindNode<float>(float const &x, float const &y, float const &z) const;
    template<> TGeoNode* GeometryCore::FindNode<double>(double const &x, double const &y, double const &z) const;

    /** **************************************************************************
    * @brief Iterator to navigate through all the nodes
    *
    * Note that this is not a fully standard forward iterator in that it lacks
    * of the postfix operator. The reason is that it's too expensive and it
    * should be avoided.
    * Also I did not bother declaring the standard type definitions
    * (that's just laziness).
    *
    * An example of iteration:
    *
    *     TGeoNode const* pCurrentNode;
    *
    *     ROOTGeoNodeForwardIterator iNode(geom->ROOTGeoManager()->GetTopNode());
    *     while ((pCurrentNode = *iNode)) {
    *       // do something with pCurrentNode
    *       ++iNode;
    *     } // while
    *
    * These iterators are one use only, and they can't be reset after a loop
    * is completed.
    */
    class ROOTGeoNodeForwardIterator {
    public:
        /// Constructor: start from this node
        ROOTGeoNodeForwardIterator(TGeoNode const* start_node){ init(start_node); }

        /// Returns the pointer to the current node, or nullptr if none
        TGeoNode const* operator* () const { return current_path.empty()? nullptr: current_path.back().self; }

        /// Points to the next node, or to nullptr if there are no more
        ROOTGeoNodeForwardIterator& operator++ ();

        /// Returns the full path of the current node
        std::vector<TGeoNode const*> get_path() const;

    protected:
        using Node_t = TGeoNode const*;
        struct NodeInfo_t {
            Node_t self; int sibling;
            NodeInfo_t(Node_t new_self, int new_sibling)
            : self(new_self), sibling(new_sibling) {}
        }; // NodeInfo_t

        /// which node, which sibling?
        std::vector<NodeInfo_t> current_path;

        void reach_deepest_descendant();

        void init(TGeoNode const* start_node);

    }; // class ROOTGeoNodeForwardIterator

} // namespace geo

} // namespace gar


//******************************************************************************
//***  template implementation
//***

//
// comparison operators between ID iterators and element iterators
//
template <typename GEOIDITER>
bool gar::geo::details::operator==
(geometry_element_iterator<GEOIDITER> const& iter, GEOIDITER const& id_iter)
{
    return iter.id_iterator() == id_iter;
} // operator==(iterator_t, id_iterator_t)

template <typename GEOIDITER>
bool gar::geo::details::operator!=
(geometry_element_iterator<GEOIDITER> const& iter, GEOIDITER const& id_iter)
{
    return iter.id_iterator() != id_iter;
} // operator!=(iterator_t, id_iterator_t)

//******************************************************************************

#endif // GEO_GEOMETRYCORE_H
