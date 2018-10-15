/**
 * @file   garsoft/Utilities/geo_vectors_utils_TVector.h
 * @brief  Specializations of `geo_vectors_utils.h` for ROOT old vector types.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   December 2, 2017
 * @see    `garsoft/Utilities/geo_vectors_utils.h`
 * @ingroup Geometry
 * 
 * This library provides facilities that can be used for both GArSoft geometry
 * vectors (`geo_vectors.h`) and ROOT `TVector3` and related, with the same
 * interface. 
 * 
 * This library depends on ROOT GenVector.
 * In the CET link list in `CMakeLists.txt`, link to `${ROOT_GENVECTOR}`.
 * 
 */

#ifndef GARSOFT_UTILITIES_GEO_VECTORS_UTILS_TVECTOR_H
#define GARSOFT_UTILITIES_GEO_VECTORS_UTILS_TVECTOR_H

#include "geo_vectors_utils.h"

// ROOT library
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"

// C/C++ standard library
#include <ostream>


namespace geo {
  namespace vect {
    
    // --- BEGIN TVector3 conversions ------------------------------------------
    /// @name TVector3 conversions
    /// @ingroup Geometry
    /// @{
    
    //--------------------------------------------------------------------------
    /// Converts a vector into a `TVector3`.
    template <typename Vector>
    TVector3 toTVector3(Vector const& v) { return convertTo<TVector3>(v); }
    
    /// @}
    //--------------------------------------------------------------------------
    
    /// Utilities to print vector types.
    namespace dump {
      
      // --- BEGIN Output of old-style ROOT vectors (TVector3 etc.) ------------
      /// @name Output of old-style ROOT vectors (TVector3 etc.)
      /// @ingroup Geometry
      /// @{
      
      /// Print a `TVector2` to an output stream.
      template <typename Stream>
      void Vector2(Stream&& out, TVector2 const& v)
        { out << "( " << v.X() << ", " << v.Y() << " )"; }
      
      /// Print a `TVector3` to an output stream.
      template <typename Stream>
      void Vector3(Stream&& out, TVector3 const& v)
        { out << "( " << v.X() << ", " << v.Y() << ", " << v.Z() << " )"; }
      
      /// Print a `TLorentzVector` to an output stream.
      template <typename Stream>
      void LorentzVector(Stream&& out, TLorentzVector const& v) {
        out 
          << "( " << v.X() << ", " << v.Y() << ", " << v.Z() << "; " << v.T()
          << " )";
      } // LorentzVector()
      
      /// Print a `TVector2` to an output stream.
      inline std::ostream& operator<< (std::ostream& out, TVector2 const& v)
        { Vector2(out, v); return out; }
      
      /// Print a `TVector3` to an output stream.
      inline std::ostream& operator<< (std::ostream& out, TVector3 const& v)
        { Vector3(out, v); return out; }
      
      /// Print a `TLorentzVector` to an output stream.
      inline std::ostream& operator<<
        (std::ostream& out, TLorentzVector const& v)
        { LorentzVector(out, v); return out; }
      
      // --- END Output of old-style ROOT vectors (TVector3 etc.) --------------
      /// @}
      
    } // namespace dump
    //--------------------------------------------------------------------------
    
    /// @}
    // --- END TVector3 conversions --------------------------------------------
    
  } // namespace vect
} // namespace geo


// The only way some generic code has to see the operator<< is for them to be
// exposed in the same namespace as the vectors they dump are; in TVector case,
// that's the global namespace... (Boost unit test checks still fail though)
using geo::vect::dump::operator<<;


//------------------------------------------------------------------------------
//--- Specialisations
//---
namespace geo {
  namespace vect {
    
    //--------------------------------------------------------------------------
    // Specialisations for: TVector2
    template <>
    inline auto mag2<TVector2>(TVector2 const& v) { return v.Mod2(); }
    
    //--------------------------------------------------------------------------
    
  } // namespace vect
} // namespace geo
//------------------------------------------------------------------------------


#endif // GARSOFT_UTILITIES_GEO_VECTORS_UTILS_TVECTOR_H
