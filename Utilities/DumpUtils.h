/**
 * @file   DumpUtils.h
 * @brief  Utilities to dump objects into a stream.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   February 28, 2017
 *
 * The library `DumpUtils.h` is currently header-only.
 *
 */

#ifndef GARSOFT_COREUTILS_DUMPUTILS_H
#define GARSOFT_COREUTILS_DUMPUTILS_H

// C++ libraries
#include <string>
#include <sstream>
#include <type_traits>


namespace gar {
  
  
  /// Namespace for GArSoft dumping utilities.
  namespace dump {
    
    namespace details {
      
      template <typename Coll>
      auto ptr_cbegin(Coll const& v) { using std::cbegin; return cbegin(v); }
      
      template <typename T>
      std::add_const_t<T>* ptr_cbegin(T* ptr) { return ptr; }
      
      
      /// Inserts `n` of elements of `a` in the specified stream.
      template <typename Stream, typename Array>
      void dumpArray(Stream&& out, Array&& a, size_t n) {
        out << "{";
        if (n == 0) { out << "}"; return; }
        auto it = ptr_cbegin(a);
        out << " " << *it;
        std::size_t i = 0;
        while (++i < n) out << "; " << (*++it);
        out << " }";
      } // dumpArray()
      
    } // namespace details
    
    /**
     * @brief Dumps the first N elements of an array.
     * @tparam Array the array type (supporting indexing operator)
     * 
     * The dump is producer by the call operator, in the form
     *     
     *     { a1, a2, ... aN }
     *     
     * 
     */
    template <typename Array>
    struct ArrayDumper {
      using Array_t = Array;
      using This_t = ArrayDumper<Array_t>;
      
      Array_t const& a; ///< A reference to the array to be printed.
      size_t n; ///< Number of elements to be printed.
      
      ArrayDumper(Array_t const& a, size_t n): a(a), n(n) {}
      
      // constructors ahead
      ArrayDumper(This_t const& from) = default;
      ArrayDumper(This_t&& from) = default;
      ArrayDumper& operator=(This_t const& from) = delete;
      ArrayDumper& operator=(This_t&& from) = delete;
      
      /// Inserts the content of the referenced array into the specified stream.
      template <typename Stream>
      void operator() (Stream&& out) const
        { details::dumpArray(std::forward<Stream>(out), a, n); }
      
      /// Converts the content of the stored vector into a string.
      explicit operator std::string() const
        { std::ostringstream sstr; this->operator()(sstr); return sstr.str(); }
      
    }; // struct ArrayDumper
    
    
    template <typename T>
    struct ArrayDumper<T*> {
      using Array_t = T*;
      using This_t = ArrayDumper<Array_t>;
      
      Array_t a; ///< A reference to the array to be printed.
      size_t n; ///< Number of elements to be printed.
      
      ArrayDumper(Array_t a, size_t n): a(a), n(n) {}
      
      /// Inserts the content of the referenced array into the specified stream.
      template <typename Stream>
      void operator() (Stream&& out) const
        { details::dumpArray(std::forward<Stream>(out), a, n); }
      
      /// Converts the content of the stored vector into a string.
      explicit operator std::string() const
        { std::ostringstream sstr; this->operator()(sstr); return sstr.str(); }
      
    }; // struct ArrayDumper<T*>
    
    
    /**
     * @brief Manipulator managing the dump of the vector content into a stream.
     * @tparam Vector the type of vector to be dumped
     * 
     * The manipulator is constructed with the object that needs to be dumped.
     * This implementation requires a vector of at least 3D with methods `X()`,
     * `Y()` and `Z()`.
     * The dump is producer by the call operator, in the form
     *     
     *     { vx, vy, vz }
     *     
     * 
     * An overloaded insertion operator will call the `operator()` method of
     * the manipulator to perform the insertion.
     * 
     * This class can be specialised to manage different vector types. 
     * For example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * #include <array>
     * 
     * namespace gar {
     *   namespace dump {
     *     
     *     // specialization for 3D std::array<>
     *     template <typename T>
     *     struct VectorDumper<std::array<T, 3U>> {
     *       
     *       using Vector_t = std::array<T, 3U>;
     *       
     *       Vector_t const& a;
     *       
     *       VectorDumper(Vector_t const& a): a(a) {}
     *       
     *       template <typename Stream>
     *       void operator() (Stream&& out) const
     *         { out << "{ " << a[0] << "; " << a[1] << "; " << a[2] << " }"; }
     *       
     *     }; // struct VectorDumper<array>
     *     
     *     
     *   } // namespace dump
     * } // namespace gar
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     */
    template <typename Vector>
    struct VectorDumper {
      using Vector_t = Vector;
      using This_t = VectorDumper<Vector_t>;
      
      Vector_t const& v; ///< A reference to the vector to be printed.
      
      explicit VectorDumper(Vector_t const& v): v(v) {}
      
      // constructors ahead
      VectorDumper(This_t const& from) = default;
      VectorDumper(This_t&& from) = default;
      VectorDumper& operator=(This_t const& from) = delete;
      VectorDumper& operator=(This_t&& from) = delete;
      
      /// Inserts the content of the stored vector into the specified stream.
      template <typename Stream>
      void operator() (Stream&& out) const
        { out << "{ " << v.X() << "; " << v.Y() << "; " << v.Z() << " }"; }
      
      /// Converts the content of the stored vector into a string.
      explicit operator std::string() const
        { std::ostringstream sstr; this->operator()(sstr); return sstr.str(); }
      
    }; // struct VectorDumper<>
    
    
    // Specialization for bare pointers.
    template <typename T>
    struct VectorDumper<T*>: public ArrayDumper<T const*> {
      explicit VectorDumper(T* v): ArrayDumper<T const*>(v, 3U) {}
    }; // VectorDumper<T*>
    
    // Specialization for C-style arrays.
    template <typename T>
    struct VectorDumper<T[3]>: public ArrayDumper<T const*> {
      explicit VectorDumper(T const* v): ArrayDumper<T const*>(v, 3U) {}
    }; // VectorDumper<T*>
    
    
    
    /**
     * @brief Returns a manipulator which will print the specified array.
     * @tparam N the number of entries to be printed
     * @tparam Array the type of array to be printed
     * @param a the array to be printed
     * @return a manipulator
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * double const data[5] = { 1., 2., 3., 4., 6. };
     * std::cout << "Data: " << gar::dump::array<5>(data) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Data: { 1; 2; 3; 4; 6 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * Addition of specific classes can be achieved via specialization of
     * the class `ArrayDumper`.
     * 
     * 
     * Requirements
     * -------------
     * 
     * * the type `Array` is required to respond to a global function `cbegin()`
     *   (like `std::cbegin()`) with a forward iterator
     * * the value type of the iterator returned by the aforementioned `Array`
     *   method must have its insertion operator into a stream defined
     * * the array `a` must not fall out of scope until the insertion into
     *   the stream happens
     * 
     */
    template <size_t N, typename Array>
    auto array(Array const& a) { return ArrayDumper<Array>(a, N); }
    
    
    /**
     * @brief Returns a manipulator which will print the specified array.
     * @tparam Vector type of vector to be printed
     * @param v the vector to be printed
     * @return a manipulator
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * std::vector<double> data = { 1., 2., 3., 4., 6. };
     * std::cout << "Data: " << gar::dump::vector(data) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Data: { 1; 2; 3; 4; 6 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * Addition of specific classes can be achieved via specialization of
     * the class `ArrayDumper`.
     * 
     * 
     * Requirements
     * -------------
     * 
     * * a method `Vector::size() const` must be available returning the number
     *   of elements in the vector
     * * the type `Vector` is required to respond to a global function
     *   `cbegin()` (like `std::cbegin()`) with a forward iterator
     * * the value type returned by the aforementioned `Vector` method must have
     *   its insertion operator into a stream defined
     * * the vector `v` must not fall out of scope until the insertion into
     *   the stream happens
     * 
     */
    template <typename Vector>
    auto vector(Vector const& v) { return ArrayDumper<Vector>(v, v.size()); }
    
    
    
    /**
     * @brief Returns a manipulator which will print the specified vector.
     * @tparam Vector3D the type of 3D vector to be printed
     * @param v the vector to be printed
     * @return a manipulator
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * TVector3 pos;
     * std::cout << "Position: " << gar::dump::vector3D(pos) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Position: { 0; 0; 0 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * Addition of specific classes can be achieved via specialization of
     * the class `VectorDumper`.
     * 
     * 
     * Requirements
     * -------------
     * 
     * * the type `Vector3D` is required to provide constant accessor methods
     *   returning the vector components, called `X()`, `Y()` and `Z()`
     * * the values returned by the aforementioned `Vector3D` methods must have
     *   their insertion operator into a stream defined
     * * the vector `v` must not fall out of scope until the insertion into
     *   the stream happens
     * 
     */
    template <typename Vector3D>
    auto vector3D(Vector3D const& v) { return VectorDumper<Vector3D>(v); }
    
    
    
    /**
     * @brief Dumps the array contained in the manipulator into a stream.
     * @tparam Stream the type of output stream
     * @tparam Array the type of array to be printed
     * @tparam NPrint the number of element of the array to be printed
     * @param out the output stream
     * @param manip the manipulator containing the array to be printed
     * @return a reference to the output stream
     * @see ArrayDumper
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * double const data[5] = { 1., 2., 3., 4., 6. };
     * std::cout << "Position: " << gar::dump::array<5>(data) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Position: { 1, 2, 3, 4, 5 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * 
     * Requirements
     * -------------
     * 
     * * the type `Vector3D` is required to provide constant accessor methods
     *   returning the vector components, called `X()`, `Y()` and `Z()`
     * * the values returned by the aforementioned `Vector3D` methods must have
     *   their insertion operator into `Stream` defined
     * 
     */
    template <typename Stream, typename Array>
    Stream& operator<< (Stream&& out, ArrayDumper<Array>&& manip)
      { manip(std::forward<Stream>(out)); return out; }
    
    
    /**
     * @brief Dumps the vector contained in the manipulator into a stream.
     * @tparam Stream the type of output stream
     * @tparam Vector3D the type of 3D vector to be printed
     * @param out the output stream
     * @param manip the manipulator containing the vector to be printed
     * @return a reference to the output stream
     * @see VectorDumper
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * TVector3 pos;
     * std::cout << "Position: " << gar::dump::vector3D(pos) << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Position: { 0; 0; 0 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * 
     * Requirements
     * -------------
     * 
     * * the type `Vector3D` is required to provide constant accessor methods
     *   returning the vector components, called `X()`, `Y()` and `Z()`
     * * the values returned by the aforementioned `Vector3D` methods must have
     *   their insertion operator into Stream defined
     * 
     */
    template <typename Stream, typename Vector>
    Stream& operator<< (Stream&& out, VectorDumper<Vector>&& manip)
      { manip(std::forward<Stream>(out)); return out; }
    
    /**
     * @brief Concatenates a vector to the specified string.
     * @tparam String a string type that can handle "addition" to `std::string`
     * @tparam Vector3D the type of 3D vector to be printed
     * @param s the string
     * @param manip the manipulator containing the vector to be printed
     * @return a new string with s concatenated to a rendering of the vector
     * @see VectorDumper
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * TVector3 pos;
     * std::runtime_error e("Position: " + gar::dump::vector3D(pos));
     * std::cout << e.what() << std::endl;
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * will produce an output like:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Position: { 0; 0; 0 }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * 
     * Requirements
     * -------------
     * 
     * * the type `Vector3D` is required to provide constant accessor methods
     *   returning the vector components, called `X()`, `Y()` and `Z()`
     * * the values returned by the aforementioned `Vector3D` methods must have
     *   their insertion operator into Stream defined
     * 
     */
    template <typename String, typename Vector>
    String operator+ (String const& s, VectorDumper<Vector> const& manip)
      { return s + std::string(manip);}
    
    /// Creates a string with s concatenated to the rendered vector.
    template <typename Vector>
    std::string operator+ (const char* s, VectorDumper<Vector> const& manip)
      { return std::string(s) + manip;}
    
    /// @see documentation of function with arguments swapped.
    template <typename String, typename Vector>
    String operator+ (VectorDumper<Vector> const& manip, String const& s)
      { return std::string(manip) + s;}
    
    /// Creates a string with the rendered vector concatenated to s.
    template <typename Vector>
    std::string operator+ (VectorDumper<Vector> const& manip, const char* s)
      { return manip + std::string(s);}
    
    /// Appends a string rendering of a vector to the specified string.
    template <typename String, typename Vector>
    String& operator+= (String& s, VectorDumper<Vector> const& manip)
      { return s += std::string(manip); }
    
  } // namespace dump
  
  
} // namespace gar


#endif // GARSOFT_COREUTILS_DUMPUTILS_H
