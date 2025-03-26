/**
 * @file   UncopiableAndUnmovableClass.h
 * @brief  Defines a class that can't be copied nor moved.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 *
 * This library is currently a pure header.
 *
 */

#ifndef COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H
#define COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H 1

namespace gar {

  /** **************************************************************************
   * @brief An empty class that can't be copied nor moved
   *
   * A class derived from this one can still be copied (or moved)
   * with an explicit effort. For example, to enable copy construction:
   *
   *     struct CopiableClass: protected UncopiableAndUnmovableClass {
   *       CopiableClass(CopiableClass const& from)
   *         : UncopiableAndUnmovableClass() // , ...
   *         {
   *           // ...
   *         }
   *     };
   *
   * the default constructor of the base class can be called explicitly instead
   * of the copy constructor. To provide an assignment operation,
   *
   *     struct MoveAssignableClass: protected UncopiableAndUnmovableClass {
   *       MoveAssignableClass& operator= (MoveAssignableClass&& from)
   *         {
   *           // ...
   *           return *this;
   *         }
   *     };
   *
   *
   */
  struct UncopiableAndUnmovableClass {

    /// Default constructor
    UncopiableAndUnmovableClass() = default;

    // @{
    /// Deleted copy and move constructors and assignments
    UncopiableAndUnmovableClass(UncopiableAndUnmovableClass const&) = delete;
    UncopiableAndUnmovableClass(UncopiableAndUnmovableClass&&) = delete;

    UncopiableAndUnmovableClass& operator=(UncopiableAndUnmovableClass const&) = delete;
    UncopiableAndUnmovableClass& operator=(UncopiableAndUnmovableClass&&) = delete;
    // @}

    /// Default destructor
    ~UncopiableAndUnmovableClass() = default;

  }; // UncopiableAndUnmovableClass

} // namespace gar

#endif // COREUTILS_UNCOPIABLEANDUNMOVEABLECLASS_H
