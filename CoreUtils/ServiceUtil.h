/**
 * @file   ServiceUtil.h
 * @brief  Utilities related to art service access
 * @author Jonathan Paley (jpaley@fnal.gov)
 *
 * This library is currently a pure header.
 *
 */

#ifndef COREUTILS_SERVICEUTIL
#define COREUTILS_SERVICEUTIL

// GArSoft libraries
#include "CoreUtils/UncopiableAndUnmovableClass.h"

// framework libraries
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/Exception.h"
#include "cetlib_except/demangle.h"

// C/C++ standard libraries
#include <type_traits>
#include <typeinfo>

namespace gar {

  namespace details {
    template <typename PROVIDER>
    struct ServiceRequirementsChecker;
  } // namespace details

  /** **************************************************************************
   * @brief Returns a constant pointer to the provider of specified service
   * @tparam T type of the service
   * @return a constant pointer to the provider of specified service
   * @throws ::art::Exception (category ::art::errors::NotFound) if pointer is null
   *
   * This function relies on the following service and provider interface:
   * - provider is not movable nor copiable
   * - service contains a type "provider_type" defined as the class of the
   *   service provider
   * - service contains a method "provider()" that returns a non-null pointer
   *   to a service provider; the service provider is owned and managed by
   *   the service, and the caller is not responsible of regulating the object
   *   lifetime, nor it should attempt to
   *
   * Violations of the protocol yield compilation errors (in case non-compliance
   * can be statically detected), or throw of exceptions.
   *
   * Example of usage:
   *
   *     auto const* geom = gar::providerFrom<geo::GeometryGAr>();
   *
   * retrieves the service provider for GArSoft geometry.
   * This requires the inclusion of "Geometry/GeometryGAr.h" header, where the
   * service is declared. Typically, both ServiceUtil.h and the header of the
   * provider class are included in the service header.
   *
   */
  template <class T>
  typename T::provider_type const* providerFrom()
  {
    details::ServiceRequirementsChecker<T>(); // instantiate a temporary...

    // retrieve the provider
    ::art::ServiceHandle<T> h;
    typename T::provider_type const* pProvider = h->provider();
    if (!pProvider) {
      throw ::art::Exception(::art::errors::NotFound)
        << "Service <" << cet::demangle_symbol(typeid(T).name()) << "> offered a null provider";
    }

    return pProvider;

  } // providerFrom()

  //----------------------------------------------------------------------------
  namespace details {
    /// Compiles only if PROVIDER class satisfied service provider requirements
    template <typename PROVIDER>
    struct ServiceProviderRequirementsChecker {

      using provider_type = PROVIDER;

      // static checks on provider class: not copiable nor movable
      static_assert(!std::is_copy_constructible<provider_type>::value,
                    "Service provider classes must not be copiable");
      static_assert(!std::is_copy_assignable<provider_type>::value,
                    "Service provider classes must not be copiable");
      static_assert(!std::is_move_constructible<provider_type>::value,
                    "Service provider classes must not be movable");
      static_assert(!std::is_move_assignable<provider_type>::value,
                    "Service provider classes must not be movable");

    }; // ServiceProviderRequirementsChecker

    template <typename SERVICE>
    struct ServiceRequirementsChecker {

      // require SERVICE::provider_type to be a type
      using provider_type = typename SERVICE::provider_type;

      // expected type for SERVICE::provider() method
      using provider_func_type = provider_type const* (SERVICE::*)() const;

      /// Checker for the provider
      ServiceProviderRequirementsChecker<provider_type> provider_checker;

      // check the provider() method
      static_assert(std::is_same<decltype(&SERVICE::provider), provider_func_type>::value,
                    "provider() method has unsupported signature");

    }; // ServiceRequirementsChecker

  } // namespace details

} // namespace gar

#endif //#COREUTILS_SERVICEUTIL
