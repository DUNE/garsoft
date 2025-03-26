/**
 * @file    ProviderPack.h
 * @brief   Data structure containing constant pointers to classes
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    November 20th, 2015
 * @version 1.0
 *
 * This is a header-only library.
 * It depends only on standard C++ and does not require additional linkage.
 */

#ifndef DETECTORINFO_PROVIDERPACK_H
#define DETECTORINFO_PROVIDERPACK_H

// C/C++ standard library
#include <limits>
#include <tuple>
#include <type_traits>

namespace gar {

  namespace details {

    template <typename... Types>
    struct has_duplicate_types;

    template <typename Key, typename... Types>
    struct index_with_type_impl;

    /**
     * @brief Hosts the index of type FindType in the list of types AmongTypes
     * @tparam FindType the type to find
     * @tparam AmongTypes the list of types
     *
     * If FindType type is not any of AmongTypes typesm a compilation error
     * message (static assertion) is issued.
     * Otherwise, the class has a `value` member pointing to the index of
     * FindType among AmongTypes.
     * This value is suct that this assertion is valid:
     *
     *     std::is_same<
     *       FindType,
     *       typename std::tuple_element<
     *         index_with_type<FindType, AmongTypes...>::value,
     *         std::tuple<AmongTypes...>
     *         >
     *       >::value
     *
     */
    template <typename FindType, typename... AmongTypes>
    struct index_with_type;

    /// Implementation detail for the extraction constructor
    template <typename DestPack, typename SourcePack, typename... ExtractProviders>
    struct SetFrom;

  } // namespace details

  /** **************************************************************************
   * @brief Container for a list of pointers to providers
   * @tparam Providers types of the providers in the parameter pack
   *
   * The pointers are stored as constant.
   * Note that this container can host any type of objects, and it has
   * "provider" in the name because the reason it was written was to provide
   * a fast way to specify a set of GArSoft service providers.
   * The only limitation is that there should be only one object per type.
   * Pointed objects are not owned by this class.
   *
   *     A a;
   *     B b;
   *     C c;
   *     D d;
   *     ProviderPack<A, B, C> pack(&a, &b, &c);
   *
   *     // obtain a constant pointerto b from pack:
   *     B const* b_ptr = pack.get<B>();
   *
   *     if (pack.has<D>()) std::cerr << "Unexpected!" << std::endl;
   *
   * (note that the latter check can be coded as a static assert)
   */
  template <typename... Providers>
  class ProviderPack {
    static_assert(!details::has_duplicate_types<Providers...>::value,
                  "Providers in ProviderPack are repeated");

    using this_type = ProviderPack<Providers...>; ///< alias of this class

    /// type used for storage of the pointers
    using tuple_type = std::tuple<Providers const*...>;

  public:
    /// Default constructor: a null provider pointer for each type
    ProviderPack() = default;

    /// Constructor: stores a provider pointer for each type
    ProviderPack(Providers const*... provider_ptrs) : providers(provider_ptrs...) {}

    /**
     * @brief Constructor: extracts the providers from anothe parameter pack
     * @tparam OtherProviders list of the providers of the source provider pack
     * @param from where to copy the information from
     *
     * This constructor requires all the providers we need to be present
     * in the source provider pack.
     */
    template <typename... OtherProviders>
    ProviderPack(ProviderPack<OtherProviders...> const& from)
    {
      details::SetFrom<this_type, ProviderPack<OtherProviders...>, Providers...>(*this, from);
    }

    /**
     * @brief Constructor: extracts the providers from anothe parameter pack
     * @tparam OtherProviders list of the providers of the source provider pack
     * @param providers all the providers needed (or more)
     *
     * This constructor requires all the providers we need to be present
     * in the source provider pack.
     */
    template <typename... OtherProviders>
    ProviderPack(OtherProviders const*... providers_tmp)
    {
      details::SetFrom<this_type, ProviderPack<OtherProviders...>, Providers...>(
        *this, ProviderPack<OtherProviders...>(providers_tmp...));
    }

    /// Returns the provider with the specified type
    template <typename Provider>
    Provider const* get() const
    {
      return std::get<details::index_with_type<Provider, Providers...>::value>(providers);
    } // get<>()

    /// Sets the provider with the specified type
    template <typename Provider>
    void set(Provider const* provider_ptr)
    {
      std::get<details::index_with_type<Provider, Providers...>::value>(providers) = provider_ptr;
    } // set<>()

    /// Returns whether there is a provider with the specified type
    template <typename Provider>
    static constexpr bool has()
    {
      return details::index_with_type_impl<Provider, Providers...>::value < sizeof...(Providers);
    } // has<>()

  private:
    tuple_type providers; ///< container of the pointers, type-safe

  }; // class ProviderPack

  /**
   * @brief Function to create a ParameterPack from the function arguments
   * @tparam Providers types of the providers in the parameter pack
   * @param providers constant pointers to the providers
   * @return a ParameterPack object containing all the specified providers
   *
   * This is an convevience function to reduce the typing needed to instantiate
   * a ParameterPack.
   */
  template <typename... Providers>
  ProviderPack<Providers...> makeProviderPack(Providers const*... providers)
  {
    return ProviderPack<Providers...>(providers...);
  }

} // namespace lar

//------------------------------------------------------------------------------
//--- Implementation details
//---
namespace gar {
  namespace details {

    //--------------------------------------------------------------------------
    //--- has_duplicate_types
    //---
    template <typename Key, typename... Types>
      struct has_duplicate_types<Key, Types...> : public std::integral_constant < bool,
      index_with_type_impl<Key, Types...>::value<sizeof...(Types) ||
                                                 has_duplicate_types<Types...>::value> {};

    template <>
    struct has_duplicate_types<> : public std::false_type {};

    //--------------------------------------------------------------------------
    //--- index_with_type
    //---
    /// Generic template, to glue the two special cases;
    /// in general, this will have in value field the index in Types of the type
    /// Key, or a number larger than the number if types in Types if Key type
    /// is not present among Types.
    template <typename Key, typename... Types>
    struct index_with_type_impl;

    // common case: one or more types
    template <typename Key, typename FirstType, typename... Types>
    struct index_with_type_impl<Key, FirstType, Types...>
      : public std::integral_constant<size_t,
                                      std::is_same<Key, FirstType>::value ?
                                        0 :
                                        index_with_type_impl<Key, Types...>::value + 1> {};

    // special case: no type; this means failure
    template <typename Key>
    struct index_with_type_impl<Key> : public std::integral_constant<size_t, 0U> {};

    template <typename FindType, typename... AmongTypes>
    struct index_with_type
      : public std::integral_constant<size_t,
                                      index_with_type_impl<FindType, AmongTypes...>::value> {
      static_assert(index_with_type<FindType, AmongTypes...>::value < sizeof...(AmongTypes),
                    "Required type is not present");
    }; // index_with_type

    //--------------------------------------------------------------------------
    //--- SetFrom
    //---
    template <typename DestPack,
              typename SourcePack,
              typename FirstProvider,
              typename... OtherProviders>
    struct SetFrom<DestPack, SourcePack, FirstProvider, OtherProviders...> {
      SetFrom(DestPack& pack, SourcePack const& from)
      {
        pack.set(from.template get<FirstProvider>());
        SetFrom<DestPack, SourcePack, OtherProviders...>(pack, from);
      }
    }; // SetFrom<First, Others...>

    template <typename DestPack, typename SourcePack>
    struct SetFrom<DestPack, SourcePack> {
      SetFrom(DestPack&, SourcePack const&) {}
    };

    //--------------------------------------------------------------------------

  } // namespace details

  //----------------------------------------------------------------------------
  //--- ProviderPack
  //---

  //----------------------------------------------------------------------------

} // namespace gar

#endif // DETECTORINFO_PROVIDERPACK_H
