#ifndef ZM_UNIT_TESTS_COMMON_HPP
#define ZM_UNIT_TESTS_COMMON_HPP

#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "zm_functions.hpp"
#include "zm_test_data.hpp"

#include <ekat_test_utils.hpp>

#include <vector>
#include <sstream>

namespace scream {
namespace zm {
namespace unit_test {

/*
 * Unit test infrastructure for zm unit tests.
 *
 * zm entities can friend scream::zm::unit_test::UnitWrap to give unit tests
 * access to private members.
 *
 * All unit test impls should be within an inner struct of UnitWrap::UnitTest for
 * easy access to useful types.
 */

struct UnitWrap {

  template <typename D=DefaultDevice>
  struct UnitTest : public KokkosTypes<D> {

    using Device      = D;
    using MemberType  = typename KokkosTypes<Device>::MemberType;
    using TeamPolicy  = typename KokkosTypes<Device>::TeamPolicy;
    using RangePolicy = typename KokkosTypes<Device>::RangePolicy;
    using ExeSpace    = typename KokkosTypes<Device>::ExeSpace;

    template <typename S>
    using view_1d = typename KokkosTypes<Device>::template view_1d<S>;
    template <typename S>
    using view_2d = typename KokkosTypes<Device>::template view_2d<S>;
    template <typename S>
    using view_3d = typename KokkosTypes<Device>::template view_3d<S>;

    template <typename S>
    using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

    using Functions          = scream::zm::Functions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

    struct Base : public UnitBase {

      Base() :
        UnitBase()
      {
        // Functions::zm_init(); // we might need this if there is ever global zm data
      }

      ~Base() = default;
    };

    // Put struct decls here
    struct Test_zm_find_mse_max;
    
  }; // UnitWrap
};

} // namespace unit_test
} // namespace zm
} // namespace scream

#endif
