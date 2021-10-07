#ifndef SPA_UNIT_TESTS_COMMON_HPP
#define SPA_UNIT_TESTS_COMMON_HPP

#include "share/scream_types.hpp"

#include "physics/spa/spa_functions.hpp"

#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {
namespace spa {
namespace unit_test {

/*
 * Unit test infrastructure for spa unit tests.
 *
 * spa entities can friend scream::spa::unit_test::UnitWrap to give unit tests
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

    using Functions          = scream::spa::SPAFunctions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    using Pack               = typename Functions::Pack;

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

    // Put struct decls here
    struct TestReadRemapData;
    struct TestReadDataFile;

  };

};


} // namespace unit_test
} // namespace shoc
} // namespace scream

#endif
