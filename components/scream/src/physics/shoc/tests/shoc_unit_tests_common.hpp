#ifndef SHOC_UNIT_TESTS_COMMON_HPP
#define SHOC_UNIT_TESTS_COMMON_HPP

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "physics/shoc/shoc_functions.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

/*
 * Unit test infrastructure for shoc unit tests.
 *
 * shoc entities can friend scream::shoc::unit_test::UnitWrap to give unit tests
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
    using uview_1d = typename ko::template Unmanaged<view_1d<S> >;

    using Functions          = scream::shoc::Functions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    using C                  = typename Functions::C;

    // Put struct decls here
    //struct TestStuff;
  };

};


} // namespace unit_test
} // namespace shoc
} // namespace scream

#endif
