#ifndef TMS_UNIT_TESTS_COMMON_HPP
#define TMS_UNIT_TESTS_COMMON_HPP

#include "tms_functions.hpp"
#include "share/eamxx_types.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

namespace scream {
namespace tms {
namespace unit_test {

/*
 * Unit test infrastructure for tms unit tests.
 *
 * tms entities can friend scream::tms::unit_test::UnitWrap to give unit tests
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

    using Functions          = scream::tms::Functions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;

    // Put struct decls here
    struct TestComputeTMS;
  };

};


} // namespace unit_test
} // namespace tms
} // namespace scream

#endif
