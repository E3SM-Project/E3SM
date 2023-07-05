#ifndef DP_UNIT_TESTS_COMMON_HPP
#define DP_UNIT_TESTS_COMMON_HPP

#include "dp_functions.hpp"
#include "share/scream_types.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {
namespace dp {
namespace unit_test {

/*
 * Unit test infrastructure for dp unit tests.
 *
 * dp entities can friend scream::dp::unit_test::UnitWrap to give unit tests
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

    using Functions          = scream::dp::Functions<Real, Device>;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    using Pack               = typename Functions::Pack;
    using IntSmallPack       = typename Functions::IntSmallPack;
    using Smask              = typename Functions::Smask;
    using C                  = typename Functions::C;

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

    // Put struct decls here
    struct TestAdvanceIopForcing;
    struct TestAdvanceIopNudging;
    struct TestAdvanceIopSubsidence;
    struct TestIopSetinitial;
    struct TestIopBroadcast;
    struct TestApplyIopForcing;
    struct TestIopDomainRelaxation;
    struct TestCrmResolvedTurb;
    struct TestIopDefaultOpts;
    struct TestIopSetopts;
    struct TestSetiopupdateInit;
    struct TestSetiopupdate;
    struct TestReadiopdata;
    struct TestIopIntht;
  };

};


} // namespace unit_test
} // namespace dp
} // namespace scream

#endif
