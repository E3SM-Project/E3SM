#ifndef P3_UNIT_TESTS_COMMON_HPP
#define P3_UNIT_TESTS_COMMON_HPP

#include "share/scream_types.hpp"
#include "p3_functions.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit test infrastructure for p3 unit tests.
 *
 * p3 entities can friend scream::p3::unit_test::UnitWrap to give unit tests
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

    using Functions          = scream::p3::Functions<Real, Device>;
    using view_ice_table     = typename Functions::view_ice_table;
    using view_collect_table = typename Functions::view_collect_table;
    using view_1d_table      = typename Functions::view_1d_table;
    using view_2d_table      = typename Functions::view_2d_table;
    using view_dnu_table     = typename Functions::view_dnu_table;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    using Pack               = typename Functions::Pack;
    using IntSmallPack       = typename Functions::IntSmallPack;
    using Smask              = typename Functions::Smask;
    using TableIce           = typename Functions::TableIce;
    using TableRain          = typename Functions::TableRain;
    using Table3             = typename Functions::Table3;
    using C                  = typename Functions::C;

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

    // Put struct decls here
    struct TestTableIce;
    struct TestTable3;
    struct TestBackToCellAverage;
    struct TestFind;
    struct TestUpwind;
    struct TestGenSed;
    struct TestP3Saturation;
    struct TestDsd2;
    struct TestP3Conservation;
    struct TestP3CloudWaterAutoconversion;
    struct TestCalcRimeDensity;
    struct TestCldliqImmersionFreezing;
    struct TestRainImmersionFreezing;
    struct TestDropletSelfCollection;
    struct TestCloudSed;
    struct TestCloudRainAccretion;
    struct TestIceSed;
    struct TestRainSed;
    struct TestGetTimeSpacePhysVariables;
    struct TestP3UpdatePrognosticIce;
    struct TestIceCollection;
    struct TestEvapSublPrecip;
    struct TestRainSelfCollection;
    struct TestP3IceMelting;
    struct TestP3SubgridVarianceScaling;
    struct TestP3UpdatePrognosticLiq;
    struct TestP3IceDepSublimation;
    struct TestP3FunctionsImposeMaxTotalNi;
    struct TestIceRelaxationTimescale;
    struct TestCalcLiqRelaxationTimescale;
    struct TestIceNucleation;
    struct TestIceCldliqWetGrowth;
    struct TestLatentHeat;
    struct TestCheckValues;
    struct TestIncloudMixing;
    struct TestP3Main;
    struct TestIceSupersatConservation;
    struct TestNcConservation;
    struct TestNrConservation;
    struct TestNiConservation;
    struct TestIceDepositionSublimation;
    struct TestPreventLiqSupersaturation;
  };

};

} // namespace unit_test
} // namespace p3
} // namespace scream

#endif
