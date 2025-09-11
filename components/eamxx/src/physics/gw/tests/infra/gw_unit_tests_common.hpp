#ifndef GW_UNIT_TESTS_COMMON_HPP
#define GW_UNIT_TESTS_COMMON_HPP

#include "gw_functions.hpp"
#include "gw_test_data.hpp"

#include "share/core/eamxx_types.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <vector>
#include <sstream>
#include <array>

namespace scream {
namespace gw {
namespace unit_test {

// Most bfb unit tests will use the same init
template <typename Engine>
inline auto get_common_init_data(Engine& engine)
{
  std::array<GwInit, 4> rv = {
    // gw_ediff::vd_lu_decomp breaks if kbot==pver

    // NOTE: All integer data is assumed to be 0-based (C style)! The
    // unit-test -> F90 GW layer needs to adjust these in d.transition if
    // it represents an index

    //     pver, pgwv,   dc, orog_only, molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv
    GwInit(  72,   20, 0.75,     false,      false,     false,         16,   8,     66,    .67, 6.28e-5),
    GwInit(  72,   20, 0.75,     true ,      false,     true ,         16,   6,     68,    .67, 6.28e-5),
    GwInit(  72,   20, 0.75,     false,      true ,     true ,         16,   3,     70,    .67, 6.28e-5),
    GwInit(  72,   20, 0.75,     true ,      true ,     false,         16,   0,     70,    .67, 6.28e-5),
  };

  for (auto& d : rv) {
    d.randomize(engine);
  }

  return rv;
}

/*
 * Unit test infrastructure for gw unit tests.
 *
 * gw entities can friend scream::gw::unit_test::UnitWrap to give unit tests
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

    using Functions          = scream::gw::Functions<Real, Device>;
    // using view_ice_table     = typename Functions::view_ice_table;
    // using view_collect_table = typename Functions::view_collect_table;
    // using view_1d_table      = typename Functions::view_1d_table;
    // using view_2d_table      = typename Functions::view_2d_table;
    // using view_dnu_table     = typename Functions::view_dnu_table;
    using Scalar             = typename Functions::Scalar;
    using Spack              = typename Functions::Spack;
    // using Pack               = typename Functions::Pack;
    // using IntSmallPack       = typename Functions::IntSmallPack;
    // using Smask              = typename Functions::Smask;
    // using TableIce           = typename Functions::TableIce;
    // using TableRain          = typename Functions::TableRain;
    // using Table3             = typename Functions::Table3;
    // using C                  = typename Functions::C;

    static constexpr Int max_pack_size = 16;
    static constexpr Int num_test_itrs = max_pack_size / Spack::n;

    struct Base : public UnitBase {

      Base() :
        UnitBase()
      {
        // Functions::gw_init(); // just in case there is ever global gw data
      }

      ~Base() = default;
    };

    // Put struct decls here
    struct TestGwdComputeTendenciesFromStressDivergence;
    struct TestGwProf;
    struct TestMomentumEnergyConservation;
    struct TestGwdComputeStressProfilesAndDiffusivities;
    struct TestGwdProjectTau;
    struct TestGwdPrecalcRhoi;
    struct TestGwDragProf;
    struct TestGwFrontProjectWinds;
    struct TestGwFrontGwSources;
    struct TestGwCmSrc;
    struct TestGwConvectProjectWinds;
    struct TestGwHeatingDepth;
    struct TestGwStormSpeed;
    struct TestGwConvectGwSources;
    struct TestGwBeresSrc;
    struct TestGwEdiff;
    struct TestGwDiffTend;
    struct TestGwOroSrc;
    struct TestVdLuDecomp;
    struct TestVdLuSolve;
  }; // UnitWrap
};

} // namespace unit_test
} // namespace gw
} // namespace scream

#endif
