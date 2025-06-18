#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

// #include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZMComputeCAPE : public UnitWrap::UnitTest<D>::Base {

  void run_zm_compute_cape()
  {
    auto engine = Base::get_engine();

    // // Set up init data
    // GwInit init_data[] = {
    //       // pver, pgwv,   dc, orog_only, molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv
    //   GwInit(  72,   20, 0.75,     false,      false,     false,         16,   60,     16,    .67, 6.28e-5),
    //   GwInit(  72,   20, 0.75,     true ,      false,     true ,         16,   60,     16,    .67, 6.28e-5),
    //   GwInit(  72,   20, 0.75,     false,      true ,     true ,         16,   60,     16,    .67, 6.28e-5),
    //   GwInit(  72,   20, 0.75,     true ,      true ,     false,         16,   60,     16,    .67, 6.28e-5),
    // };

    // for (auto& d : init_data) {
    //   d.randomize(engine);
    // }

    // // Set up inputs
    // GwdComputeTendenciesFromStressDivergenceData baseline_data[] = {
    //   //                                        ncol, ngwv, do_taper,   dt, effgw, init
    //   GwdComputeTendenciesFromStressDivergenceData(2,   10,    false,  0.4,   0.3, init_data[0]),
    //   GwdComputeTendenciesFromStressDivergenceData(3,   10,    false,  0.4,   0.3, init_data[0]),
    //   GwdComputeTendenciesFromStressDivergenceData(4,   10,    true ,  0.4,   0.3, init_data[0]),
    //   GwdComputeTendenciesFromStressDivergenceData(5,   10,    true ,  0.4,   0.3, init_data[0]),
    // };
    
    // Get data from test
    for (auto& d : test_data) {
      gwd_compute_tendencies_from_stress_divergence(d);
    }

    // Verify results, all data should be in C layout
    // ???

  } // run_zm_compute_cape

};

} // namespace unit_test
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("zm_compute_cape", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZMComputeCAPE;

  TestStruct t;
  t.run_zm_compute_cape();
}

} // empty namespace
