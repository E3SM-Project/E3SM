#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "doubly-periodic/dp_functions.hpp"
#include "doubly-periodic/dp_functions_f90.hpp"

#include "dp_unit_tests_common.hpp"

namespace scream {
namespace dp {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestIopSetinitial {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    IopSetinitialData f90_data[] = {
      //                plev, pcnst, nelemd, np, nstep, psobs, use_replay, dynproc, have_t, have_q, have_ps, have_u, have_v, have_numliq, have_cldliq, have_numice, have_cldice, scm_zero_non_iop_tracers, is_first_restart_step
      IopSetinitialData(72  , 5    , 10    , 2 , 10,    0.1,   true      , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 10,    0.1,   false     , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 1 ,    0.1,   true      , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 1 ,    0.1,   false     , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 1 ,    0.1,   true      , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , false),
      IopSetinitialData(72  , 5    , 10    , 2 , 1 ,    0.1,   false     , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , false),
      IopSetinitialData(72  , 5    , 10    , 2 , 0 ,    0.1,   true      , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 0 ,    0.1,   false     , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , true                    , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 0 ,    0.1,   true      , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , false                   , true),
      IopSetinitialData(72  , 5    , 10    , 2 , 0 ,    0.1,   false     , true   , true  , true  , true   , true  , true  , true       , true       , true       , true       , false                   , true),
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(IopSetinitialData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    IopSetinitialData cxx_data[] = {
      IopSetinitialData(f90_data[0]),
      IopSetinitialData(f90_data[1]),
      IopSetinitialData(f90_data[2]),
      IopSetinitialData(f90_data[3]),
      IopSetinitialData(f90_data[4]),
      IopSetinitialData(f90_data[5]),
      IopSetinitialData(f90_data[6]),
      IopSetinitialData(f90_data[7]),
      IopSetinitialData(f90_data[8]),
      IopSetinitialData(f90_data[9]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      iop_setinitial(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.init();
      iop_setinitial_f(d.plev, d.pcnst, d.nelemd, d.np, d.nstep, d.psobs, d.use_replay, d.dynproc, d.have_t, d.have_q, d.have_ps, d.have_u, d.have_v, d.have_numliq, d.have_cldliq, d.have_numice, d.have_cldice, d.scm_zero_non_iop_tracers, d.is_first_restart_step, d.qmin, d.uobs, d.vobs, d.numliqobs, d.numiceobs, d.cldliqobs, d.cldiceobs, d.dx_short, &d.tracers, &d.elem, &d.dyn_dx_size, d.tobs, d.qobs);
    }

    // We can't call into fortran. Due to all the dependencies it has, it's not possible
    // to build it in standalone eamxx. Without fortran, we cannot do BFB tests.
#if 0
    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        IopSetinitialData& d_f90 = f90_data[i];
        IopSetinitialData& d_cxx = cxx_data[i];

      }
    }
#endif
  } // run_bfb

};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("iop_setinitial_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIopSetinitial;

  TestStruct::run_bfb();
}

} // empty namespace
