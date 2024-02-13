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
struct UnitWrap::UnitTest<D>::TestIopSetopts {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    IopSetoptsData f90_data[max_pack_size] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(IopSetoptsData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    using host_view = typename view_1d<IopSetoptsData>::host_mirror_type;
    host_view cxx_host("cxx_host", max_pack_size);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());

    // Get data from fortran
    for (auto& d : f90_data) {
      iop_setopts(d);
    }

    // Get data from cxx. Run iop_setopts from a kernel and copy results back to host
    for (Int i = 0; i < num_test_itrs; ++i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack iop_nudge_tq_high_in, iop_nudge_tq_low_in, iop_nudge_tscale_in, iop_perturb_high_in, scmlat_in, scmlon_in;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        iop_nudge_tq_high_in[s] = cxx_host(vs).iop_nudge_tq_high_in;
        iop_nudge_tq_low_in[s] = cxx_host(vs).iop_nudge_tq_low_in;
        iop_nudge_tscale_in[s] = cxx_host(vs).iop_nudge_tscale_in;
        iop_perturb_high_in[s] = cxx_host(vs).iop_perturb_high_in;
        scmlat_in[s] = cxx_host(vs).scmlat_in;
        scmlon_in[s] = cxx_host(vs).scmlon_in;
      }

      Functions::iop_setopts(scmlat_in, scmlon_in, cxx_host(0).iopfile_in, cxx_host(0).single_column_in, cxx_host(0).scm_iop_srf_prop_in, cxx_host(0).iop_nudge_tq_in, cxx_host(0).iop_nudge_uv_in, iop_nudge_tq_low_in, iop_nudge_tq_high_in, iop_nudge_tscale_in, cxx_host(0).scm_observed_aero_in, cxx_host(0).iop_dosubsidence_in, cxx_host(0).scm_multcols_in, cxx_host(0).dp_crm_in, iop_perturb_high_in, cxx_host(0).precip_off_in, cxx_host(0).scm_zero_non_iop_tracers_in);
    }

    // Verify BFB results
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        IopSetoptsData& d_f90 = f90_data[i];
        IopSetoptsData& d_cxx = cxx_host[i];
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("iop_setopts_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIopSetopts;

  TestStruct::run_bfb();
}

} // empty namespace
