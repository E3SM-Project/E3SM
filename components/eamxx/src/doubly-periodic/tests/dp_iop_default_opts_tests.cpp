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
struct UnitWrap::UnitTest<D>::TestIopDefaultOpts {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    IopDefaultOptsData f90_data[max_pack_size] = {
      // TODO
    };

    static constexpr Int num_runs = 0; //sizeof(f90_data) / sizeof(IopDefaultOptsData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    using host_view = typename view_1d<IopDefaultOptsData>::host_mirror_type;
    host_view cxx_host("cxx_host", max_pack_size);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());

    // Get data from fortran
    for (auto& d : f90_data) {
      iop_default_opts(d);
    }

    // Get data from cxx. Run iop_default_opts from a kernel and copy results back to host
    for (Int i = 0; i < num_test_itrs; ++i) {
      const Int offset = i * Spack::n;

      // Init outputs
      Spack iop_nudge_tq_high_out(0), iop_nudge_tq_low_out(0), iop_nudge_tscale_out(0), iop_perturb_high_out(0), scmlat_out(0), scmlon_out(0);

      Functions::iop_default_opts(scmlat_out, scmlon_out, cxx_host(0).iopfile_out, cxx_host(0).single_column_out, cxx_host(0).scm_iop_srf_prop_out, cxx_host(0).iop_nudge_tq_out, cxx_host(0).iop_nudge_uv_out, iop_nudge_tq_low_out, iop_nudge_tq_high_out, iop_nudge_tscale_out, cxx_host(0).scm_observed_aero_out, cxx_host(0).iop_dosubsidence_out, cxx_host(0).scm_multcols_out, cxx_host(0).dp_crm_out, iop_perturb_high_out, cxx_host(0).precip_off_out, cxx_host(0).scm_zero_non_iop_tracers_out);

      // Copy spacks back into cxx_host view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_host(vs).iop_nudge_tq_high_out = iop_nudge_tq_high_out[s];
        cxx_host(vs).iop_nudge_tq_low_out = iop_nudge_tq_low_out[s];
        cxx_host(vs).iop_nudge_tscale_out = iop_nudge_tscale_out[s];
        cxx_host(vs).iop_perturb_high_out = iop_perturb_high_out[s];
        cxx_host(vs).scmlat_out = scmlat_out[s];
        cxx_host(vs).scmlon_out = scmlon_out[s];
      }
    }

    // Verify BFB results
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        IopDefaultOptsData& d_f90 = f90_data[i];
        IopDefaultOptsData& d_cxx = cxx_host[i];
        REQUIRE(d_f90.scmlat_out == d_cxx.scmlat_out);
        REQUIRE(d_f90.scmlon_out == d_cxx.scmlon_out);
        REQUIRE(d_f90.iopfile_out == d_cxx.iopfile_out);
        REQUIRE(d_f90.single_column_out == d_cxx.single_column_out);
        REQUIRE(d_f90.scm_iop_srf_prop_out == d_cxx.scm_iop_srf_prop_out);
        REQUIRE(d_f90.iop_nudge_tq_out == d_cxx.iop_nudge_tq_out);
        REQUIRE(d_f90.iop_nudge_uv_out == d_cxx.iop_nudge_uv_out);
        REQUIRE(d_f90.iop_nudge_tq_low_out == d_cxx.iop_nudge_tq_low_out);
        REQUIRE(d_f90.iop_nudge_tq_high_out == d_cxx.iop_nudge_tq_high_out);
        REQUIRE(d_f90.iop_nudge_tscale_out == d_cxx.iop_nudge_tscale_out);
        REQUIRE(d_f90.scm_observed_aero_out == d_cxx.scm_observed_aero_out);
        REQUIRE(d_f90.iop_dosubsidence_out == d_cxx.iop_dosubsidence_out);
        REQUIRE(d_f90.scm_multcols_out == d_cxx.scm_multcols_out);
        REQUIRE(d_f90.dp_crm_out == d_cxx.dp_crm_out);
        REQUIRE(d_f90.iop_perturb_high_out == d_cxx.iop_perturb_high_out);
        REQUIRE(d_f90.precip_off_out == d_cxx.precip_off_out);
        REQUIRE(d_f90.scm_zero_non_iop_tracers_out == d_cxx.scm_zero_non_iop_tracers_out);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("iop_default_opts_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIopDefaultOpts;

  TestStruct::run_bfb();
}

} // empty namespace
