#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestNcConservation : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    NcConservationData baseline_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
      d.dt = baseline_data[0].dt; // Hold this fixed, this is not packed data
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before reads so that
    // inout data is in original state
    view_1d<NcConservationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&baseline_data[0], &baseline_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        baseline_data[i].read(Base::m_fid);
      }
    }

    // Get data from cxx. Run nc_conservation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack nc, nc2ni_immers_freeze_tend, nc2nr_autoconv_tend, nc_accret_tend, nc_collect_tend, nc_selfcollect_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        nc[s] = cxx_device(vs).nc;
        nc2ni_immers_freeze_tend[s] = cxx_device(vs).nc2ni_immers_freeze_tend;
        nc2nr_autoconv_tend[s] = cxx_device(vs).nc2nr_autoconv_tend;
        nc_accret_tend[s] = cxx_device(vs).nc_accret_tend;
        nc_collect_tend[s] = cxx_device(vs).nc_collect_tend;
        nc_selfcollect_tend[s] = cxx_device(vs).nc_selfcollect_tend;
      }

      Functions::nc_conservation(nc, nc_selfcollect_tend, cxx_device(offset).dt, nc_collect_tend, nc2ni_immers_freeze_tend, nc_accret_tend, nc2nr_autoconv_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).nc2ni_immers_freeze_tend = nc2ni_immers_freeze_tend[s];
        cxx_device(vs).nc2nr_autoconv_tend = nc2nr_autoconv_tend[s];
        cxx_device(vs).nc_accret_tend = nc_accret_tend[s];
        cxx_device(vs).nc_collect_tend = nc_collect_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        NcConservationData& d_baseline = baseline_data[i];
        NcConservationData& d_cxx = cxx_host[i];
        REQUIRE(d_baseline.nc_collect_tend == d_cxx.nc_collect_tend);
        REQUIRE(d_baseline.nc2ni_immers_freeze_tend == d_cxx.nc2ni_immers_freeze_tend);
        REQUIRE(d_baseline.nc_accret_tend == d_cxx.nc_accret_tend);
        REQUIRE(d_baseline.nc2nr_autoconv_tend == d_cxx.nc2nr_autoconv_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        cxx_host(s).write(Base::m_fid);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("nc_conservation_bfb", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNcConservation;

  T t;
  t.run_bfb();
}

} // empty namespace
