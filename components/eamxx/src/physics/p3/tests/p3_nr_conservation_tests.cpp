#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestNrConservation : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    NrConservationData baseline_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
      d.dt = baseline_data[0].dt; // hold dt fixed, it is not packed data
      d.nmltratio = baseline_data[0].nmltratio; // hold nmltratio fixed, it is not packed data
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before reads so that
    // inout data is in original state
    view_1d<NrConservationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&baseline_data[0], &baseline_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        baseline_data[i].read(Base::m_fid);
      }
    }

    // Get data from cxx. Run nr_conservation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack nc2nr_autoconv_tend, ncshdc, ni2nr_melt_tend, nr, nr2ni_immers_freeze_tend, nr_collect_tend, nr_evap_tend, nr_ice_shed_tend, nr_selfcollect_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        nc2nr_autoconv_tend[s] = cxx_device(vs).nc2nr_autoconv_tend;
        ncshdc[s] = cxx_device(vs).ncshdc;
        ni2nr_melt_tend[s] = cxx_device(vs).ni2nr_melt_tend;
        nr[s] = cxx_device(vs).nr;
        nr2ni_immers_freeze_tend[s] = cxx_device(vs).nr2ni_immers_freeze_tend;
        nr_collect_tend[s] = cxx_device(vs).nr_collect_tend;
        nr_evap_tend[s] = cxx_device(vs).nr_evap_tend;
        nr_ice_shed_tend[s] = cxx_device(vs).nr_ice_shed_tend;
        nr_selfcollect_tend[s] = cxx_device(vs).nr_selfcollect_tend;
      }

      Functions::nr_conservation(nr, ni2nr_melt_tend, nr_ice_shed_tend, ncshdc, nc2nr_autoconv_tend, cxx_device(offset).dt, cxx_device(offset).nmltratio, nr_collect_tend, nr2ni_immers_freeze_tend, nr_selfcollect_tend, nr_evap_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).nr2ni_immers_freeze_tend = nr2ni_immers_freeze_tend[s];
        cxx_device(vs).nr_collect_tend = nr_collect_tend[s];
        cxx_device(vs).nr_evap_tend = nr_evap_tend[s];
        cxx_device(vs).nr_selfcollect_tend = nr_selfcollect_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        NrConservationData& d_baseline = baseline_data[i];
        NrConservationData& d_cxx = cxx_host[i];
        REQUIRE(d_baseline.nr_collect_tend == d_cxx.nr_collect_tend);
        REQUIRE(d_baseline.nr2ni_immers_freeze_tend == d_cxx.nr2ni_immers_freeze_tend);
        REQUIRE(d_baseline.nr_selfcollect_tend == d_cxx.nr_selfcollect_tend);
        REQUIRE(d_baseline.nr_evap_tend == d_cxx.nr_evap_tend);
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

TEST_CASE("nr_conservation_bfb", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNrConservation;

  T t;
  t.run_bfb();
}

} // empty namespace
