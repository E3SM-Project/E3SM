#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestNrConservation {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    NrConservationData f90_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
      d.dt = f90_data[0].dt; // hold dt fixed, it is not packed data
      d.nmltratio = f90_data[0].nmltratio; // hold nmltratio fixed, it is not packed data
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    view_1d<NrConservationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Get data from fortran
    for (auto& d : f90_data) {
      nr_conservation(d);
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
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < max_pack_size; ++i) {
        NrConservationData& d_f90 = f90_data[i];
        NrConservationData& d_cxx = cxx_host[i];
        REQUIRE(d_f90.nr_collect_tend == d_cxx.nr_collect_tend);
        REQUIRE(d_f90.nr2ni_immers_freeze_tend == d_cxx.nr2ni_immers_freeze_tend);
        REQUIRE(d_f90.nr_selfcollect_tend == d_cxx.nr_selfcollect_tend);
        REQUIRE(d_f90.nr_evap_tend == d_cxx.nr_evap_tend);
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
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNrConservation;

  TestStruct::run_bfb();
}

} // empty namespace
