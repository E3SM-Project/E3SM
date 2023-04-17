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
struct UnitWrap::UnitTest<D>::TestNiConservation {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    NiConservationData f90_data[max_pack_size];

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
      d.dt = f90_data[0].dt; // hold dt fixed, it is not packed data
    }

    // Create copies of data for use by cxx and sync it to device. Needs to happen before fortran calls so that
    // inout data is in original state
    view_1d<NiConservationData> cxx_device("cxx_device", max_pack_size);
    const auto cxx_host = Kokkos::create_mirror_view(cxx_device);
    std::copy(&f90_data[0], &f90_data[0] + max_pack_size, cxx_host.data());
    Kokkos::deep_copy(cxx_device, cxx_host);

    // Get data from fortran
    for (auto& d : f90_data) {
      ni_conservation(d);
    }

    // Get data from cxx. Run ni_conservation from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack nc2ni_immers_freeze_tend, ni, ni2nr_melt_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend, nr2ni_immers_freeze_tend;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        nc2ni_immers_freeze_tend[s] = cxx_device(vs).nc2ni_immers_freeze_tend;
        ni[s] = cxx_device(vs).ni;
        ni2nr_melt_tend[s] = cxx_device(vs).ni2nr_melt_tend;
        ni_nucleat_tend[s] = cxx_device(vs).ni_nucleat_tend;
        ni_selfcollect_tend[s] = cxx_device(vs).ni_selfcollect_tend;
        ni_sublim_tend[s] = cxx_device(vs).ni_sublim_tend;
        nr2ni_immers_freeze_tend[s] = cxx_device(vs).nr2ni_immers_freeze_tend;
      }

      Functions::ni_conservation(ni, ni_nucleat_tend, nr2ni_immers_freeze_tend, nc2ni_immers_freeze_tend, cxx_device(offset).dt, ni2nr_melt_tend, ni_sublim_tend, ni_selfcollect_tend);

      // Copy spacks back into cxx_device view
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cxx_device(vs).ni2nr_melt_tend = ni2nr_melt_tend[s];
        cxx_device(vs).ni_selfcollect_tend = ni_selfcollect_tend[s];
        cxx_device(vs).ni_sublim_tend = ni_sublim_tend[s];
      }

    });

    Kokkos::deep_copy(cxx_host, cxx_device);

    // Verify BFB results
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < max_pack_size; ++i) {
        NiConservationData& d_f90 = f90_data[i];
        NiConservationData& d_cxx = cxx_host[i];
        REQUIRE(d_f90.ni2nr_melt_tend == d_cxx.ni2nr_melt_tend);
        REQUIRE(d_f90.ni_sublim_tend == d_cxx.ni_sublim_tend);
        REQUIRE(d_f90.ni_selfcollect_tend == d_cxx.ni_selfcollect_tend);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("ni_conservation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNiConservation;

  TestStruct::run_bfb();
}

} // empty namespace
