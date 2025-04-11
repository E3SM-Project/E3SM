#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit-tests for p3 ice collection functions.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestRainSelfCollection : public UnitWrap::UnitTest<D>::Base {

  void run_rain_self_collection_bfb_tests() {

    RainSelfCollectionData dc[max_pack_size] = {
      //  rho, qr_incld, nr_incld, nr_selfcollect_tend
      {1.060E+00, 1.354E-03, 1.401E+04, 0.000E+00},
      {1.136E+00, 1.319E-04, 4.210E+03, 0.000E+00},
      {1.166E+00, 5.339E-03, 1.804E+04, 0.000E+00},
      {1.120E+00, 5.650E-05, 8.277E+05, 0.000E+00},

      {1.060E+00, 1.354E-03, 1.401E+04, 0.000E+00},
      {1.136E+00, 1.319E-04, 4.210E+03, 0.000E+00},
      {1.166E+00, 5.339E-03, 1.804E+04, 0.000E+00},
      {1.120E+00, 5.650E-05, 8.277E+05, 0.000E+00},

      {1.060E+00, 1.354E-03, 1.401E+04, 0.000E+00},
      {1.136E+00, 1.319E-04, 4.210E+03, 0.000E+00},
      {1.166E+00, 5.339E-03, 1.804E+04, 0.000E+00},
      {1.120E+00, 5.650E-05, 8.277E+05, 0.000E+00},

      {1.060E+00, 1.354E-03, 1.401E+04, 0.000E+00},
      {1.136E+00, 1.319E-04, 4.210E+03, 0.000E+00},
      {1.166E+00, 5.339E-03, 1.804E+04, 0.000E+00},
      {1.120E+00, 5.650E-05, 8.277E+05, 0.000E+00}
    };

    //Sync to device
    view_1d<RainSelfCollectionData> dc_device("dc", max_pack_size);
    auto dc_host = Kokkos::create_mirror_view(dc_device);

    std::copy(&dc[0], &dc[0] + max_pack_size, dc_host.data());
    Kokkos::deep_copy(dc_device, dc_host);

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        dc[i].read(Base::m_fid);
      }
    }

    //Run function from a kernal and copy results back to the host
    Kokkos::parallel_for(RangePolicy(0, num_test_itrs), KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho_local, qr_incld_local, nr_incld_local, nr_selfcollect_tend_local;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho_local[s]      = dc_device(vs).rho;
        qr_incld_local[s] = dc_device(vs).qr_incld;
        nr_incld_local[s] = dc_device(vs).nr_incld;
        nr_selfcollect_tend_local[s]    = dc_device(vs).nr_selfcollect_tend;
      }

      Functions::rain_self_collection(
          rho_local, qr_incld_local, nr_incld_local, nr_selfcollect_tend_local,
          p3::Functions<Real,DefaultDevice>::P3Runtime());

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        dc_device(vs).rho                  = rho_local[s];
        dc_device(vs).qr_incld             = qr_incld_local[s];
        dc_device(vs).nr_incld             = nr_incld_local[s];
        dc_device(vs).nr_selfcollect_tend  = nr_selfcollect_tend_local[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(dc_host, dc_device);

    // Validate results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(dc[s].rho                 == dc_host(s).rho);
        REQUIRE(dc[s].qr_incld            == dc_host(s).qr_incld);
        REQUIRE(dc[s].nr_incld            == dc_host(s).nr_incld);
        REQUIRE(dc[s].nr_selfcollect_tend == dc_host(s).nr_selfcollect_tend);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        dc_host(s).write(Base::m_fid);
      }
    }
  }

  void run_bfb() {
    run_rain_self_collection_bfb_tests();
  }

}; //TestRainselfCollection

} // namespace unit_test
} // p3
} // scream

namespace {

TEST_CASE("p3_rain_self_collection_test", "[p3_rain_self_collection_test"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainSelfCollection;

  T t;
  t.run_bfb();
}

} // namespace
