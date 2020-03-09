#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

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
struct UnitWrap::UnitTest<D>::TestRainSelfCollection {

  static void run_rain_self_collection_bfb_tests(){

    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);
    
    RainSelfCollectionData dc[max_pack_size] = {
      //  rho, qr_incld, nr_incld, nrslf
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
    KTH::view_1d<RainSelfCollectionData> dc_device("dc", Spack::n);
    auto dc_host = Kokkos::create_mirror_view(dc_device);

    std::copy(&dc[0], &dc[0] + Spack::n, dc_host.data());
    Kokkos::deep_copy(dc_device, dc_host);

    //Get data from fortran
  	for (Int i = 0; i < Spack::n; ++i) {
      rain_self_collection(dc[i]);
  	}

    //Run function from a kernal and copy results back to the host
  	Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      // Init pack inputs
      Spack rho_local, qr_incld_local, nr_incld_local, nrslf_local;
      for (Int s = 0; s < Spack::n; ++s) {
        rho_local[s] = dc_device(s).rho;
        qr_incld_local[s] = dc_device(s).qr_incld;
        nr_incld_local[s] = dc_device(s).nr_incld;
        nrslf_local[s] = dc_device(s).nrslf;
       }

      Functions::rain_self_collection(rho_local, qr_incld_local, nr_incld_local, nrslf_local);
      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        dc_device(s).rho = rho_local[s];
        dc_device(s).qr_incld = qr_incld_local[s];
        dc_device(s).nr_incld = nr_incld_local[s];
        dc_device(s).nrslf = nrslf_local[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(dc_host, dc_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(dc[s].rho == dc_host(s).rho);
      REQUIRE(dc[s].qr_incld == dc_host(s).qr_incld);
      REQUIRE(dc[s].nr_incld == dc_host(s).nr_incld);
      REQUIRE(dc[s].nrslf == dc_host(s).nrslf);
    }
  }

  static void run_bfb(){
    run_rain_self_collection_bfb_tests();
  }

}; //TestRainselfCollection


} // namespace unit_test
} // p3
} // scream

namespace{
  TEST_CASE("p3_rain_self_collection_test", "[p3_rain_self_collection_test"){
    scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRainSelfCollection::run_bfb();
  }
} // namespace