#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit-tests for p3 ice table functions.
 */

template <typename D>
struct UnitWrap::UnitTest<D>::TestDsd2 {

  static void run_cloud_bfb()
  {
    // Read in tables
    view_2d_table vn_table_vals; view_2d_table vm_table_vals; view_2d_table revap_table_vals;
    view_1d_table mu_r_table_vals; view_dnu_table dnu;
    Functions::init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu);

    // Load some lookup inputs, need at least one per pack value
    GetCloudDsd2Data gcdd[max_pack_size] = {
      {0.998086E-06, 0.114746E+01, 0.100000E+01},
      {0.292481E-07, 0.758028E+00, 0.100000E+01},
      {0.510000E-02, 0.970269E+00, 0.100000E+01},
      {0.510000E-02, 0.980946E+00, 0.100000E+01},

      {0.998086E-06, 0.114746E+01, 0.100000E+01},
      {0.0,          0.758028E+00, 0.100000E+01},
      {0.510000E-02, 0.970269E+00, 0.100000E+01},
      {0.0,          0.980946E+00, 0.100000E+01},

      {0.000000E+00, 1.144941E+00, 1.000000E+00},
      {0.292481E-07, 0.758028E+00, 0.100000E+01},
      {0.510000E-02, 0.970269E+00, 0.100000E+01},
      {0.510000E-02, 0.980946E+00, 0.100000E+01},

      {0.998086E-06, 0.114746E+01, 0.100000E+01},
      {0.292481E-07, 0.758028E+00, 0.100000E+01},
      {0.510000E-02, 0.970269E+00, 0.100000E+01},
      {0.510000E-02, 0.980946E+00, 0.100000E+01}
    };

    // Sync to device
    view_1d<GetCloudDsd2Data> gcdd_device("gcdd", max_pack_size);
    const auto gcdd_host = Kokkos::create_mirror_view(gcdd_device);
    std::copy(&gcdd[0], &gcdd[0] + max_pack_size, gcdd_host.data());
    Kokkos::deep_copy(gcdd_device, gcdd_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      get_cloud_dsd2(gcdd[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qc, rho, nc;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qc[s]    = gcdd_device(vs).qc;
        rho[s]   = gcdd_device(vs).rho;
        nc[s]    = gcdd_device(vs).nc_in;
      }

      Spack mu_c(0.0), nu(0.0), lamc(0.0), cdist(0.0), cdist1(0.0);
      Functions::get_cloud_dsd2(qc, nc, mu_c, rho, nu, dnu, lamc, cdist, cdist1);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        gcdd_device(vs).nc_out = nc[s];
        gcdd_device(vs).mu_c   = mu_c[s];
        gcdd_device(vs).nu     = nu[s];
        gcdd_device(vs).lamc   = lamc[s];
        gcdd_device(vs).cdist  = cdist[s];
        gcdd_device(vs).cdist1 = cdist1[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(gcdd_host, gcdd_device);

    // Validate results
    if (SCREAM_BFB_TESTING) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(gcdd[s].nc_out == gcdd_host(s).nc_out);
        REQUIRE(gcdd[s].mu_c   == gcdd_host(s).mu_c);
        REQUIRE(gcdd[s].nu     == gcdd_host(s).nu);
        REQUIRE(gcdd[s].lamc   == gcdd_host(s).lamc);
        REQUIRE(gcdd[s].cdist  == gcdd_host(s).cdist);
        REQUIRE(gcdd[s].cdist1 == gcdd_host(s).cdist1);
      }
    }
  }

  static void run_cloud_phys()
  {
    // TODO
  }

  static void run_rain_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    GetRainDsd2Data grdd[max_pack_size] = {
      {0.100000E-01, 0.100000E+01},
      {0.100000E-01, 0.100000E+01},
      {0.156316E-03, 0.100000E+01},
      {0.148647E-03, 0.100000E+01},

      {0.100000E-01, 0.100000E+01},
      {0.0,          0.100000E+01},
      {0.156316E-03, 0.100000E+01},
      {0.0,          0.100000E+01},

      {0.100000E-01, 0.100000E+01},
      {0.100000E-01, 0.100000E+01},
      {0.156316E-03, 0.100000E+01},
      {0.148647E-03, 0.100000E+01},

      {0.100000E-01, 0.100000E+01},
      {0.100000E-01, 0.100000E+01},
      {0.156316E-03, 0.100000E+01},
      {0.148647E-03, 0.100000E+01}
    };

    // Sync to device
    KTH::view_1d<GetRainDsd2Data> grdd_host("grdd_host", max_pack_size);
    view_1d<GetRainDsd2Data> grdd_device("grdd_host", max_pack_size);
    std::copy(&grdd[0], &grdd[0] + max_pack_size, grdd_host.data());
    Kokkos::deep_copy(grdd_device, grdd_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      get_rain_dsd2(grdd[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack qr, cld_frac_r, nr;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        qr[s]    = grdd_device(vs).qr;
        nr[s]    = grdd_device(vs).nr_in;
      }

      Spack mu_r(0.0), lamr(0.0), cdistr(0.0), logn0r(0.0);
      Functions::get_rain_dsd2(qr, nr, mu_r, lamr,
                               p3::Functions<Real,DefaultDevice>::P3Runtime());
      Functions::get_cdistr_logn0r(qr, nr, mu_r, lamr, cdistr, logn0r);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        grdd_device(vs).nr_out = nr[s];
        grdd_device(vs).mu_r = mu_r[s];
        grdd_device(vs).lamr = lamr[s];
        grdd_device(vs).cdistr = cdistr[s];
        grdd_device(vs).logn0r = logn0r[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(grdd_host, grdd_device);

    // Validate results
    if (SCREAM_BFB_TESTING) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(grdd[s].nr_out == grdd_host(s).nr_out);
        REQUIRE(grdd[s].mu_r   == grdd_host(s).mu_r);
        REQUIRE(grdd[s].lamr   == grdd_host(s).lamr);
        REQUIRE(grdd[s].cdistr == grdd_host(s).cdistr);
        REQUIRE(grdd[s].logn0r == grdd_host(s).logn0r);
      }
    }
  }

  static void run_rain_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_cloud_dsd2", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDsd2;

  TD::run_cloud_phys();
  TD::run_cloud_bfb();
}

TEST_CASE("p3_rain_dsd2", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDsd2;

  TD::run_rain_phys();
  TD::run_rain_bfb();
}

}
