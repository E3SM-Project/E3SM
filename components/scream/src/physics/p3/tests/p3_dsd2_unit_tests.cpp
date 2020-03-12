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
    view_2d_table vn_table; view_2d_table vm_table; view_2d_table revap_table;
    view_1d_table mu_r_table; view_dnu_table dnu;
    Functions::init_kokkos_tables(vn_table, vm_table, revap_table, mu_r_table, dnu);

    constexpr Scalar qsmall = C::QSMALL;
    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    // Load some lookup inputs, need at least one per pack value
    GetCloudDsd2Data gcdd[max_pack_size] = {
      {0.998086E-06, 0.114746E+01, 0.100000E+01, 0.174298E+09},
      {0.292481E-07, 0.758028E+00, 0.100000E+01, 0.374772E+08},
      {0.510000E-02, 0.970269E+00, 0.100000E+01, 0.123550E+09},
      {0.510000E-02, 0.980946E+00, 0.100000E+01, 0.123550E+09},

      {0.998086E-06, 0.114746E+01, 0.100000E+01, 0.174298E+09},
      {0.0,          0.758028E+00, 0.100000E+01, 0.374772E+08},
      {0.510000E-02, 0.970269E+00, 0.100000E+01, 0.123550E+09},
      {0.0,          0.980946E+00, 0.100000E+01, 0.123550E+09},

      {0.000000E+00, 1.144941E+00, 1.000000E+00, 0.000000E+00},
      {0.292481E-07, 0.758028E+00, 0.100000E+01, 0.374772E+08},
      {0.510000E-02, 0.970269E+00, 0.100000E+01, 0.123550E+09},
      {0.510000E-02, 0.980946E+00, 0.100000E+01, 0.123550E+09},

      {0.998086E-06, 0.114746E+01, 0.100000E+01, 0.174298E+09},
      {0.292481E-07, 0.758028E+00, 0.100000E+01, 0.374772E+08},
      {0.510000E-02, 0.970269E+00, 0.100000E+01, 0.123550E+09},
      {0.510000E-02, 0.980946E+00, 0.100000E+01, 0.123550E+09}
    };

    // Get data from fortran
    for (Int i = 0; i < Spack::n; ++i) {
      get_cloud_dsd2(gcdd[i]);
    }

    // Sync to device
    view_1d<GetCloudDsd2Data> gcdd_device("gcdd", Spack::n);
    const auto gcdd_host = Kokkos::create_mirror_view(gcdd_device);
    std::copy(&gcdd[0], &gcdd[0] + Spack::n, gcdd_host.data());
    Kokkos::deep_copy(gcdd_device, gcdd_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      // Init pack inputs
      Spack qc, rho, lcldm, nc;
      for (Int s = 0; s < Spack::n; ++s) {
        qc[s]    = gcdd_device(s).qc;
        rho[s]   = gcdd_device(s).rho;
        lcldm[s] = gcdd_device(s).lcldm;
        nc[s]    = gcdd_device(s).nc_in;
      }

      Smask gt_small(qc > qsmall);
      Spack mu_c(0.0), nu(0.0), lamc(0.0), cdist(0.0), cdist1(0.0);
      Functions::get_cloud_dsd2(gt_small, qc, nc, mu_c, rho, nu, dnu, lamc, cdist, cdist1, lcldm);

      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        gcdd_device(s).nc_out = nc[s];
        gcdd_device(s).mu_c = mu_c[s];
        gcdd_device(s).nu = nu[s];
        gcdd_device(s).lamc = lamc[s];
        gcdd_device(s).cdist = cdist[s];
        gcdd_device(s).cdist1 = cdist1[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(gcdd_host, gcdd_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(gcdd[s].nc_out == gcdd_host(s).nc_out);
      REQUIRE(gcdd[s].mu_c   == gcdd_host(s).mu_c);
      REQUIRE(gcdd[s].nu     == gcdd_host(s).nu);
      REQUIRE(gcdd[s].lamc   == gcdd_host(s).lamc);
      REQUIRE(gcdd[s].cdist  == gcdd_host(s).cdist);
      REQUIRE(gcdd[s].cdist1 == gcdd_host(s).cdist1);
    }
  }

  static void run_cloud_phys()
  {
    // TODO
  }

  static void run_rain_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    constexpr Scalar qsmall = C::QSMALL;
    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    GetRainDsd2Data grdd[max_pack_size] = {
      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.156316E-03, 0.100000E+01, 0.194363E+03},
      {0.148647E-03, 0.100000E+01, 0.184827E+03},

      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.0,          0.100000E+01, 0.124340E+05},
      {0.156316E-03, 0.100000E+01, 0.194363E+03},
      {0.0,          0.100000E+01, 0.184827E+03},

      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.156316E-03, 0.100000E+01, 0.194363E+03},
      {0.148647E-03, 0.100000E+01, 0.184827E+03},

      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.100000E-01, 0.100000E+01, 0.124340E+05},
      {0.156316E-03, 0.100000E+01, 0.194363E+03},
      {0.148647E-03, 0.100000E+01, 0.184827E+03}
    };

    // Get data from fortran
    for (Int i = 0; i < Spack::n; ++i) {
      get_rain_dsd2(grdd[i]);
    }

    // Sync to device
    KTH::view_1d<GetRainDsd2Data> grdd_host("grdd_host", Spack::n);
    view_1d<GetRainDsd2Data> grdd_device("grdd_host", Spack::n);
    std::copy(&grdd[0], &grdd[0] + Spack::n, grdd_host.data());
    Kokkos::deep_copy(grdd_device, grdd_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      // Init pack inputs
      Spack qr, rcldm, nr;
      for (Int s = 0; s < Spack::n; ++s) {
        qr[s]    = grdd_device(s).qr;
        rcldm[s] = grdd_device(s).rcldm;
        nr[s]    = grdd_device(s).nr_in;
      }

      Smask gt_small(qr > qsmall);
      Spack mu_r(0.0), lamr(0.0), cdistr(0.0), logn0r(0.0);
      Functions::get_rain_dsd2(gt_small, qr, nr, mu_r, lamr, cdistr, logn0r, rcldm);

      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        grdd_device(s).nr_out = nr[s];
        grdd_device(s).mu_r = mu_r[s];
        grdd_device(s).lamr = lamr[s];
        grdd_device(s).cdistr = cdistr[s];
        grdd_device(s).logn0r = logn0r[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(grdd_host, grdd_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(grdd[s].nr_out == grdd_host(s).nr_out);
      REQUIRE(grdd[s].mu_r   == grdd_host(s).mu_r);
      REQUIRE(grdd[s].lamr   == grdd_host(s).lamr);
      REQUIRE(grdd[s].cdistr == grdd_host(s).cdistr);
      REQUIRE(grdd[s].logn0r == grdd_host(s).logn0r);
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
