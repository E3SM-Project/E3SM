#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

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
struct UnitWrap::UnitTest<D>::TestIceCollection {

  static void run_ice_cldliq_bfb()
  {
    // Read in tables
    view_2d_table vn_table;
    view_2d_table vm_table;
    view_2d_table revap_table;
    view_1d_table mu_r_table; view_dnu_table dnu;
    Functions::init_kokkos_tables(vn_table, vm_table, revap_table, mu_r_table, dnu);

    // Load some lookup inputs, need at least one per pack value
    IceCldliqCollectionData cldliq[max_pack_size] = {
      //  rho      temp      rhofaci     f1pr04     qitot      qc           nitot      nc
      {4.056E-03, 4.02E+01, 8.852E-01, 0.174E+00, 1.221E-14, 5.100E-03, 9.558E+04, 9.952E+05},
      {6.852E-02, 5.01E+01, 8.852E-01, 0.374E+00, 1.221E-15, 4.100E-15, 9.558E+04, 9.952E+05},
      {8.852E-02, 6.00E+01, 8.900E-01, 0.123E+00, 1.221E-12, 3.100E-03, 9.558E+04, 9.952E+05},
      {1.902E-01, 1.00E+02, 8.900E-01, 0.123E+00, 1.221E-15, 2.100E-03, 9.558E+04, 9.952E+05},

      {2.201E-01, 2.00E+02, 0.100E+01, 0.174E+00, 1.221E-10, 1.100E-15, 2.558E+05, 9.952E+06},
      {3.502E-01, 3.00E+02, 0.100E+01, 0.374E+00, 1.221E-15, 8.100E-15, 2.558E+05, 9.952E+06},
      {4.852E-01, 5.00E+02, 0.100E+01, 0.123E+00, 1.221E-08, 4.100E-04, 2.558E+05, 9.952E+06},
      {5.852E-01, 8.00E+02, 0.100E+01, 0.123E+00, 1.221E-15, 2.100E-04, 2.558E+05, 9.952E+06},

      {6.852E-01, 1.00E+03, 0.950E+00, 0.150E+00, 1.221E-06, 9.952E-05, 4.596E+05, 1.734E+07},
      {7.852E-01, 2.00E+03, 0.950E+00, 0.374E+00, 1.221E-15, 4.952E-05, 4.596E+05, 1.734E+07},
      {8.852E-01, 4.00E+03, 0.950E+00, 0.123E+00, 1.221E-04, 1.952E-15, 4.596E+05, 1.734E+07},
      {9.852E-01, 6.00E+03, 0.950E+00, 0.123E+00, 1.221E-03, 9.952E-06, 4.596E+05, 1.734E+07},

      {1.002E+01, 1.00E+04, 1.069E+00, 0.174E+00, 1.221E-15, 6.952E-06, 6.596E+05, 1.734E+08},
      {1.152E+01, 2.00E+04, 1.069E+00, 0.374E+00, 1.221E-02, 3.952E-06, 6.596E+05, 1.734E+08},
      {1.252E+01, 4.00E+04, 1.069E+00, 0.123E+00, 1.221E-02, 1.952E-06, 6.596E+05, 1.734E+08},
      {1.352E+01, 8.00E+04, 1.069E+00, 0.123E+00, 1.221E-02, 9.952E-15, 6.596E+05, 1.734E+08}
    };

    // Sync to device
    view_1d<IceCldliqCollectionData> cldliq_device("cldliq", max_pack_size);
    const auto cldliq_host = Kokkos::create_mirror_view(cldliq_device);
    std::copy(&cldliq[0], &cldliq[0] + max_pack_size, cldliq_host.data());
    Kokkos::deep_copy(cldliq_device, cldliq_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      ice_cldliq_collection(cldliq[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho, temp, rhofaci, f1pr04, qitot_incld, qc_incld, nitot_incld, nc_incld;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]            = cldliq_device(vs).rho;
        temp[s]           = cldliq_device(vs).temp;
        rhofaci[s]        = cldliq_device(vs).rhofaci;
        f1pr04[s]         = cldliq_device(vs).f1pr04;
        qitot_incld[s]    = cldliq_device(vs).qitot_incld;
        qc_incld[s]       = cldliq_device(vs).qc_incld;
        nitot_incld[s]    = cldliq_device(vs).nitot_incld;
        nc_incld[s]       = cldliq_device(vs).nc_incld;
      }

      Spack qccol{0.0};
      Spack nccol{0.0};
      Spack qcshd{0.0};
      Spack ncshdc{0.0};

      Functions::ice_cldliq_collection(rho, temp, rhofaci, f1pr04, qitot_incld,
                                       qc_incld, nitot_incld, nc_incld,
                                       qccol, nccol, qcshd, ncshdc);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        cldliq_device(vs).qccol  = qccol[s];
        cldliq_device(vs).nccol  = nccol[s];
        cldliq_device(vs).qcshd  = qcshd[s];
        cldliq_device(vs).ncshdc = ncshdc[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(cldliq_host, cldliq_device);

    // Validate results
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cldliq[s].qccol   == cldliq_host(s).qccol);
      REQUIRE(cldliq[s].nccol   == cldliq_host(s).nccol);
      REQUIRE(cldliq[s].qcshd   == cldliq_host(s).qcshd);
      REQUIRE(cldliq[s].ncshdc  == cldliq_host(s).ncshdc);
    }
  }

  static void run_ice_cldliq_phys()
  {
    // TODO
  }

  static void run_ice_rain_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    IceRainCollectionData rain[max_pack_size] = {
      //  rho      temp      rhofaci     logn0r     f1pr07    f1pr08        qitot      nitot    qr (required)
      {4.056E-03, 4.02E+01, 8.852E-01, 0.174E+00, 1.221E-14, 5.100E-03, 9.558E-04, 9.952E+02, 5.100E-03},
      {6.852E-02, 5.01E+01, 8.852E-01, 0.374E+00, 1.221E-13, 4.100E-03, 9.558E-15, 9.952E+02, 5.100E-15},
      {8.852E-02, 6.00E+01, 8.900E-01, 0.123E+00, 1.221E-12, 3.100E-03, 9.558E-04, 9.952E+02, 5.100E-03},
      {1.902E-01, 1.00E+02, 8.900E-01, 0.123E+00, 1.221E-11, 2.100E-03, 9.558E-04, 9.952E+02, 5.100E-15},

      {2.201E-01, 2.00E+02, 0.100E+01, 0.174E+00, 1.221E-10, 1.100E-03, 2.558E-05, 9.952E+02, 5.100E-15},
      {3.502E-01, 3.00E+02, 0.100E+01, 0.374E+00, 1.221E-09, 8.100E-04, 2.558E-15, 9.952E+02, 5.100E-15},
      {4.852E-01, 5.00E+02, 0.100E+01, 0.123E+00, 1.221E-08, 4.100E-04, 2.558E-05, 9.952E+02, 5.100E-03},
      {5.852E-01, 8.00E+02, 0.100E+01, 0.123E+00, 1.221E-07, 2.100E-04, 2.558E-05, 9.952E+02, 5.100E-03},

      {6.852E-01, 1.00E+03, 0.950E+00, 0.150E+00, 1.221E-06, 9.952E-05, 4.596E-05, 1.734E+03, 5.100E-15},
      {7.852E-01, 2.00E+03, 0.950E+00, 0.374E+00, 1.221E-05, 4.952E-05, 4.596E-15, 1.734E+03, 5.100E-15},
      {8.852E-01, 4.00E+03, 0.950E+00, 0.123E+00, 1.221E-04, 1.952E-05, 4.596E-05, 1.734E+03, 5.100E-03},
      {9.852E-01, 6.00E+03, 0.950E+00, 0.123E+00, 1.221E-03, 9.952E-06, 4.596E-05, 1.734E+03, 5.100E-03},

      {1.002E+01, 1.00E+04, 1.069E+00, 0.174E+00, 1.221E-02, 6.952E-06, 6.596E-05, 1.734E+03, 5.100E-15},
      {1.152E+01, 2.00E+04, 1.069E+00, 0.374E+00, 1.221E-02, 3.952E-06, 6.596E-15, 1.734E+03, 5.100E-03},
      {1.252E+01, 4.00E+04, 1.069E+00, 0.123E+00, 1.221E-02, 1.952E-06, 6.596E-05, 1.734E+03, 5.100E-15},
      {1.352E+01, 8.00E+04, 1.069E+00, 0.123E+00, 1.221E-02, 9.952E-07, 6.596E-05, 1.734E+03, 5.100E-03}
    };

    // Sync to device
    KTH::view_1d<IceRainCollectionData> rain_host("rain_host", max_pack_size);
    view_1d<IceRainCollectionData> rain_device("rain_host", max_pack_size);
    std::copy(&rain[0], &rain[0] + max_pack_size, rain_host.data());
    Kokkos::deep_copy(rain_device, rain_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      ice_rain_collection(rain[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho, temp, rhofaci, logn0r, f1pr07, f1pr08, qitot_incld, nitot_incld, qr_incld;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]         = rain_device(vs).rho;
        temp[s]        = rain_device(vs).temp;
        rhofaci[s]     = rain_device(vs).rhofaci;
        logn0r[s]      = rain_device(vs).logn0r;
        f1pr07[s]      = rain_device(vs).f1pr07;
        f1pr08[s]      = rain_device(vs).f1pr08;
        qitot_incld[s] = rain_device(vs).qitot_incld;
        nitot_incld[s] = rain_device(vs).nitot_incld;
        qr_incld[s]    = rain_device(vs).qr_incld;
      }

      Spack qrcol(0.0), nrcol(0.0);
      Functions::ice_rain_collection(rho, temp, rhofaci, logn0r, f1pr07, f1pr08,
                                     qitot_incld, nitot_incld, qr_incld,
                                     qrcol, nrcol);

      // Copy results back into views
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rain_device(vs).qrcol = qrcol[s];
        rain_device(vs).nrcol = nrcol[s];
      }
    });

    // Sync back to host
    Kokkos::deep_copy(rain_host, rain_device);

    // Validate results
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(rain[s].qrcol == rain_host(s).qrcol);
      REQUIRE(rain[s].nrcol == rain_host(s).nrcol);
    }
  }

  static void run_ice_rain_phys()
  {
    // TODO
  }

  static void run_ice_self_bfb()
  {
    using KTH = KokkosTypes<HostDevice>;

    IceSelfCollectionData self[max_pack_size] = {
     //   rho      rhofaci    f1pr03     eii       qirim      qitot        nitot
      {4.056E-03, 8.852E-01, 0.174E+00, 1.221E-14, 5.100E-03, 9.558E-04, 9.952E+03},
      {6.852E-02, 8.852E-01, 0.374E+00, 1.221E-13, 0.000E+00, 9.558E-15, 9.952E+03},
      {8.852E-02, 8.900E-01, 0.123E+00, 1.221E-12, 3.100E-03, 9.558E-04, 9.952E+03},
      {1.902E-01, 9.900E-01, 0.123E+00, 1.221E-11, 2.100E-03, 9.558E-15, 9.952E+03},

      {2.201E-01, 0.100E+01, 0.174E+00, 1.221E-10, 1.100E-03, 2.558E-15, 9.952E+04},
      {3.502E-01, 0.100E+01, 0.374E+00, 1.221E-09, 0.000E+00, 2.558E-15, 9.952E+04},
      {4.852E-01, 0.100E+01, 0.123E+00, 1.221E-08, 4.100E-04, 2.558E-05, 9.952E+04},
      {5.852E-01, 0.100E+01, 0.123E+00, 1.221E-07, 2.100E-04, 2.558E-05, 9.952E+04},

      {6.852E-01, 0.950E+00, 0.150E+00, 1.221E-06, 0.000E+00, 4.596E-05, 1.734E+04},
      {7.852E-01, 0.950E+00, 0.374E+00, 1.221E-05, 4.952E-05, 4.596E-05, 1.734E+04},
      {8.852E-01, 0.950E+00, 0.123E+00, 1.221E-04, 1.952E-05, 4.596E-05, 1.734E+04},
      {9.852E-01, 0.950E+00, 0.123E+00, 1.221E-03, 9.952E-06, 4.596E-15, 1.734E+04},

      {1.002E+01, 1.069E+00, 0.174E+00, 1.221E-02, 6.952E-06, 6.596E-15, 1.734E+04},
      {1.152E+01, 1.069E+00, 0.374E+00, 1.221E-02, 3.952E-06, 6.596E-05, 1.734E+04},
      {1.252E+01, 1.069E+00, 0.123E+00, 1.221E-02, 1.952E-06, 6.596E-05, 1.734E+04},
      {1.352E+01, 1.069E+00, 0.123E+00, 1.221E-02, 9.952E-07, 6.596E-05, 1.734E+04}
    };

    // Sync to device
    KTH::view_1d<IceSelfCollectionData> self_host("self_host", max_pack_size);
    view_1d<IceSelfCollectionData> self_device("self_host", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      ice_self_collection(self[i]);
    }

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Spack::n;

      // Init pack inputs
      Spack rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld;
      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        rho[s]         = self_device(vs).rho;
        rhofaci[s]     = self_device(vs).rhofaci;
        f1pr03[s]      = self_device(vs).f1pr03;
        eii[s]         = self_device(vs).eii;
        qirim_incld[s] = self_device(vs).qirim_incld;
        qitot_incld[s] = self_device(vs).qitot_incld;
        nitot_incld[s] = self_device(vs).nitot_incld;
      }

      Spack nislf{0.0};
      Functions::ice_self_collection(rho, rhofaci, f1pr03, eii, qirim_incld, qitot_incld, nitot_incld,
                                     nislf);

      for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
        self_device(vs).nislf = nislf[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(self[s].nislf == self_host(s).nislf);
    }
  }


  static void run_ice_self_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_ice_cldliq", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceCollection;

  TD::run_ice_cldliq_phys();
  TD::run_ice_cldliq_bfb();
}

TEST_CASE("p3_ice_rain", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceCollection;

  TD::run_ice_rain_phys();
  TD::run_ice_rain_bfb();
}

TEST_CASE("p3_ice_self", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceCollection;
  TD::run_ice_self_phys();
  TD::run_ice_self_bfb();
}
}
