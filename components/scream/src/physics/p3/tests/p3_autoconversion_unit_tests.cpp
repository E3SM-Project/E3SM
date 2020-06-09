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

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3CloudWaterAutoconversion
{

static void  cloud_water_autoconversion_unit_bfb_tests(){

  CloudWaterAutoconversionData cwadc[max_pack_size] = {
    // rho, qc_incld, nc_incld, qc_relvar
    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},

    {9.703E-01, 5.100E-03, 2.061E+08, 1.0},
    {1.006E+00,  5.100E-03, 1.988E+08,1.0},
    //{1.139E+00 },
    {1.139E+00,  0.0,       0.0,      1.0},
    {1.151E+00,  1.000E-06, 1.737E+08,1.0},
  };

  // Sync to device
  view_1d<CloudWaterAutoconversionData> cwadc_device("cwadc", max_pack_size);
  auto cwadc_host = Kokkos::create_mirror_view(cwadc_device);

  // This copy only copies the input variables.
  std::copy(&cwadc[0], &cwadc[0] + max_pack_size, cwadc_host.data());
  Kokkos::deep_copy(cwadc_device, cwadc_host);

  // Get data from fortran
  for (Int i = 0; i < max_pack_size; ++i) {
    cloud_water_autoconversion(cwadc[i]);
  }

    // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, nc_incld, qr_incld, mu_c, nu, qcaut, ncautc, ncautr, qc_relvar;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]       = cwadc_device(vs).rho;
      qc_incld[s]  = cwadc_device(vs).qc_incld;
      nc_incld[s]  = cwadc_device(vs).nc_incld;
      qc_relvar[s] = cwadc_device(vs).qc_relvar;
      qcaut[s]     = cwadc_device(vs).qcaut;
      ncautc[s]    = cwadc_device(vs).ncautc;
      ncautr[s]    = cwadc_device(vs).ncautr;
    }

    Functions::cloud_water_autoconversion(rho, qc_incld, nc_incld,
      qc_relvar, qcaut, ncautc, ncautr);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      cwadc_device(vs).rho       = rho[s];
      cwadc_device(vs).qc_incld  = qc_incld[s];
      cwadc_device(vs).nc_incld  = nc_incld[s];
      cwadc_device(vs).qc_relvar = qc_relvar[s];
      cwadc_device(vs).qcaut     = qcaut[s];
      cwadc_device(vs).ncautc    = ncautc[s];
      cwadc_device(vs).ncautr    = ncautr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cwadc_host, cwadc_device);

  // Validate results
  for (Int s = 0; s < max_pack_size; ++s) {
    REQUIRE(cwadc[s].rho       == cwadc_host(s).rho);
    REQUIRE(cwadc[s].qc_incld  == cwadc_host(s).qc_incld);
    REQUIRE(cwadc[s].nc_incld  == cwadc_host(s).nc_incld);
    REQUIRE(cwadc[s].qc_relvar == cwadc_host(s).qc_relvar);
    REQUIRE(cwadc[s].qcaut     == cwadc_host(s).qcaut);
    REQUIRE(cwadc[s].ncautc    == cwadc_host(s).ncautc);
    REQUIRE(cwadc[s].ncautr    == cwadc_host(s).ncautr);
  }
}

  static void run_bfb(){
    cloud_water_autoconversion_unit_bfb_tests();
  }

  KOKKOS_FUNCTION  static void autoconversion_is_positive(const Int &i, Int &errors){

    const Spack rho(1.0), qc_relvar(1.0);
    Spack qc_incld, nc_incld(1e7), qcaut(0.0), ncautc(0.0), ncautr(0.0);
    for(int si=0; si<Spack::n; ++si){
        qc_incld[si] = 1e-6 * i * Spack::n + si;
      }
    Functions::cloud_water_autoconversion(rho, qc_incld, nc_incld, qc_relvar, qcaut, ncautc, ncautr);
        if((qcaut < 0.0).any()){errors++;}
    }

  static void run_physics(){

    int nerr = 0;

    Kokkos::parallel_reduce("TestAutoConversionPositive", 1000, KOKKOS_LAMBDA(const Int& i,  Int& errors) {
      autoconversion_is_positive(i, errors);
    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);

  }

}; //  TestP3CloudWaterAutoconversion

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace{

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_physics();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_bfb();
}

} // namespace

