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

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3CloudWaterAutoconversion : public UnitWrap::UnitTest<D>::Base
{

void cloud_water_autoconversion_unit_bfb_tests() {

  CloudWaterAutoconversionData cwadc[max_pack_size] = {
    // rho, qc_incld, nc_incld, inv_qc_relvar
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

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      cwadc[i].read(Base::m_fid);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Spack::n;

    // Init pack inputs
    Spack rho, inv_rho, qc_incld, nc_incld, qr_incld, mu_c, nu, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, inv_qc_relvar;
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      rho[s]                 = cwadc_device(vs).rho;
      qc_incld[s]            = cwadc_device(vs).qc_incld;
      nc_incld[s]            = cwadc_device(vs).nc_incld;
      inv_qc_relvar[s]       = cwadc_device(vs).inv_qc_relvar;
      qc2qr_autoconv_tend[s] = cwadc_device(vs).qc2qr_autoconv_tend;
      nc2nr_autoconv_tend[s] = cwadc_device(vs).nc2nr_autoconv_tend;
      ncautr[s]              = cwadc_device(vs).ncautr;
    }

    Functions::cloud_water_autoconversion(
        rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend,
        nc2nr_autoconv_tend, ncautr,
        p3::Functions<Real,DefaultDevice>::P3Runtime());

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Spack::n; ++s, ++vs) {
      cwadc_device(vs).rho                 = rho[s];
      cwadc_device(vs).qc_incld            = qc_incld[s];
      cwadc_device(vs).nc_incld            = nc_incld[s];
      cwadc_device(vs).inv_qc_relvar       = inv_qc_relvar[s];
      cwadc_device(vs).qc2qr_autoconv_tend = qc2qr_autoconv_tend[s];
      cwadc_device(vs).nc2nr_autoconv_tend = nc2nr_autoconv_tend[s];
      cwadc_device(vs).ncautr              = ncautr[s];
    }
  });

  // Sync back to host
  Kokkos::deep_copy(cwadc_host, cwadc_device);

  // Validate results
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(cwadc[s].rho                  == cwadc_host(s).rho);
      REQUIRE(cwadc[s].qc_incld             == cwadc_host(s).qc_incld);
      REQUIRE(cwadc[s].nc_incld             == cwadc_host(s).nc_incld);
      REQUIRE(cwadc[s].inv_qc_relvar        == cwadc_host(s).inv_qc_relvar);
      REQUIRE(cwadc[s].qc2qr_autoconv_tend  == cwadc_host(s).qc2qr_autoconv_tend);
      REQUIRE(cwadc[s].nc2nr_autoconv_tend  == cwadc_host(s).nc2nr_autoconv_tend);
      REQUIRE(cwadc[s].ncautr               == cwadc_host(s).ncautr);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      cwadc_host(s).write(Base::m_fid);
    }
  }
}

  void run_bfb() {
    cloud_water_autoconversion_unit_bfb_tests();
  }

  KOKKOS_FUNCTION static void autoconversion_is_positive(const Int &i, Int &errors){

    const Spack rho(1.0), inv_qc_relvar(1.0);
    Spack qc_incld, nc_incld(1e7), qc2qr_autoconv_tend(0.0), nc2nr_autoconv_tend(0.0), ncautr(0.0);
    for(int si=0; si<Spack::n; ++si){
        qc_incld[si] = 1e-6 * i * Spack::n + si;
    }
    Functions::cloud_water_autoconversion(
        rho, qc_incld, nc_incld, inv_qc_relvar, qc2qr_autoconv_tend,
        nc2nr_autoconv_tend, ncautr,
        p3::Functions<Real, DefaultDevice>::P3Runtime());
    if((qc2qr_autoconv_tend < 0.0).any()) {
      errors++;
    }
  }

  void run_physics(){

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

namespace {

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion;

  T t;
  t.run_physics();
  t.run_bfb();
}

} // namespace
