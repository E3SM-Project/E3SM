#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"
#include "share/util/scream_arch.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit-tests for p3_functions.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestP3Func
{

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure, const Scalar& correct_sat_ice_p,
    const Scalar& correct_sat_liq_p, const Scalar&  correct_mix_ice_r, const Scalar& correct_mix_liq_r, int& errors ){

    const Spack temps(temperature);
    const Spack pres(pressure);

    Spack sat_ice_p = Functions::polysvp1(temps, true);
    Spack sat_liq_p = Functions::polysvp1(temps, false);

    Spack mix_ice_r = Functions::qv_sat(temps, pres, true);
    Spack mix_liq_r = Functions::qv_sat(temps, pres, false);

    // The correct results were computed with double precision, so we need
    // significantly greater tolerance for single precision.
    Scalar tol = (util::is_single_precision<Scalar>::value || util::OnGpu<ExeSpace>::value) ? C::Tol*100 : C::Tol;

    for(int s = 0; s < sat_ice_p.n; ++s){
      // Test vapor pressure
      if (abs(sat_ice_p[s] - correct_sat_ice_p) > tol ) {errors++;}
      if (abs(sat_liq_p[s] - correct_sat_liq_p) > tol)  {errors++;}
      //Test mixing-ratios
      if (abs(mix_ice_r[s] -  correct_mix_ice_r) > tol ) {errors++;}
      if (abs(mix_liq_r[s] -  correct_mix_liq_r) > tol ) {errors++;}
    }
  }

  static void run()
  {
    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors = 0;
      const auto tmelt = C::Tmelt;
      // Test values @ the melting point of H20 @ 1e5 Pa
      saturation_tests(tmelt, 1e5, 610.7960763188032, 610.7960763188032,
         0.003822318507864685,  0.003822318507864685, errors);

      //Test vaules @ 243.15K @ 1e5 Pa
      saturation_tests(243.15, 1e5, 37.98530141245404, 50.98455924912173,
         0.00023634717905493638,  0.0003172707211143376, errors);

      //Test values @ 303.15 @ 1e5 Pa
      saturation_tests(303.15, 1e5, 4242.757341329608, 4242.757341329608,
        0.0275579183092878, 0.0275579183092878, errors);

    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
};


template <typename D>
struct UnitWrap::UnitTest<D>::TestP3CloudWaterAutoconversion
{

static void  cloud_water_autoconversion_unit_bfb_tests(){
  using KTH = KokkosTypes<HostDevice>;

  static constexpr Int max_pack_size = 16;
  REQUIRE(Spack::n <= max_pack_size);

  CloudWaterAutoconversionData cwadc[max_pack_size] = {
    { 0.97026902585098274, 5.1000000000000004e-3, 206128398.07453227},
    { 1.0061301158991891,  5.1000000000000004e-3, 198781446.69316244},
    { 1.1393248270523915 },
    { 1.1512545299884895,  9.9999999999999995e-7, 173723529.23727444},

    { 0.97026902585098274, 5.1000000000000004e-3, 206128398.07453227},
    { 1.0061301158991891,  5.1000000000000004e-3, 198781446.69316244},
    { 1.1393248270523915 },
    { 1.1512545299884895,  9.9999999999999995e-7, 173723529.23727444},

    { 0.97026902585098274, 5.1000000000000004e-3, 206128398.07453227},
    { 1.0061301158991891,  5.1000000000000004e-3, 198781446.69316244},
    { 1.1393248270523915 },
    { 1.1512545299884895,  9.9999999999999995e-7, 173723529.23727444},

    { 0.97026902585098274, 5.1000000000000004e-3, 206128398.07453227},
    { 1.0061301158991891,  5.1000000000000004e-3, 198781446.69316244},
    { 1.1393248270523915 },
    { 1.1512545299884895,  9.9999999999999995e-7, 173723529.23727444},
  };

  // Sync to device
  KTH::view_1d<CloudWaterAutoconversionData> cwadc_host("cwadc_host", Spack::n);
  view_1d<CloudWaterAutoconversionData> cwadc_device("cwadc_host", Spack::n);

  // This copy only copies the input variables.
  std::copy(&cwadc[0], &cwadc[0] + Spack::n, cwadc_host.data());
  Kokkos::deep_copy(cwadc_device, cwadc_host);

  // Get data from fortran
  for (Int i = 0; i < Spack::n; ++i) {
    cloud_water_autoconversion(cwadc[i]);
  }

  // This copy also copies the output from the fortran function into the host view. These values
  // are need to check the values returned from
  std::copy(&cwadc[0], &cwadc[0] + Spack::n, cwadc_host.data());

    // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
    // Init pack inputs
    Spack rho, inv_rho, qc_incld, nc_incld, qr_incld, mu_c, nu, qcaut, ncautc, ncautr;
    for (Int s = 0; s < Spack::n; ++s) {
      rho[s] = cwadc_device(s).rho;
      qc_incld[s] = cwadc_device(s).qc_incld;
      nc_incld[s] = cwadc_device(s).nc_incld;
      qcaut[s] = cwadc_device(s).qcaut;
      ncautc[s] = cwadc_device(s).ncautc;
      ncautr[s] = cwadc_device(s).ncautr;
    }

    Functions::cloud_water_autoconversion(rho, qc_incld, nc_incld,
      qcaut, ncautc, ncautr);
    // Copy results back into views
    for (Int s = 0; s < Spack::n; ++s) {
      cwadc_device(s).rho = rho[s];
      cwadc_device(s).qc_incld = qc_incld[s];
      cwadc_device(s).nc_incld = nc_incld[s];
      cwadc_device(s).qcaut = qcaut[s];
      cwadc_device(s).ncautc = ncautc[s];
      cwadc_device(s).ncautr = ncautr[s];
    }

  });

    // Sync back to host
    Kokkos::deep_copy(cwadc_host, cwadc_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
       REQUIRE(cwadc[s].rho == cwadc_host(s).rho);
       REQUIRE(cwadc[s].qc_incld == cwadc_host(s).qc_incld);
       REQUIRE(cwadc[s].nc_incld == cwadc_host(s).nc_incld);
       REQUIRE(cwadc[s].qcaut == cwadc_host(s).qcaut);
       REQUIRE(cwadc[s].ncautc == cwadc_host(s).ncautc);
       REQUIRE(cwadc[s].ncautr == cwadc_host(s).ncautr);
     }
}

  static void run_bfb(){
    cloud_water_autoconversion_unit_bfb_tests();
  }

  KOKKOS_FUNCTION  static void autoconversion_is_positive(const Int &i, Int &errors){

    const Spack rho(1.0);
    Spack qc_incld, nc_incld(1e7), qcaut(0.0), ncautc(0.0), ncautr(0.0);
    for(int si=0; si<Spack::n; ++si){
        qc_incld[si] = 1e-6 * i * Spack::n + si;
      }
        Functions::cloud_water_autoconversion(rho, qc_incld, nc_incld, qcaut, ncautc, ncautr);
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

}
}
}

namespace {

TEST_CASE("p3_functions", "[p3_functions]")
{
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Func::run();
}

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_physics();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_bfb();
}

} // namespace
