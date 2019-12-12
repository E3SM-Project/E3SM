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
  for (Int i = 0; i < max_pack_size; ++i) {
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

 
  ///BALLI-new code
  
  template <typename D>
  struct UnitWrap::UnitTest<D>::TestP3UpdatePrognosticIce
  {
  
     
    KOKKOS_FUNCTION static void  update_prognostic_ice_unit_bfb_tests(){
      using KTH = KokkosTypes<HostDevice>;
      
      static constexpr Int max_pack_size = 16;
      
      REQUIRE(Spack::n <= max_pack_size);
      
      //fortran generated data is input to the following 
      // Jim's print statement for fortran code
      //print '("JGFF get_cloud_dsd2: nc=",ES20.12," mu_c=",ES20.12," nu=",ES20.12," lamc=",ES20.12," cdist=",ES20.12," cdist1=",ES20.12)', nc, mu_c, nu, lamc, cdist, cdist1

      P3UpdatePrognosticIceData pupidc[max_pack_size] = {
	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},


	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},


	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},


	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, 1, 0, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

	{ 4.9078102252657532E-019, 5.9698869870615695E-009, 0.0000000000000000, 3796125.0713374666, 1.7736894054953034E-004,
	  0.0000000000000000, 1.4300045023665814E-033, 8.1978258403749801E-026, 3.0430945682517439E-011, 3.4778223637162733E-004,
	  0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 5.1385646807762583E-007, 
	  0.0000000000000000, 0.0000000000000000, 2.7052763117748554E-002, 0.0000000000000000, 1.9208860554273338E-010,
	  1.0686413520835873, 333700.00000000000, 2834700.0000000000, true, false, 1800.0000000000000, 0.20000000000000001,
	  453.11687137025223, 287.20419568524414, 5.0000000000000001E-003, 6.4285714285714288E-005, 123444710.42445230,
	  7.3684210526315801E-006, 1000000.0000000000, 1.0000000000000000E-004, 1000000.0000000000, 6.4285714285714288E-005,
	  1.0000000000000000E-002},

      };
      
      // Sync to device
      KTH::view_1d<P3UpdatePrognosticIceData> pupidc_host("pupidc_host", Spack::n);
      view_1d<P3UpdatePrognosticIceData> pupidc_device("pupidc_host", Spack::n);
      
      // This copy only copies the input variables.
      std::copy(&pupidc[0], &pupidc[0] + Spack::n, pupidc_host.data());
      Kokkos::deep_copy(pupidc_device, pupidc_host);

      // Get data from fortran
      for (Int i = 0; i < max_pack_size; ++i) {
	update_prognostic_ice(pupidc[i]);
      }

      // This copy also copies the output from the fortran function into the host view. These values
      // are need to check the values returned from
      std::copy(&pupidc[0], &pupidc[0] + Spack::n, pupidc_host.data());      
      
      // Run the lookup from a kernel and copy results back to host
      Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
	  // Init pack inputs
	  Spack qcheti, qccol, qcshd, nccol, ncheti, ncshdc, qrcol, nrcol, qrheti, nrheti, nrshdr,
            qimlt, nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, exner, xlf, xxls,
            nmltratio, rhorime_c, th, qv, qc, nc, qr, nr, qitot, nitot, qirim, birim;
	  Scalar dt;
	  bool log_predictNc, log_wetgrowth;

	  // variables with single values assigned outside of the for loop
	  dt            = pupidc_device(0).dt;
	  log_predictNc = pupidc_device(0).log_predictNc;
	  log_wetgrowth = pupidc_device(0).log_wetgrowth;

	  for (Int s = 0; s < Spack::n; ++s) {
	    
	    qcheti[s] = pupidc_device(s).qcheti;
	    qccol[s]  = pupidc_device(s).qccol; 
	    qcshd[s]  = pupidc_device(s).qcshd; 
	    nccol[s]  = pupidc_device(s).nccol; 
	    ncheti[s] = pupidc_device(s).ncheti;
	    ncshdc[s] = pupidc_device(s).ncshdc;
	    qrcol[s]  = pupidc_device(s).qrcol;
	    nrcol[s]  = pupidc_device(s).nrcol;
	    qrheti[s] = pupidc_device(s).qrheti;
	    nrheti[s] = pupidc_device(s).nrheti;
	    nrshdr[s] = pupidc_device(s).nrshdr;	      
	    qimlt[s]  = pupidc_device(s).qimlt;
	    nimlt[s]  = pupidc_device(s).nimlt;
	    qisub[s]  = pupidc_device(s).qisub;
	    qidep[s]  = pupidc_device(s).qidep;
	    qinuc[s]  = pupidc_device(s).qinuc;
	    ninuc[s]  = pupidc_device(s).ninuc;
	    nislf[s]  = pupidc_device(s).nislf;
	    nisub[s]  = pupidc_device(s).nisub;
	    qiberg[s] = pupidc_device(s).qiberg;
	    exner[s]  = pupidc_device(s).exner;
	    xlf[s]    = pupidc_device(s).xlf;
	    xxls[s]   = pupidc_device(s).xxls;

	    //dt[s]    = pupidc_device(s).dt;
	    nmltratio[s] = pupidc_device(s).nmltratio;
	    rhorime_c[s] = pupidc_device(s).rhorime_c;
	    th[s]    = pupidc_device(s).th;
	    qv[s]    = pupidc_device(s).qv;
	    qc[s]    = pupidc_device(s).qc;
	    nc[s]    = pupidc_device(s).nc;
	    qr[s]    = pupidc_device(s).qr;
	    nr[s]    = pupidc_device(s).nr;	      
	    qitot[s] = pupidc_device(s).qitot;
	    nitot[s] = pupidc_device(s).nitot;
	    qirim[s] = pupidc_device(s).qirim;
	    birim[s] = pupidc_device(s).birim;
	  }

	  Functions::update_prognostic_ice(qcheti, qccol, qcshd, nccol, ncheti,ncshdc,
					   qrcol,   nrcol,  qrheti,  nrheti,  nrshdr,
					   qimlt,  nimlt,  qisub,  qidep,  qinuc,  ninuc,
					   nislf,  nisub,  qiberg,  exner,  xxls,  xlf,
					   log_predictNc, log_wetgrowth,  dt,  nmltratio,
					   rhorime_c, th, qv, qitot, nitot, qirim,
					   birim, qc, nc, qr, nr);
	  
	  // Copy results back into views
	  pupidc_device(0).dt            = dt;
	  pupidc_device(0).log_predictNc = log_predictNc;
	  pupidc_device(0).log_wetgrowth = log_wetgrowth;
	  for (Int s = 0; s < Spack::n; ++s) {
	    
	    pupidc_device(s).qcheti = qcheti[s]; 	
	    pupidc_device(s).qccol  = qccol[s]; 		   
	    pupidc_device(s).qcshd  = qcshd[s];
	    pupidc_device(s).nccol  = nccol[s];		   
	    pupidc_device(s).ncheti = ncheti[s];		   
	    pupidc_device(s).ncshdc = ncshdc[s];		   
	    pupidc_device(s).qrcol  = qrcol[s];			   
	    pupidc_device(s).nrcol  = nrcol[s];			   
	    pupidc_device(s).qrheti = qrheti[s];		   
	    pupidc_device(s).nrheti = nrheti[s];		   
	    pupidc_device(s).nrshdr = nrshdr[s];		   
	    pupidc_device(s).qimlt  = qimlt[s];			   
	    pupidc_device(s).nimlt  = nimlt[s];			   
	    pupidc_device(s).qisub  = qisub[s];			   
	    pupidc_device(s).qidep  = qidep[s];			   
	    pupidc_device(s).qinuc  = qinuc[s];			   
	    pupidc_device(s).ninuc  = ninuc[s];			   
	    pupidc_device(s).nislf  = nislf[s];			   
	    pupidc_device(s).nisub  = nisub[s];			   
	    pupidc_device(s).qiberg = qiberg[s];		   
	    pupidc_device(s).exner  = exner[s];			   
	    pupidc_device(s).xlf    = xlf[s];			   
	    pupidc_device(s).xxls   = xxls[s];  		   
	    
	    pupidc_device(s).nmltratio = nmltratio[s];
	    pupidc_device(s).rhorime_c = rhorime_c[s];
	    pupidc_device(s).th    = th[s];
	    pupidc_device(s).qv	   = qv[s];		 
	    pupidc_device(s).qc	   = qc[s];		 
	    pupidc_device(s).nc	   = nc[s];		 
	    pupidc_device(s).qr	   = qr[s];		 
	    pupidc_device(s).nr    = nr[s];	  	 
	    pupidc_device(s).qitot = qitot[s];
	    pupidc_device(s).nitot = nitot[s];
	    pupidc_device(s).qirim = qirim[s];
	    pupidc_device(s).birim = birim[s];
	  }
	  
	});
      
      // Sync back to host
      Kokkos::deep_copy(pupidc_host, pupidc_device);

      // Validate results
      //First verify the single value variables and then the ones in a pack
      REQUIRE(pupidc[0].dt            == pupidc_host(0).dt);
      REQUIRE(pupidc[0].log_predictNc == pupidc_host(0).log_predictNc);
      REQUIRE(pupidc[0].log_wetgrowth == pupidc_host(0).log_wetgrowth);

      for (Int s = 0; s < Spack::n; ++s) {

	REQUIRE(pupidc[s].qcheti== pupidc_host(s).qcheti);
	REQUIRE(pupidc[s].qccol == pupidc_host(s).qccol);
	REQUIRE(pupidc[s].qcshd == pupidc_host(s).qcshd);
	REQUIRE(pupidc[s].nccol == pupidc_host(s).nccol);
	REQUIRE(pupidc[s].ncheti== pupidc_host(s).ncheti);
	REQUIRE(pupidc[s].ncshdc== pupidc_host(s).ncshdc);
	REQUIRE(pupidc[s].qrcol == pupidc_host(s).qrcol);
	REQUIRE(pupidc[s].nrcol == pupidc_host(s).nrcol);
	REQUIRE(pupidc[s].qrheti== pupidc_host(s).qrheti);
	REQUIRE(pupidc[s].nrheti== pupidc_host(s).nrheti);
	REQUIRE(pupidc[s].nrshdr== pupidc_host(s).nrshdr);
	REQUIRE(pupidc[s].qimlt == pupidc_host(s).qimlt);
	REQUIRE(pupidc[s].nimlt == pupidc_host(s).nimlt);
	REQUIRE(pupidc[s].qisub == pupidc_host(s).qisub);
	REQUIRE(pupidc[s].qidep == pupidc_host(s).qidep);
	REQUIRE(pupidc[s].qinuc == pupidc_host(s).qinuc);
	REQUIRE(pupidc[s].ninuc == pupidc_host(s).ninuc);
	REQUIRE(pupidc[s].nislf == pupidc_host(s).nislf);
	REQUIRE(pupidc[s].nisub == pupidc_host(s).nisub);
	REQUIRE(pupidc[s].qiberg== pupidc_host(s).qiberg);
	REQUIRE(pupidc[s].exner == pupidc_host(s).exner);
	REQUIRE(pupidc[s].xlf   == pupidc_host(s).xlf);
	REQUIRE(pupidc[s].xxls  == pupidc_host(s).xxls);

	REQUIRE(pupidc[s].nmltratio == pupidc_host(s).nmltratio);
	REQUIRE(pupidc[s].rhorime_c == pupidc_host(s).rhorime_c);
	REQUIRE(pupidc[s].qc        == pupidc_host(s).qc); 
	REQUIRE(pupidc[s].nr        == pupidc_host(s).nr);
	REQUIRE(pupidc[s].qr        == pupidc_host(s).qr);     
	REQUIRE(pupidc[s].qv        == pupidc_host(s).qv);  
	REQUIRE(pupidc[s].nc        == pupidc_host(s).nc);
	REQUIRE(pupidc[s].qitot     == pupidc_host(s).qitot); 
	REQUIRE(pupidc[s].nitot     == pupidc_host(s).nitot);
	REQUIRE(pupidc[s].qirim     == pupidc_host(s).qirim);
	REQUIRE(pupidc[s].birim     == pupidc_host(s).birim );
	REQUIRE(pupidc[s].th        == pupidc_host(s).th);   

	}
    }
  
    static void run_bfb(){
      update_prognostic_ice_unit_bfb_tests();
    }

  };//TestP3UpdatePrognosticIce
}//namespace unit_test 
}//namespace p3 
}//namespace scream 

namespace {

TEST_CASE("p3_functions", "[p3_functions]")
{
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Func::run();
}

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_physics();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_bfb();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3UpdatePrognosticIce::run_bfb();
}

} // namespace
