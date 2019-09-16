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

    const auto sat_ice_p = Functions::polysvp1(temps, true);
    const auto sat_liq_p = Functions::polysvp1(temps, false);

    const auto mix_ice_r = Functions::qv_sat(temps, pres, true);
    const auto mix_liq_r = Functions::qv_sat(temps, pres, false);

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
struct UnitWrap::UnitTest<D>::TestP3Conservation
{

  KOKKOS_FUNCTION static void cloud_water_conservation_tests(int& errors){

    Spack qc(1e-5);
    Spack qcnuc(0.0);
    Scalar dt = 1.1;
    Spack qcaut(1e-4);
    Spack qcacc(0.0);
    Spack qccol(0.0);
    Spack qcheti(0.0);
    Spack qcshd(0.0);
    Spack qiberg(0.0);
    Spack qisub(1.0);
    Spack qidep(1.0);

    const auto ratio = qc/(qcaut * dt);
    Functions::cloud_water_conservation(qc, qcnuc, dt,
    qcaut,qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);


    //Here we check the case where sources > sinks
    if(((qcaut - 1e-4*ratio) != 0.0).any()){errors++;}
    if((qcacc != 0).any()){errors++;}
    if((qccol != 0).any()){errors++;}
    if((qcheti != 0).any()){errors++;}
    if((qcshd != 0).any()){errors++;}
    if((qiberg != 0).any()){errors++;}
    if(((qisub - (1.0 - ratio))!=0.0).any()){errors++;}
    if(((qidep - (1.0 - ratio))!=0.0).any()){errors++;}

    // Now actually check conservation. We are basically checking here that
    // qcaut, the only non-zero sink, is corrected so that within a dt
    // it does not overshoot qc.
    if(((qcaut*dt - qc) != 0.0).any()){errors++;}

    // Check the case where sources > sinks with sinks = 0
    qcaut = 0.0;
    Functions::cloud_water_conservation(qc, qcnuc, dt,
    qcaut,qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);

    //In this case qidep and qisub should be set to zero
    if((qisub != 0.0).any()){errors++;}
    if((qidep != 0.0).any()){errors++;}

  }

  KOKKOS_FUNCTION static void rain_water_conservation_tests(int& errors){

    Spack qr(1e-5);
    Spack qcaut(0.0);
    Spack qcacc(0.0);
    Spack qimlt(0.0);
    Spack qcshd(0.0);
    Scalar dt = 1.1;
    Spack qrevp(1e-4);
    Spack qrcol(0.0);
    Spack qrheti(0.0);

    const auto ratio = qr/(qrevp * dt);
    //Call function being tested
    Functions::rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, dt, qrevp, qrcol, qrheti);

    //Here we check cases where source > sinks and sinks > 1e-20
    if((qcaut!=0.0).any()){errors++;}
    if((qcacc!=0.0).any()){errors++;}
    if((qimlt!=0.0).any()){errors++;}
    if((qcshd!=0.0).any()){errors++;}
    if(((qrevp - 1e-4*ratio) != 0.0).any()){errors++;}
    if((qrcol!=0.0).any()){errors++;}
    if((qrheti!=0.0).any()){errors++;}

    //Now test that conservation has actually been enforced
    if((abs(qrevp * dt - qr)> 0.0).any()){errors++;}

  }


  KOKKOS_FUNCTION static void ice_water_conservation_tests(int& errors){

    Spack qitot(1e-5);
    Spack qidep(0.0);
    Spack qinuc(0.0);
    Spack qrcol(0.0);
    Spack qccol(0.0);
    Spack qrheti(0.0);
    Spack qcheti(0.0);
    Spack qiberg(0.0);
    Scalar dt = 1.1;
    Spack qisub(1e-4);
    Spack qimlt(0.0);

    const auto ratio = qitot/(qisub * dt);

    //Call function being tested
    Functions::ice_water_conservation(qitot, qidep, qinuc, qrcol, qccol, qrheti, qcheti, qiberg, dt, qisub, qimlt);
    //Here we check cases where source > sinks and sinks > 1e-20
    if((qidep!=0.0).any()){errors++;}
    if((qinuc!=0.0).any()){errors++;}
    if((qrcol!=0.0).any()){errors++;}
    if((qccol!=0.0).any()){errors++;}
    if((qrheti!=0.0).any()){errors++;}
    if((qcheti!=0.0).any()){errors++;}
    if((qiberg!=0.0).any()){errors++;}
    if(((qisub - 1e-4 * ratio)!=0).any()){errors++;}
    if((qimlt!=0).any()){errors++;}

    //Now test that conservation has actually been enforced
    if((abs(qisub * dt - qitot)>0.0).any()){errors++;}
  }

  static void run()
  {
    int nerr = 0;

    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("ConservationTests", policy, KOKKOS_LAMBDA(const MemberType& tem, int& errors){

      errors = 0;

      cloud_water_conservation_tests(errors);

      rain_water_conservation_tests(errors);

      ice_water_conservation_tests(errors);

    }, nerr);

    Kokkos::fence();

    REQUIRE(nerr==0);
  }

  KOKKOS_FUNCTION static void cloud_water_conservation_unit_bfb_tests(){

  }

  KOKKOS_FUNCTION static void ice_water_conservation_unit_bfb_tests()
  {
    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    IceWaterConservationData iwdc[max_pack_size] = {
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},

      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 1.9205467584100191e-4},
      {5.0e-8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 1.8234653652173277e-7, 0.0},
      {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 2.3237448636383435e-3},
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0},
    };

    // Sync to device
    KTH::view_1d<IceWaterConservationData> iwdc_host("iwdc_host", Spack::n);
    view_1d<IceWaterConservationData> iwdc_device("iwdc_host", Spack::n);

    // This copy only copies the input variables.
    std::copy(&iwdc[0], &iwdc[0] + Spack::n, iwdc_host.data());
    Kokkos::deep_copy(iwdc_device, iwdc_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      ice_water_conservation(iwdc[i]);
    }

    // This copy also copies the output from the fortran function into the host view. These values
    // are need to check the values returned from
    std::copy(&iwdc[0], &iwdc[0] + Spack::n, iwdc_host.data());

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      // Init pack inputs
      Spack qitot,qidep,qinuc,qiberg,qrcol,qccol,qrheti,qcheti,qisub,qimlt;
      for (Int s = 0; s < Spack::n; ++s) {
        qitot[s]  = iwdc_device(s).qitot;
        qidep[s]  = iwdc_device(s).qidep;
        qinuc[s]  = iwdc_device(s).qinuc;
        qiberg[s] = iwdc_device(s).qiberg;
        qrcol[s]  = iwdc_device(s).qrcol;
        qccol[s]  = iwdc_device(s).qccol;
        qrheti[s] = iwdc_device(s).qrheti;
        qcheti[s] = iwdc_device(s).qcheti;
        qisub[s] = iwdc_device(s).qisub;
        qimlt[s] = iwdc_device(s).qimlt;
      }

      Functions::ice_water_conservation(qitot,qidep,qinuc,qiberg,qrcol,qccol,qrheti,qcheti,iwdc_device(0).dt,qisub,qimlt);
      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        iwdc_device(s).qitot = qitot[s];
        iwdc_device(s).qidep = qidep[s];
        iwdc_device(s).qinuc = qinuc[s];
        iwdc_device(s).qiberg = qiberg[s];
        iwdc_device(s).qrcol = qrcol[s];
        iwdc_device(s).qccol = qccol[s];
        iwdc_device(s).qrheti = qrheti[s];
        iwdc_device(s).qcheti = qcheti[s];
        iwdc_device(s).qisub = qisub[s];
        iwdc_device(s).qimlt = qimlt[s];
      }

    });

    // Sync back to host
    Kokkos::deep_copy(iwdc_host, iwdc_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(iwdc[s].qitot == iwdc_host(s).qitot);
      REQUIRE(iwdc[s].qidep == iwdc_host(s).qidep );
      REQUIRE(iwdc[s].qinuc == iwdc_host(s).qinuc);
      REQUIRE(iwdc[s].qiberg == iwdc_host(s).qiberg);
      REQUIRE(iwdc[s].qrcol  == iwdc_host(s).qrcol);
      REQUIRE(iwdc[s].qccol == iwdc_host(s).qccol);
      REQUIRE(iwdc[s].qrheti == iwdc_host(s).qrheti);
      REQUIRE(iwdc[s].qcheti == iwdc_host(s).qcheti);
      REQUIRE(iwdc[s].qisub == iwdc_host(s).qisub);
      REQUIRE(iwdc[s].qimlt == iwdc_host(s).qimlt);
    }

  }

  KOKKOS_FUNCTION static void rain_water_conservation_unit_bfb_tests(){

  }



  static void run_bfb(){
      //cloud_water_conservation_tests(errors);

      //rain_water_conservation_tests(errors);

      ice_water_conservation_unit_bfb_tests();

  }

};

}
}
}

namespace {

TEST_CASE("p3_functions", "[p3_functions]")
{
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Func::run();
}

TEST_CASE("p3_conservation_test", "[p3_conservation_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Conservation::run();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Conservation::run_bfb(); 
}

} // namespace
