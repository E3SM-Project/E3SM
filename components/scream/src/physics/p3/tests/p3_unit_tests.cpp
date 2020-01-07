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

  KOKKOS_FUNCTION static Scalar condition_number_polysvp1(const Spack& temperature, const bool ice){

    const auto sat0 = Functions::polysvp1(temperature, ice);
    const auto sat1 = Functions::polysvp1(temperature + 1e-4, ice);

    const auto cond = (abs(sat1 - sat0)/abs(sat0)) /(abs(1e-4)/abs(temperature));


    return 1.25 * cond[0];
  }

  KOKKOS_FUNCTION static Spack condition_number_qv_sat(const Spack& temperature, const Spack& pressure, const bool ice){

    const Scalar deltap = 1e0;
    const Scalar deltat = 1e-2;

    const auto qv_sat0 = Functions::qv_sat(temperature, pressure, ice);
    const auto qv_sat_t1 = Functions::qv_sat(temperature + deltat, pressure, ice);
    const auto qv_sat_p1 = Functions::qv_sat(temperature, pressure + deltap , ice);

    Spack  jac[2];
    jac[0] = (qv_sat_t1 - qv_sat0)/deltat;
    jac[1] = (qv_sat_p1 - qv_sat0)/deltap;

    const auto norm_jac = sqrt(jac[0]*jac[0] + jac[1]*jac[1]);

    auto cond = jac[0][0] * temperature/qv_sat0;
    const auto  set_mask = jac[0][0] * temperature/qv_sat0  < jac[1][0] * pressure/qv_sat0; 
    cond.set(set_mask, jac[1][0] * pressure[0]/qv_sat0[0]);

    return cond*1.25;

  }

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure, const Scalar& correct_sat_ice_p,
    const Scalar& correct_sat_liq_p, const Scalar&  correct_mix_ice_r, const Scalar& correct_mix_liq_r, int& errors ){

    const Spack temps(temperature);
    const Spack pres(pressure);

    const auto sat_ice_p = Functions::polysvp1(temps, true);
    const auto cond_ice_p = condition_number_polysvp1(temps, true);
    const auto sat_liq_p = Functions::polysvp1(temps, false);
    const auto cond_liq_p = condition_number_polysvp1(temps, false);

    const auto mix_ice_r = Functions::qv_sat(temps, pres, true);
    const auto cond_ice_r  = condition_number_qv_sat(temps, pres, true);
    const auto mix_liq_r = Functions::qv_sat(temps, pres, false);
    const auto cond_liq_r = condition_number_qv_sat(temps, pres, false);

    for(int s = 0; s < sat_ice_p.n; ++s){
      // Test vapor pressure note that we multipy by numerically computed realtive condition number 
      if (abs(sat_ice_p[s] - correct_sat_ice_p) >  cond_ice_p*C::Tol * abs(correct_sat_ice_p)){errors++;}
      if (abs(sat_liq_p[s] - correct_sat_liq_p) >  cond_liq_p*C::Tol * abs(correct_sat_liq_p)){errors++;}
      //Test mixing-ratios
      if (abs(mix_ice_r[s] -  correct_mix_ice_r) >   cond_ice_r[s] * C::Tol * abs(correct_mix_ice_r)) {errors++; std::cout << mix_ice_r[s] -  correct_mix_ice_r <<  "\t" << C::Tol * abs(correct_mix_ice_r) << "\n";}
      if (abs(mix_liq_r[s] -  correct_mix_liq_r) >   cond_liq_r[s] * C::Tol * abs(correct_mix_liq_r)) {errors++;}
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
      saturation_tests(tmelt, sp(1e5), sp(610.7960763188032), sp(610.7960763188032),
         sp(0.003822318507864685),  sp(0.003822318507864685), errors);

      //Test vaules @ 243.15K @ 1e5 Pa
      saturation_tests(sp(243.15), sp(1e5), sp(37.98530141245404), sp(50.98455924912173),
         sp(0.00023634717905493638),  sp(0.0003172707211143376), errors);

      //Test values @ 303.15 @ 1e5 Pa
      saturation_tests(sp(303.15), sp(1e5), sp(4242.757341329608), sp(4242.757341329608),
        sp(0.0275579183092878), sp(0.0275579183092878), errors);

    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
};

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3Conservation
{


  KOKKOS_FUNCTION static void cloud_water_conservation_tests_device(){

    using KTH = KokkosTypes<HostDevice>;

    CloudWaterConservationData cwdc[1] = {{sp(1e-5), 0.0, sp(1.1), sp(1e-4), 0.0, 0.0, 0.0, 0.0, 0.0, sp(1.0), sp(1.0)}};

    // Sync to device
    KTH::view_1d<CloudWaterConservationData> cwdc_host("cwdc_host", 1);
    view_1d<CloudWaterConservationData> cwdc_device("cwdc_host", 1);

    // This copy only copies the input variables.
    std::copy(&cwdc[0], &cwdc[0] + 1, cwdc_host.data());
    Kokkos::deep_copy(cwdc_device, cwdc_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Spack qc(cwdc_device(0).qc);
      Spack qcnuc(cwdc_device(0).qcnuc);
      Spack qcaut(cwdc_device(0).qcaut);
      Spack qcacc(cwdc_device(0).qcacc);
      Spack qccol(cwdc_device(0).qccol);
      Spack qcheti(cwdc_device(0).qcheti);
      Spack qcshd(cwdc_device(0).qcshd);
      Spack qiberg(cwdc_device(0).qiberg);
      Spack qisub(cwdc_device(0).qisub);
      Spack qidep(cwdc_device(0).qidep);

      Functions::cloud_water_conservation(qc, qcnuc, cwdc_device(0).dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);

      cwdc_device(0).qc = qc[0];
      cwdc_device(0).qcnuc = qcnuc[0];
      cwdc_device(0).qcaut = qcaut[0];
      cwdc_device(0).qcacc = qcacc[0];
      cwdc_device(0).qccol = qccol[0];
      cwdc_device(0).qcheti = qcheti[0];
      cwdc_device(0).qcshd = qcshd[0];
      cwdc_device(0).qiberg = qiberg[0];
      cwdc_device(0).qisub = qisub[0];
      cwdc_device(0).qidep = qidep[0];
    });

    // Sync back to host
    Kokkos::deep_copy(cwdc_host, cwdc_device);

    const auto ratio = cwdc[0].qc/(cwdc[0].qcaut * cwdc[0].dt);
    REQUIRE(abs(cwdc_host(0).qcaut - cwdc[0].qcaut*ratio) <= C::Tol);
    REQUIRE(cwdc_host(0).qcacc == 0.0);
    REQUIRE(cwdc_host(0).qccol == 0.0);
    REQUIRE(cwdc_host(0).qcheti == 0.0);
    REQUIRE(cwdc_host(0).qcshd == 0.0);
    REQUIRE(cwdc_host(0).qiberg == 0.0);
    REQUIRE(abs(cwdc_host(0).qisub -(1.0 - ratio)) <= C::Tol);
    REQUIRE(abs(cwdc_host(0).qidep - (1.0 - ratio)) <= C::Tol);
    REQUIRE(cwdc_host[0].qcaut * cwdc[0].dt <= cwdc_host[0].qc);


  }


  KOKKOS_FUNCTION static void cloud_water_conservation_tests(int& errors){

    Spack qc(sp(1e-5));
    Spack qcnuc(0.0);
    Scalar dt = sp(1.1);
    Spack qcaut(sp(1e-4));
    Spack qcacc(0.0);
    Spack qccol(0.0);
    Spack qcheti(0.0);
    Spack qcshd(0.0);
    Spack qiberg(0.0);
    Spack qisub(sp(1.0));
    Spack qidep(sp(1.0));

    const auto ratio = qc/(qcaut * dt);
    Functions::cloud_water_conservation(qc, qcnuc, dt,
    qcaut,qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);


    //Here we check the case where sources > sinks
    if(((qcaut - sp(1.0e-4)*ratio) != 0.0).any()){errors++;}
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
    if((qcaut*dt != qc).any()){errors++;}

    // Check the case where sources > sinks with sinks = 0
    qcaut = 0.0;
    Functions::cloud_water_conservation(qc, qcnuc, dt,
    qcaut,qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);

    //In this case qidep and qisub should be set to zero
    if((qisub != 0.0).any()){errors++;}
    if((qidep != 0.0).any()){errors++;}

  }


  KOKKOS_FUNCTION static void rain_water_conservation_tests_device(){
     using KTH = KokkosTypes<HostDevice>;

     RainWaterConservationData rwdc[1] = {{sp(1e-5), 0.0, 0.0, 0.0, 0.0, sp(1.1), sp(1e-4), 0.0, 0.0 }};

     // Sync to device
     KTH::view_1d<RainWaterConservationData> rwdc_host("rwdc_host", 1);
     view_1d<RainWaterConservationData> rwdc_device("rwdc_host", 1);

     // This copy only copies the input variables.
     std::copy(&rwdc[0], &rwdc[0] + 1, rwdc_host.data());
     Kokkos::deep_copy(rwdc_device, rwdc_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      Spack qr(rwdc_device(0).qr);
      Spack qcaut(rwdc_device(0).qcaut);
      Spack qcacc(rwdc_device(0).qcacc);
      Spack qimlt(rwdc_device(0).qimlt);
      Spack qcshd(rwdc_device(0).qcshd);
      Spack qrevp(rwdc_device(0).qrevp);
      Spack qrcol(rwdc_device(0).qrcol);
      Spack qrheti(rwdc_device(0).qrheti);

      Functions::rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, rwdc_device(0).dt, qrevp, qrcol, qrheti);

      rwdc_device(0).qr = qr[0];
      rwdc_device(0).qcaut = qcaut[0];
      rwdc_device(0).qcacc = qcacc[0];
      rwdc_device(0).qimlt = qimlt[0];
      rwdc_device(0).qcshd = qcshd[0];
      rwdc_device(0).qrevp = qrevp[0];
      rwdc_device(0).qrcol = qrcol[0];
      rwdc_device(0).qrheti = qrheti[0];
    });


    // Sync back to host
    Kokkos::deep_copy(rwdc_host, rwdc_device);
    const auto ratio = rwdc[0].qr/(rwdc[0].qrevp * rwdc[0].dt);

    //Here we check cases where source > sinks and sinks > 1e-20
    REQUIRE(rwdc_host(0).qcaut == 0.0);
    REQUIRE(rwdc_host(0).qcacc == 0.0);
    REQUIRE(rwdc_host(0).qimlt == 0.0);
    REQUIRE(rwdc_host(0).qcshd == 0.0);

    //Check the value of qrevp
    REQUIRE(abs(rwdc_host(0).qrevp- rwdc[0].qrevp*ratio)<= C::Tol);

    //Now test that conservation has actually been enforced
    REQUIRE( rwdc_host(0).qr - rwdc_host(0).qrevp * rwdc_host(0).dt  ==  0.0);

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

    cloud_water_conservation_tests_device();

    rain_water_conservation_tests_device();
    //int nerr = 0;

    //TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    //Kokkos::parallel_reduce("ConservationTests", policy, KOKKOS_LAMBDA(const MemberType& tem, int& errors){

     // errors = 0;

      //cloud_water_conservation_tests(errors);

      //rain_water_conservation_tests(errors);

      //ice_water_conservation_tests(errors);

   // }, nerr);

   // Kokkos::fence();

   // REQUIRE(nerr==0);
  }

  static void cloud_water_conservation_unit_bfb_tests(){

    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    CloudWaterConservationData cwdc[max_pack_size] = {
      {9.9999999999999995e-7, 0.0, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 0.0, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 0.0, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 0.0, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 0.0, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 0.0, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 0.0, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 0.0, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 0.0, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7},

      {9.9999999999999995e-7, 0.0, 1800.0, 1.5832574016248739e-12, 1.0630996907148179e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {6.4285714285714288e-5, 0.0, 1800.0, 5.0577951315583066e-7, 7.7585489624948031e-4, 1.5683327213659326E-4, 1.2893174331809564e-14, 0.0, 5.0463073442953805e-6, 0.0, 5.1387602886199180e-7},
      {0.0, 0.0, 1800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {7.1428571428571434e-5, 0.0, 1800.0, 5.1480988828550771e-7, 7.7585489624948031e-4, 1.5597668529004373e-4, 4.9926620576534573e-14, 0.0, 6.7718890050008472e-6, 0.0, 7.1052455549903861e-7}
    };

    // Sync to device
    KTH::view_1d<CloudWaterConservationData> cwdc_host("cwdc_host", Spack::n);
    view_1d<CloudWaterConservationData> cwdc_device("cwdc_host", Spack::n);

    // This copy only copies the input variables.
    std::copy(&cwdc[0], &cwdc[0] + Spack::n, cwdc_host.data());
    Kokkos::deep_copy(cwdc_device, cwdc_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      cloud_water_conservation(cwdc[i]);
    }

    // This copy also copies the output from the fortran function into the host view. These values
    // are need to check the values returned from
    std::copy(&cwdc[0], &cwdc[0] + Spack::n, cwdc_host.data());

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      // Init pack inputs
      Spack qc, qcnuc, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep;
      for (Int s = 0; s < Spack::n; ++s) {
        qc[s] = cwdc_device(s).qc;
        qcnuc[s] = cwdc_device(s).qcnuc;
        qcaut[s] = cwdc_device(s).qcaut;
        qcacc[s] = cwdc_device(s).qcacc;
        qccol[s] = cwdc_device(s).qccol;
        qcheti[s] = cwdc_device(s).qcheti;
        qcshd[s] = cwdc_device(s).qcshd;
        qiberg[s] = cwdc_device(s).qiberg;
        qisub[s] = cwdc_device(s).qisub;
        qidep[s] = cwdc_device(s).qidep;
      }

      Functions::cloud_water_conservation(qc, qcnuc, cwdc_device(0).dt, qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);
      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        cwdc_device(s).qc = qc[s];
        cwdc_device(s).qcnuc = qcnuc[s];
        cwdc_device(s).qcaut = qcaut[s];
        cwdc_device(s).qcacc = qcacc[s];
        cwdc_device(s).qccol = qccol[s];
        cwdc_device(s).qcheti = qcheti[s];
        cwdc_device(s).qiberg = qiberg[s];
        cwdc_device(s).qisub = qisub[s];
        cwdc_device(s).qidep = qidep[s];
      }

    });
    // Sync back to host
    Kokkos::deep_copy(cwdc_host, cwdc_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(cwdc[s].qc == cwdc_host(s).qc);
      REQUIRE(cwdc[s].qcnuc == cwdc_host(s).qcnuc);
      REQUIRE(cwdc[s].qcaut == cwdc_host(s).qcaut);
      REQUIRE(cwdc[s].qcacc == cwdc_host(s).qcacc);
      REQUIRE(cwdc[s].qccol == cwdc_host(s).qccol);
      REQUIRE(cwdc[s].qcheti == cwdc_host(s).qcheti);
      REQUIRE(cwdc[s].qiberg == cwdc_host(s).qiberg);
      REQUIRE(cwdc[s].qisub == cwdc_host(s).qisub);
      REQUIRE(cwdc[s].qidep == cwdc_host(s).qidep);
    }

  }

  static void ice_water_conservation_unit_bfb_tests()
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
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0}
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

  static void rain_water_conservation_unit_bfb_tests(){

    using KTH = KokkosTypes<HostDevice>;

    static constexpr Int max_pack_size = 16;
    REQUIRE(Spack::n <= max_pack_size);

    RainWaterConservationData rwdc[max_pack_size] = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10},

      {0.0, 0.0, 0.0, 0.0, 0.0, 1800.0, 0.0, 0.0, 0.0},
      {3.6842105263157901e-6, 1.8910609577335389e-12, 6.5659507736611415e-9, 2.0267066625093075e-3, 1.3686661018890648e-9, 1800.0, 0.0, 0.0, 0.0},
      {1.0000000000000001e-5, 1.3239078166546396e-11, 4.5967389456540289e-8, 0.0, 0.0, 1800.0, 0.0, 1.4619847302347994e-33, 1.3104200383028957e-8},
      {8.9473684210526319e-6, 1.1338778389922441e-11, 3.9369360589471763e-8, 0.0, 0.0, 1800.0, 0.0, 1.4495908589465900e-33, 8.5051489557327688e-10}
    };

    // Sync to device
    KTH::view_1d<RainWaterConservationData> rwdc_host("rwdc_host", Spack::n);
    view_1d<RainWaterConservationData> rwdc_device("rwdc_host", Spack::n);

    // This copy only copies the input variables.
    std::copy(&rwdc[0], &rwdc[0] + Spack::n, rwdc_host.data());
    Kokkos::deep_copy(rwdc_device, rwdc_host);

    // Get data from fortran
    for (Int i = 0; i < max_pack_size; ++i) {
      rain_water_conservation(rwdc[i]);
    }

    // This copy also copies the output from the fortran function into the host view. These values
    // are need to check the values returned from
    std::copy(&rwdc[0], &rwdc[0] + Spack::n, rwdc_host.data());

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(RangePolicy(0, 1), KOKKOS_LAMBDA(const Int& i) {
      // Init pack inputs
      Spack qr, qcaut, qcacc, qimlt, qcshd, qrevp, qrcol, qrheti;
      for (Int s = 0; s < Spack::n; ++s) {
        qr[s] = rwdc_device(s).qr;
        qcaut[s] = rwdc_device(s).qcaut;
        qcacc[s] = rwdc_device(s).qcacc;
        qimlt[s] = rwdc_device(s).qimlt;
        qcshd[s] = rwdc_device(s).qcshd;
        qrevp[s] = rwdc_device(s).qrevp;
        qrcol[s] = rwdc_device(s).qrcol;
        qrheti[s] = rwdc_device(s).qrheti;
      }

      Functions::rain_water_conservation(qr, qcaut, qcacc, qimlt, qcshd, rwdc_device(0).dt, qrevp, qrcol, qrheti);
      // Copy results back into views
      for (Int s = 0; s < Spack::n; ++s) {
        rwdc_device(s).qr = qr[s];
        rwdc_device(s).qcaut = qcaut[s];
        rwdc_device(s).qcacc = qcacc[s];
        rwdc_device(s).qimlt = qimlt[s];
        rwdc_device(s).qcshd = qcshd[s];
        rwdc_device(s).qrevp = qrevp[s];
        rwdc_device(s).qrcol = qrcol[s];
        rwdc_device(s).qrheti = qrheti[s];
      }

    });

    // Sync back to host
    Kokkos::deep_copy(rwdc_host, rwdc_device);

    // Validate results
    for (Int s = 0; s < Spack::n; ++s) {
      REQUIRE(rwdc[s].qr == rwdc_host(s).qr);
      REQUIRE(rwdc[s].qcaut == rwdc_host(s).qcaut);
      REQUIRE(rwdc[s].qcacc == rwdc_host(s).qcacc);
      REQUIRE(rwdc[s].qimlt == rwdc_host(s).qimlt);
      REQUIRE(rwdc[s].qcshd == rwdc_host(s).qcshd);
      REQUIRE(rwdc[s].qrevp == rwdc_host(s).qrevp);
      REQUIRE(rwdc[s].qrcol == rwdc_host(s).qrcol);
      REQUIRE(rwdc[s].qrheti == rwdc_host(s).qrheti);
    }

  }

  static void run_bfb(){

      cloud_water_conservation_unit_bfb_tests();

      rain_water_conservation_unit_bfb_tests();

      ice_water_conservation_unit_bfb_tests();

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

TEST_CASE("p3_conservation_test", "[p3_conservation_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Conservation::run();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Conservation::run_bfb(); 
}

TEST_CASE("p3_cloud_water_autoconversion_test", "[p3_cloud_water_autoconversion_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_physics();
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3CloudWaterAutoconversion::run_bfb();
}

} // namespace
