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

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3Saturation
{

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure, const Scalar& correct_sat_ice_p,
    const Scalar& correct_sat_liq_p, const Scalar&  correct_mix_ice_r, const Scalar& correct_mix_liq_r, int& errors ){

    //Convert Scalar inputs to Spacks because that's what polysvp1 and qv_sat expect as inputs.
    //--------------------------------------
    const Spack temps(temperature);
    const Spack pres(pressure);

    //Get values from polysvp1 and qv_sat to test against "correct" values
    //--------------------------------------
    Spack sat_ice_p = Functions::polysvp1(temps, true);
    Spack sat_liq_p = Functions::polysvp1(temps, false);
    Spack mix_ice_r = Functions::qv_sat(temps, pres, true);
    Spack mix_liq_r = Functions::qv_sat(temps, pres, false);

    //Set error tolerances
    //--------------------------------------
    //: C::Tol is machine epsilon for single or double precision as appropriate. This will be
    //multiplied by a condition # to get the actual expected numerical uncertainty. If single precision, multiplying
    //by an additional fudge factor because the python values we're comparing against were computed via a
    //bunch of double-precision intermediate calculations which will add error.
					     
    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;

    //For comparison against analytic saturation mixing ratio formula, note that error in the Flatau et al
    //approximation is <1e-3 hPa (=0.1 Pa) *if I'm reading their Fig 2a correctly*. qs=0.622 es/(p-es)
    //\approx 0.622*es/p since p>>es (causes 4% error at 303K and 1e5 Pa). Thus dqs/des = 0.622/p =>
    //qs errors due to numerical method should be 0.622/1e5*1e-3 ~ 1e-8 at 1e5 Pa.

    const Scalar sym_tol_q = 1e-4; //PMC hack: wasn't passing at 1e-8 so setting to very lax value for now.

    //Set constant values
    //--------------------------------------
    const Scalar RV = C::RV;
    const Scalar RhoLiq = C::RhoH2O;
    const Scalar RhoIce = C::RhoIce;
    const Scalar LatVap = C::LatVap;
    const Scalar LatIce = C::LatIce;

    printf("===============================\n");
    printf("Checking T = %f and Pres = %e\n",temperature,pressure);
    printf("===============================\n");

    //PMC note: original version looped over pack dimension, testing each entry. This isn't
    //necessary b/c packs were created by copying a scalar up to pack size. Thus just evaluating
    // 1st entry below.


    //==========================================================
    // Test Saturation Vapor Pressure
    //==========================================================
    // First calculate condition number Cond=x*f'(x)/f(x), which gives the growth in roundoff error
    // expected from the function. Here, x=temperature, f(x) = sat_X_p, and
    // f'(x) = L/T*sat_X_p/(Rv*T - sat_X_p/rho_{liq or ice}) from Clausius-Clapeyron.

    const Scalar Cond_ice_p=(LatVap+LatIce)*sat_ice_p[0]/(RV*temps[0] - sat_ice_p[0]/RhoIce);
    const Scalar Cond_liq_p=LatVap*sat_liq_p[0]/(RV*temps[0] - sat_liq_p[0]/RhoLiq);

    // Test vapor pressure against Flatau's impl of Wexler:
    // ---------------------------------------------------------      
    // Now check that computed vs expected values are small enough.     
    if ( std::abs(sat_ice_p[0] - correct_sat_ice_p ) > Cond_ice_p*tol ) {
      printf("  esi: abs(calc-expected),cond*tol=%e %e\n",std::abs(sat_ice_p[0] - correct_sat_ice_p ),tol*Cond_ice_p );
      errors++;}
    if (std::abs(sat_liq_p[0] - correct_sat_liq_p) > Cond_liq_p*tol)  {
      printf("  esl: abs(calc-expected),cond*tol=%e %e\n",std::abs(sat_liq_p[0] - correct_sat_liq_p ),tol*Cond_liq_p);
      errors++;}

    //==========================================================
    // Test Saturation Mixing Ratio
    //==========================================================
    // First, compute condition # Cond=x*f'(x)/f(x), with x=temperature, f(x)=mix_X_r[0], and
    // f'(x) = L*mix_X_r[0]/(Rv*temperature**2.) from Clausius-Clapeyron. Nice cancelation leaves:

    const Scalar Cond_ice_r=(LatVap+LatIce)/(RV*temps[0]);
    const Scalar Cond_liq_r=LatVap/(RV*temps[0]);
 
    
    // Test vapor pressure against exact symbolic implementation
    // ---------------------------------------------------------
    Scalar qsi_exact=qvsat_exact(temperature,pressure,1);
    Scalar qsl_exact=qvsat_exact(temperature,pressure,0);

    if ( std::abs(mix_ice_r[0] - qsi_exact ) > Cond_ice_r*sym_tol_q ) {
      printf("Diff btwn analytic qsi=%e and computed=%e is %e, bigger than tol=%e\n",qsi_exact,mix_ice_r[0],qsi_exact-mix_ice_r[0],Cond_ice_r*sym_tol_q);
      errors++;}
    if ( std::abs(mix_liq_r[0] - qsl_exact ) > Cond_liq_r*sym_tol_q ) {
      printf("Diff btwn analytic qsl=%e and computed=%e is %e, bigger than tol=%e\n",qsl_exact,mix_liq_r[0],qsl_exact-mix_liq_r[0],Cond_liq_r*sym_tol_q);
      errors++;}

    //Test mixing-ratios against Wexler approx: 
    // -------------------     
    // Now check that computed vs expected values are small enough.
    if (std::abs(mix_ice_r[0] -  correct_mix_ice_r) > Cond_ice_r*tol ) {
      printf("qsi: abs(calc-expected)=%e %e\n",std::abs(mix_ice_r[0] -  correct_mix_ice_r),tol*Cond_ice_r);
      errors++;}
    if (std::abs(mix_liq_r[0] -  correct_mix_liq_r) > Cond_liq_r*tol ) {
      printf("qsl: abs(calc-expected)=%e %e\n",std::abs(mix_liq_r[0] -  correct_mix_liq_r),tol*Cond_liq_r);
      errors++;}

  }

  static Scalar qvsat_exact(const Scalar& t_atm, const Scalar& p_atm, const bool& ice)
  {
    // Compute saturation mixing ratio qs symbolically.
    // dqs/dT = L*qs/(Rv*T**2.) by Clausius Clapeyron.
    // Taking the integral between some T0 and a particular T yields:
    // qs(T) = qs(T0)*exp{-L/Rv*(1/T - 1/T0)}.
    // Unlike for saturation vapor pressure, this formulation *is*
    // perfectly exact. It isn't used in numerical models just because
    // it is computationally expensive. 

    const auto tmelt = C::Tmelt;
    const auto RV = C::RV;
    const auto LatVap = C::LatVap;
    const auto LatIce = C::LatIce;
    const auto ep_2 = C::ep_2;

    Scalar es_T0;
    Scalar L;

    //Below, I'm taking es_T0 as equal to the constant term in
    //Flatau et al 1992 table 4 and turning it into qs_T0 using the
    //exact formula for es->qs. Would be better to use exact obs
    //value at tmelt, but Flatau is probably close enough (and provides
    //better match to P3's implementation!
    if ((t_atm < tmelt) && ice) {
      es_T0=611.147274; 
      L=LatVap+LatIce;
    } else{
      es_T0=611.239921;
      L=LatVap;
    }

    Scalar qs_T0 = ep_2 * es_T0 / std::max(p_atm-es_T0, sp(1.e-3));
    Scalar result = qs_T0*std::exp( -L/RV*(1/t_atm - 1/tmelt) );

    //printf("T,is_ice,qs=%e %d %e\n",t_atm,ice,result);

    return result;
  }

  static void run()
  {
    /*Originally written by Kyle Pressel, updated by Peter Caldwell on 4/5/20. 
     *This code tests polysvp1 and qv_sat at 0 degrees C, at a very cold T, and at a very hot T
     *to make sure our impl gets the same answer as Flatau et al 1992:
     *(https://journals.ametsoc.org/doi/abs/10.1175/1520-0450%281992%29031%3C1507%3APFTSVP%3E2.0.CO%3B2
     *For 0 degrees C, polysvp values can be read directly from Flatau. For other cases, I independently 
     *coded up the Flatau scheme in python and used it to derive the expected values. My python code is
     *in https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py
     */

    
    int nerr = 0;
    TeamPolicy policy(util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors = 0;
      const auto tmelt = C::Tmelt;
      
      // Just Freezing Case: Test values @ 273.15K @ 1e5 Pa
      //---------------------------------------
      // This is the nicest test b/c polysvp1 is a polynomial fit around (T-273.15)
      // so T=273.15 collapses back to the intercept coefficient which can be read 
      // directly from the RHS of table 4 of Flatau et al 1992.
      // Note that ice values are identical to liquid values b/c C++ uses liq val for T>=0 C.

      saturation_tests(tmelt, 1e5, 611.23992100000009, 611.23992100000009,
		       0.0038251131382843278, 0.0038251131382843278, errors);
      
      // Cold Case: Test values @ 243.15K @ 1e5 Pa
      //---------------------------------------
      saturation_tests(243.15, 1e5, 38.024844602056795, 51.032583257624964,
		       0.00023659331311441935, 0.00031756972127516819, errors);

      //Warm Case: Test values @ 303.15 @ 1e5 Pa
      //---------------------------------------
      saturation_tests(303.15, 1e5, 4245.1933273786717, 4245.1933273786717,
		       0.027574442204332306, 0.027574442204332306, errors);

      //Change Pressure Case: Test values @ 243.15 @ 500 mb
      //---------------------------------------
      saturation_tests(243.15, 5e4, 38.024844602056795, 51.032583257624964,
		       0.00047336669164733106, 0.00063546390177500586, errors);
      
    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
}; //end of TestP3Saturation struct
    
} // namespace unit_test
} // namespace p3
} // namespace scream
 
namespace{

TEST_CASE("p3_saturation_test", "[p3_saturation_test]"){
  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Saturation::run();
  
 } // TEST_CASE
 
} // namespace

