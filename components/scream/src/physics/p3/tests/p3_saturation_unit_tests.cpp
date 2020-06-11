#include "catch2/catch.hpp"

#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/common/physics_functions.hpp"
#include "physics/common/physics_saturation_impl.hpp"

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

    //Allow usage of saturation functions
    using physics = scream::physics::Functions<Scalar, Device>;
    
    //Convert Scalar inputs to Spacks because that's what polysvp1 and qv_sat expect as inputs.
    //--------------------------------------
    const Spack temps(temperature);
    const Spack pres(pressure);

    //Get values from polysvp1 and qv_sat to test against "correct" values
    //--------------------------------------
    Spack sat_ice_p = physics::polysvp1(temps, true);
    Spack sat_liq_p = physics::polysvp1(temps, false);
    Spack mix_ice_r = physics::qv_sat(temps, pres, true);
    Spack mix_liq_r = physics::qv_sat(temps, pres, false);

    //Set error tolerances
    //--------------------------------------
    //: C::Tol is machine epsilon for single or double precision as appropriate. This will be
    //multiplied by a condition # to get the actual expected numerical uncertainty. If single precision, multiplying
    //by an additional fudge factor because the python values we're comparing against were computed via a
    //bunch of double-precision intermediate calculations which will add error.

    Scalar tol = (util::is_single_precision<Scalar>::value ) ? C::Tol*2 : C::Tol;

    //Set constant values
    //--------------------------------------
    const Scalar RV = C::RV;
    const Scalar RhoLiq = C::RhoH2O;
    const Scalar RhoIce = C::RhoIce;
    const Scalar LatVap = C::LatVap;
    const Scalar LatIce = C::LatIce;

    //PMC note: original version looped over pack dimension, testing each entry. This isn't
    //necessary b/c packs were created by copying a scalar up to pack size. Thus just evaluating
    // 1st entry below.

    //==========================================================
    // Test Saturation Vapor Pressure
    //==========================================================
    // First calculate condition number Cond=x*f'(x)/f(x), which gives the growth in roundoff error
    // expected from the function. Here, x=temperature, f(x) = sat_X_p, and
    // f'(x) = L/T*sat_X_p/(Rv*T - sat_X_p/rho_{liq or ice}) from Clausius-Clapeyron.
    // Note this analytical solution isn't exact because it misses things like surface tension
    // so it can't replace curve fits like Flatau. But Clausius-Clapeyron is good enough for
    // getting a ballpark value, which is all we're doing here.

    const Scalar Cond_ice_p=(LatVap+LatIce)*sat_ice_p[0]/(RV*temps[0] - sat_ice_p[0]/RhoIce);
    const Scalar Cond_liq_p=LatVap*sat_liq_p[0]/(RV*temps[0] - sat_liq_p[0]/RhoLiq);

    // Test vapor pressure against Flatau's impl of Wexler:
    // ---------------------------------------------------------
    // Now check that computed vs expected values are small enough.
    if ( std::abs(sat_ice_p[0] - correct_sat_ice_p ) > Cond_ice_p*tol ) {
      printf("esi for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_ice_p[0] - correct_sat_ice_p ),tol*Cond_ice_p );
      errors++;}
    if (std::abs(sat_liq_p[0] - correct_sat_liq_p) > Cond_liq_p*tol)  {
      printf("esl  for T = %f abs diff is %e but max allowed is %e\n",temperature,std::abs(sat_liq_p[0] - correct_sat_liq_p ),tol*Cond_liq_p);
      errors++;}

    //==========================================================
    // Test Saturation Mixing Ratio
    //==========================================================
    // First, compute condition # Cond=x*f'(x)/f(x), with x=temperature, f(x)=mix_X_r[0], and
    // f'(x) = L*mix_X_r[0]/(Rv*temperature**2.) from Clausius-Clapeyron. Nice cancelation leaves:

    const Scalar Cond_ice_r=(LatVap+LatIce)/(RV*temps[0]);
    const Scalar Cond_liq_r=LatVap/(RV*temps[0]);

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

