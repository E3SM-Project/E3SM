#include "catch2/catch.hpp"

#include "physics/share/physics_functions.hpp"
#include "physics/share/physics_saturation_impl.hpp"
#include "physics_unit_tests_common.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace physics {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSaturation
{

  KOKKOS_FUNCTION static Scalar condNum(const Scalar& svp, const Scalar& temp, const bool isIce=false){

    //computes condition number for saturation vapor pressure calc.

    //Set constant values
    //--------------------------------------
    static constexpr Scalar RV = C::RV;
    static constexpr Scalar RhoLiq = C::RHO_H2O;
    static constexpr Scalar RhoIce = C::RhoIce;
    static constexpr Scalar LatVap = C::LatVap;
    static constexpr Scalar LatIce = C::LatIce;

    //==========================================================
    // Test Saturation Vapor Pressure
    //==========================================================
    // First calculate condition number Cond=x*f'(x)/f(x), which gives the growth in roundoff error
    // expected from the function. Here, x=temperature, f(x) = sat_X_p, and
    // f'(x) = L/T*sat_X_p/(Rv*T - sat_X_p/rho_{liq or ice}) from Clausius-Clapeyron.
    // Note this analytical solution isn't exact because it misses things like surface tension
    // so it can't replace curve fits like Flatau. But Clausius-Clapeyron is good enough for
    // getting a ballpark value, which is all we're doing here.

    //for liquid
    auto latent = LatVap;
    auto rho    = RhoLiq;

    if(isIce){
      //for ice
      latent += LatIce;
      rho     = RhoIce;
    }

    return latent*svp/(RV*temp - svp/rho); //return condition number

  }

  KOKKOS_FUNCTION  static void saturation_tests(const Scalar& temperature, const Scalar& pressure,
          const Scalar *expected_vals, int& errors ){

    //Nomenclature:
    //subscript "_fp"  stands for "Flatau Pressure"
    //subscript "_fr"  stands for "Flatau mixing Ratios"
    //subscript "_mkp" stands for "Murphy Koop Pressure"
    //subscript "_mkr" stands for "Murphy Koop mixing Ratios"

    //extract values from the array
    const auto expected_sat_ice_fp = expected_vals[0];
    const auto expected_sat_liq_fp = expected_vals[1];
    const auto expected_mix_ice_fr = expected_vals[2];
    const auto expected_mix_liq_fr = expected_vals[3];

    const auto expected_sat_ice_mkp = expected_vals[4];
    const auto expected_sat_liq_mkp = expected_vals[5];
    const auto expected_mix_ice_mkr = expected_vals[6];
    const auto expected_mix_liq_mkr = expected_vals[7];

    //Allow usage of saturation functions
    using physics = scream::physics::Functions<Scalar, Device>;

    //Convert Scalar inputs to Spacks because that's what polysvp1 and qv_sat expect as inputs.
    //--------------------------------------
    const Spack temps(temperature);
    const Spack pres(pressure);

    //Get values from polysvp1 and qv_sat (qv_sat calls polysvp1 here) to test against "expected" values
    //--------------------------------------
    const Spack sat_ice_fp  = physics::polysvp1(temps, true, Smask(true));
    const Spack sat_liq_fp  = physics::polysvp1(temps, false, Smask(true));
    //last argument "0" of qv_sat function below forces qv_sat to call "polysvp1"
    const Spack mix_ice_fr = physics::qv_sat(temps, pres, true, Smask(true), physics::Polysvp1);
    const Spack mix_liq_fr = physics::qv_sat(temps, pres, false,Smask(true), physics::Polysvp1);

    //Get values from MurphyKoop_svp and qv_sat (qv_sat calls MurphyKoop_svp here) to test against "expected" values
    const Spack sat_ice_mkp   = physics::MurphyKoop_svp(temps, true, Smask(true));
    const Spack sat_liq_mkp   = physics::MurphyKoop_svp(temps, false, Smask(true));
    //last argument "1" of qv_sat function below forces qv_sat to call "MurphyKoop_svp"
    const Spack mix_ice_mkr  = physics::qv_sat(temps, pres, true,  Smask(true), physics::MurphyKoop);
    const Spack mix_liq_mkr  = physics::qv_sat(temps, pres, false, Smask(true), physics::MurphyKoop);

    //Set error tolerances
    //--------------------------------------
    //: C::macheps is machine epsilon for single or double precision as appropriate. This will be
    //multiplied by a condition # to get the actual expected numerical uncertainty.

    static constexpr Scalar tol = C::macheps
#ifdef NDEBUG
      * 10000
#endif
;

    //PMC note: original version looped over pack dimension, testing each entry. This isn't
    //necessary b/c packs were created by copying a scalar up to pack size. Thus just evaluating
    // 1st entry below.

    const Scalar Cond_ice_fp = condNum(sat_ice_fp[0],temps[0], true); //"isIce=true"
    const Scalar Cond_liq_fp = condNum(sat_liq_fp[0],temps[0]);

    const Scalar Cond_ice_mkp = condNum(sat_ice_mkp[0],temps[0], true);//"isIce=true"
    const Scalar Cond_liq_mkp = condNum(sat_liq_mkp[0],temps[0]);

    // Test vapor pressure against Flatau's impl of Wexler:
    // ---------------------------------------------------------
    // Now check that computed vs expected values are small enough.
    if ( std::abs(sat_ice_fp[0] - expected_sat_ice_fp ) > Cond_ice_fp*tol ) {
      printf("esi_fp for T = %f abs diff is %e but max allowed is %e\n",
             temperature,std::abs(sat_ice_fp[0] - expected_sat_ice_fp ),tol*Cond_ice_fp );
      errors++;
    }
    if (std::abs(sat_liq_fp[0] - expected_sat_liq_fp) > Cond_liq_fp*tol)  {
      printf("esl_fp  for T = %f abs diff is %e but max allowed is %e\n",
             temperature,std::abs(sat_liq_fp[0] - expected_sat_liq_fp ),tol*Cond_liq_fp);
      errors++;
    }

    // Test vapor pressure against Murphy and Koop:
    // ---------------------------------------------------------
    // Now check that computed vs expected values are small enough.
    if ( std::abs(sat_ice_mkp[0] - expected_sat_ice_mkp ) > Cond_ice_mkp*tol ) {
      printf("esi_mkp for T = %f abs diff is %e but max allowed is %e\n",
             temperature,std::abs(sat_ice_mkp[0] - expected_sat_ice_mkp ),tol*Cond_ice_mkp );
      errors++;
    }
    if (std::abs(sat_liq_mkp[0] - expected_sat_liq_mkp) > Cond_liq_mkp*tol)  {
      printf("esl_mkp  for T = %f abs diff is %e but max allowed is %e\n",
             temperature,std::abs(sat_liq_mkp[0] - expected_sat_liq_mkp ),tol*Cond_liq_mkp);
      errors++;
    }

    //Set constant values
    //--------------------------------------
    static constexpr Scalar RV = C::RV;
    static constexpr Scalar LatVap = C::LatVap;
    static constexpr Scalar LatIce = C::LatIce;

    //==========================================================
    // Test Saturation Mixing Ratio
    //==========================================================
    // First, compute condition # Cond=x*f'(x)/f(x), with x=temperature, f(x)=mix_X_r[0], and
    // f'(x) = L*mix_X_r[0]/(Rv*temperature**2.) from Clausius-Clapeyron. Nice cancelation leaves:

    const Scalar Cond_ice_r=(LatVap+LatIce)/(RV*temps[0]);
    const Scalar Cond_liq_r=LatVap/(RV*temps[0]);

    //Test mixing-ratios against Wexler approx:
    // -------------------
    // Now check that computed vs expected values are small enough (Flatau).
    if (std::abs(mix_ice_fr[0] -  expected_mix_ice_fr) > Cond_ice_r*tol ) {
      printf("qsi_fp: abs(calc-expected)=%e %e\n",std::abs(mix_ice_fr[0] -  expected_mix_ice_fr),tol*Cond_ice_r);
      errors++;
    }
    if (std::abs(mix_liq_fr[0] -  expected_mix_liq_fr) > Cond_liq_r*tol ) {
      printf("qsl_fp: abs(calc-expected)=%e %e\n",std::abs(mix_liq_fr[0] -  expected_mix_liq_fr),tol*Cond_liq_r);
      errors++;
    }

    // Now check that computed vs expected values are small enough (Murphy and Koop).
    if (std::abs(mix_ice_mkr[0] -  expected_mix_ice_mkr) > Cond_ice_r*tol ) {
      printf("qsi_mkp: abs(calc-expected)=%e %e\n",std::abs(mix_ice_mkr[0] -  expected_mix_ice_mkr),tol*Cond_ice_r);
      errors++;
    }
    if (std::abs(mix_liq_mkr[0] -  expected_mix_liq_mkr) > Cond_liq_r*tol ) {
      printf("qsl_mkp: abs(calc-expected)=%e %e\n",std::abs(mix_liq_mkr[0] -  expected_mix_liq_mkr),tol*Cond_liq_r);
      errors++;
    }
  }

  static void run()
  {
    /*Originally written by Kyle Pressel, updated by Peter Caldwell on 4/5/20.
     *This code tests polysvp1 and qv_sat at 0 degrees C, at a very cold T, and at a very hot T
     *to make sure our impl gets the same answer as Flatau et al 1992:
     *(https://journals.ametsoc.org/jamc/article/31/12/1507/14870/Polynomial-Fits-to-Saturation-Vapor-Pressure)
     *For 0 degrees C, polysvp values can be read directly from Flatau. For other cases, I independently
     *coded up the Flatau scheme (polysvp1) in python and used it to derive the expected values. My python code is
     *in https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py
     */

    int nerr = 0;
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("TestTableIce::run", policy, KOKKOS_LAMBDA(const MemberType&, int& errors) {

      errors = 0;

      //Following struct stores data for lauching the saturation_tests function
      //for polysvp1 and MurphyKoop(MK) saturation schemes. We are supplying temperature,
      //pressure, expected values in double and single precision. Polysvp1 is very
      //sensitive (due to shear number of calcs required at low temps)to single
      //precision values at low temperatures, therefore we explicity supply single
      //precision expected values. Single precision expected values are computed
      //from a python code using numpy's float32. The code uses the same functions
      //for polysvp1 and MK as in the code present at:
      // https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py

      //total number expected vals to store for each precision for both MK and polysvp1
      static constexpr int ncrv = 8; // 4 for each saturation scheme

      struct sat_test_args {
        Scalar t_atm; //temperature [K]
        Scalar pres; //pressure [Pa]
        Scalar dp_data[ncrv]; // double precision expected vals for MK and polysvp1
        Scalar sp_data[ncrv]; // single precision expected vals for MK and polysvp1
      };

      static constexpr auto ncases   = 10; //total number of test cases

      sat_test_args stargs[ncases]; // variable to store all the data for launching tests

      // Just Freezing Case: Test values @ 273.15K (tmelt) @ 1e5 Pa
      //---------------------------------------
      // This is the nicest test b/c polysvp1 is a polynomial fit around (T-273.15)
      // so T=273.15 collapses back to the intercept coefficient which can be read
      // directly from the RHS of table 4 of Flatau et al 1992.
      // Note that ice values are identical to liquid values b/c C++ uses liq val for T>=0 C.

      static constexpr auto tmelt = C::Tmelt;
      static constexpr auto atm_pres = 1e5;

      stargs[0] = {tmelt,
                   atm_pres,
                   //dp_data
                  {611.23992100000009, 611.23992100000009,
                   0.0038251131382843278, 0.0038251131382843278,
                   611.2126978267946, 611.2126978267946,
                   0.0038249417291628678, 0.0038249417291628678},
                   //sp_data
                  {611.2399, 611.2399, 0.0038251132, 0.0038251132,
                   611.2123, 611.2123, 0.003824939, 0.003824939}
      };

      // Cold Case: Test values @ 243.15K @ 1e5 Pa
      //---------------------------------------
      stargs[1] = {243.15,
                   atm_pres,
                   //dp_data
                  {38.024844602056795, 51.032583257624964,
                   0.00023659331311441935, 0.00031756972127516819,
                   38.01217277745647, 50.93561537896607,
                   0.00023651443812988484, 0.0003169659941390894},
                   //sp_data
                  {38.024666, 51.03264, 0.00023659221, 0.00031757005,
                   38.012142, 50.935585, 0.00023651426, 0.0003169658}
      };

      //Warm Case: Test values @ 303.15K @ 1e5 Pa
      //---------------------------------------
      stargs[2] = {303.15,
                   atm_pres,
                   //dp_data
                  {4245.1933273786717, 4245.1933273786717,
                   0.027574442204332306, 0.027574442204332306,
                   4246.814076877233, 4246.814076877233,
                   0.027585436614272162, 0.027585436614272162},
                   //sp_data
                  {4245.1934, 4245.1934, 0.027574444, 0.027574444,
                   4246.8115, 4246.8115, 0.027585419, 0.027585419}
      };

      //Following values are picked from Murphy and Koop (2005)
      //Table C1 titled: "VALUES RECOMMENDED FOR CHECKING COMPUTER CODES"
      //Saturation vapor pressure (SVP) values in the table were upto only 5 significant digits.
      //Python code at https://github.com/E3SM-Project/scream-docs.git analysis-scripts/test_qv_sat.py
      //was extended to print Murphy and Koop SVP. "expected values" below are from that python code
      //for both polysvp1 and MK. Python code's computed "expected values" are exactly
      //same as compared to the values in MK table for upto 5 significant digits.

      //Test values @ 150K @ 1e5 Pa
      stargs[3] = {150,
                   atm_pres,
                   //dp_data
                  {0.0565113360640801, 0.17827185346988017,
                   3.514840868181163e-07, 1.1088004687434155e-06,
                   6.1061006509816675e-06, 1.5621037177920032e-05,
                   3.79781500154674e-11, 9.715825652191167e-11},
                   //sp_data
                  {0.05517006, 0.17905235, 3.4314175e-07, 1.113655e-06,
                   6.1060973e-06, 1.5621057e-05, 3.797813e-11, 9.715838e-11}
      };

      //Test values @ 180K @ 1e5 Pa
      stargs[4] = {180,
                   atm_pres,
                   //dp_data
                  {0.0565113360640801, 0.17827185346988017,
                   3.514840868181163e-07, 1.1088004687434155e-06,
                   0.005397500125274297, 0.01123923029036248,
                   3.3570864981545485e-08, 6.990471437859482e-08},
                   //sp_data
                  {0.05517006, 0.17905235, 3.4314175e-07, 1.113655e-06,
                   0.0053975, 0.011239249, 3.3570867e-08, 6.990483e-08}
      };

      //Test values @ 210K  @ 1e5 Pa
      stargs[5] = {210,
                   atm_pres,
                   //dp_data
                  {0.7021857060894199, 1.2688182238880685,
                   4.3674192198021676e-06, 7.891776277281936e-06,
                   0.7020234713180218, 1.2335424085746476,
                   4.366410153074117e-06, 7.672365591576349e-06},
                   //sp_data
                  {0.7016659, 1.2685776, 4.364186e-06, 7.890279e-06,
                   0.70202345, 1.2335398, 4.36641e-06, 7.67235e-06}
      };

      //Test values @ 240K @ 1e5 Pa
      stargs[6] = {240,
                   atm_pres,
                   //dp_data
                  {27.280908658710246, 37.77676490183603,
                   0.00016972553017335565, 0.0002350491600776575,
                   27.272365420780556, 37.66700070557609,
                   0.00016967236474679822, 0.0002343659437158037},
                   //sp_data
                  {27.28076, 37.776802, 0.0001697246, 0.00023504937,
                   27.27235, 37.666973, 0.00016967228, 0.00023436577}
      };

      //Test values @ 273.16K @ 1e5 Pa
      stargs[7] = {273.16,
                   atm_pres,
                   //dp_data
                  {611.6840516537769, 611.6840516537769,
                   0.003827909594290528, 0.003827909594290528,
                   611.6570436443282, 611.6570436443282,
                   0.0038277395384149105, 0.0038277395384149105},
                   //sp_data
                  {611.68445, 611.68445, 0.0038279123, 0.0038279123,
                   611.65607, 611.65607, 0.0038277335, 0.0038277335}
      };
      //Test values @ 300K @ 1e5 Pa
      stargs[8] = {300,
                   atm_pres,
                   //dp_data
                  {3535.4066341569387, 3535.4066341569387,
                   0.022795088436007804, 0.022795088436007804,
                   3536.7644130514645, 3536.7644130514645,
                   0.022804163906259393, 0.022804163906259393},
                   //sp_data
                  {3535.4077, 3535.4077, 0.022795096, 0.022795096,
                   3536.7559, 3536.7559, 0.022804108, 0.022804108}
      };

      //Change Pressure Case: Test values @ 243.15 @ 500 mb
      //---------------------------------------
      stargs[9] = {243.15,
                   5e4,
                   //dp_data
                  {38.024844602056795, 51.032583257624964,
                   0.00047336669164733106, 0.00063546390177500586,
                   38.01217277745647, 50.93561537896607,
                   0.00047320882161578, 0.0006342552147122389},
                   //sp_data
                  {38.024666, 51.03264, 0.00047336446, 0.00063546456,
                   38.012142, 50.935585, 0.00047320846, 0.00063425483}
      };

      //Launch Tests:
      //---------------------------------------------

      //find out which precision to use
      constexpr auto is_single_prec = ekat::is_single_precision<Scalar>::value;

      if(is_single_prec){
        for (size_t ista = 0; ista < ncases; ista++) {
          //use sp_data values
          saturation_tests(stargs[ista].t_atm,stargs[ista].pres,stargs[ista].sp_data,errors);
        }
      } else{
        for (size_t ista = 0; ista < ncases; ista++) {
          //use dp_data values
          saturation_tests(stargs[ista].t_atm,stargs[ista].pres,stargs[ista].dp_data,errors);
        }
      }


    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);
  }
}; //end of TestSaturation struct

} // namespace unit_test
} // namespace physics
} // namespace scream

namespace{

TEST_CASE("physics_saturation_test", "[physics_saturation_test]"){
  scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSaturation::run();

 } // TEST_CASE

} // namespace
