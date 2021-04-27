
#include "catch2/catch.hpp"

#include "physics/share/physics_functions.hpp"
#include "physics/share/physics_universal_impl.hpp"
#include "physics_unit_tests_common.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

namespace scream {
namespace physics {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestUniversal
{

//-----------------------------------------------------------------------------------------------//
  KOKKOS_FUNCTION static void exner_tests(const Scalar& pressure, const Scalar& expected_val, int& errors){

    // Allow usage of universal functions
    using physics = scream::physics::Functions<Scalar, Device>;
    // Gather the machine epsilon for the error tolerance
    static constexpr Scalar eps = C::macheps;
    Real tol = 1000*eps;

    //========================================================
    // Test calculation of exners function given a specific pressure
    //========================================================
    // This function tests the "get_exner" universal physics function.  Comparing the output against
    // an expected exner result.
    //
    // Inputs:
    //   pressure:     A pressure value, Pa
    //   expected_val: What exners function should produce for the given pressure value, unitless.
    // Outputs:
    //   errors:       A tally of any errors that this test detects.
    //========================================================

    const Spack pres(pressure);
    const Spack exner = physics::get_exner(pres,Smask(true));

    if (std::abs(exner[0] - expected_val)>tol) {
      printf("exner_test: abs(exner-expected_exner)=%e is larger than the tol=%e\n",std::abs(exner[0] - expected_val),tol);
      errors++;
    }

  } // exner_tests

//-----------------------------------------------------------------------------------------------//
  KOKKOS_FUNCTION static void T_th_conversion_test(const Scalar& T_in, const Scalar& pres_in, int& errors){

    // Allow usage of universal functions
    using physics = scream::physics::Functions<Scalar, Device>;
    // Gather the test tolerance
    static constexpr Scalar eps = C::macheps;
    Real tol = 1000*eps;

    //========================================================
    // Test conversion of temperature to potential temperature, and vice versa
    //========================================================
    // This function tests both the T_to_th and th_to_T universal conversion functions.
    //
    // Inputs:
    //   T_in:    An example temperature, K.
    //   pres_in: A pressure value to use for conversion, Pa.
    // Outputs:
    //   errors:  A tally of any errors that this test detects.
    //========================================================

    const Spack T(T_in);
    const Spack pres(pres_in);
    // Retrieve exners function for this pressure level.  Note, get_exner is tested in a separate test.
    const Spack exner = physics::get_exner(pres,Smask(true));

    // Test conversion from temperature (T) to potential temperature (th)
    // Use T_in to convert from T to th
    // Note: T to th conversion is just a simple formula of th = T/exner
    const Spack th_atm = physics::T_to_th(T,exner,Smask(true));
    Real expected_th = T_in/exner[0];
    if (std::abs(th_atm[0]-expected_th)>tol) {
      printf("T to th test: abs(th_atm-expected_th)=%e is larger than the tol=%e\n",std::abs(th_atm[0]-expected_th),tol);
      errors++;
    }

    // Test conversion from potential temperature (th) to temperature (T)
    // Use T_in to convert from th to T, note, we use T_in again because we are just testing the formula and T_in should be sufficient for that.
    // Note: th to T conversion is just a simple formula of T = th*exner
    const Spack T_atm = physics::th_to_T(T,exner,Smask(true));
    Real expected_T = T_in*exner[0];
    if (std::abs(T_atm[0]-expected_T)>tol) {
      printf("T to th test: abs(th_atm-expected_th)=%e is larger than the tol=%e\n",std::abs(T_atm[0]-expected_T),tol);
      errors++;
    }

    // Test that T_to_th and th_to_T are inverses of each other.
    // Test T to th as inverses of each other
    // Converting T to th and then back to T should return the same value again.
    const Spack th_temp = physics::T_to_th(T,exner,Smask(true));
    const Spack T_new   = physics::th_to_T(th_temp,exner,Smask(true));
    if (std::abs(T_new[0]-T_in)>tol) {
      printf("T to th test: abs[ T_new (%.3e) - T_in (%.3e) ]=%e is larger than the tol=%e\n",T_new[0],T_in,std::abs(T_new[0]-T_in),tol);
      errors++;
    }
    // and vice versa
    const Spack T_temp = physics::th_to_T(T,exner,Smask(true));
    const Spack th_new = physics::T_to_th(T_temp,exner,Smask(true));
    if (std::abs(th_new[0]-T_in)>tol) {
      printf("T to th test: abs[ th_new (%.3e) - T_in (%.3e) ]=%e is larger than the tol=%e\n",th_new[0],T_in,std::abs(th_new[0]-T_in),tol);
      errors++;
    }
  } // T_th_conversion_test
//-----------------------------------------------------------------------------------------------//
  KOKKOS_FUNCTION static void dz_tests(const Scalar& zi_top, const Scalar& zi_bot, int& errors){

    // Allow usage of universal functions
    using physics = scream::physics::Functions<Scalar, Device>;
    // Gather the test tolerance
    static constexpr Scalar eps = C::macheps;
    Real tol = 1000*eps;

    //========================================================
    // Test calculation of layer thickness using get_dz
    //========================================================
    // This function tests the function "get_dz" for a set of interface layer heights.
    //
    // Inputs:
    //   zi_top: The above surface height of the top of the verical layer, m.
    //   zi_bot: The above surface height of the bottom of the verical layer, m.
    // Outputs:
    //   errors:  A tally of any errors that this test detects.
    //========================================================

    const Spack z_top(zi_top);
    const Spack z_bot(zi_bot);
    Real expected_dz = zi_top-zi_bot;
    const Spack dz = physics::get_dz(z_top, z_bot, Smask(true));
    if (std::abs(dz[0]-expected_dz)>tol) {
      printf("get_dz test: abs(dz-expected_dz)=%e is larger than the tol=%e\n",std::abs(dz[0]-expected_dz),tol);
      errors++;
    }

  } // dz_test
//-----------------------------------------------------------------------------------------------//
  KOKKOS_FUNCTION static void dse_tests(const Scalar& temp, const Scalar& height, const Scalar surface_height, int& errors){

    // Allow usage of universal functions
    using physics = scream::physics::Functions<Scalar, Device>;
    // Gather the test tolerance
    static constexpr Scalar eps = C::macheps;
    Real tol = 1000*eps;

    //========================================================
    // Test calculation of dry static energy using get_dse
    //========================================================
    // This function tests the function "get_dse".
    //
    // Inputs:
    //   temp is the atmospheric temperature. Units in K.
    //   height is the geopotential height above surface at midpoints. Units in m.
    //   surface_height is the surface geopotential height. Units in m.
    // Outputs:
    //   errors:  A tally of any errors that this test detects.
    //========================================================

    static constexpr Scalar cp = C::CP;
    static constexpr Scalar ggr = C::gravit;
    const Spack T_mid(temp);
    const Spack z_mid(height);
    Real expected_dse = cp*temp+ggr*height+surface_height;
    const Spack dse = physics::get_dse(T_mid, z_mid, surface_height, Smask(true));
    if (std::abs(dse[0]-expected_dse)>tol) {
      printf("get_dse test: abs(dse-expected_dse)=%e is larger than the tol=%e\n",std::abs(dse[0]-expected_dse),tol);
      errors++;
    }

  } // dse_test
//-----------------------------------------------------------------------------------------------//
  KOKKOS_FUNCTION static void virtual_temperature_test(const Scalar& T_mid_in, const Scalar& qv_in, int& errors){

    // Allow usage of universal functions
    using physics = scream::physics::Functions<Scalar, Device>;
    // Gather the machine epsilon for the error tolerance
    static constexpr Scalar eps = C::macheps;
    static constexpr Scalar ep_2 = C::ep_2;
    Real tol = 2000*eps;

    //========================================================
    // Test calculation of virtual temperature using get_virtual_temperature
    //========================================================
    // This function tests the function "get_virtual_temperature".
    //
    // Inputs:
    //   T_mid_in is the atmospheric temperature. Units in K.
    //   qv_in is the atmospheric water vapor mass mixing ratio.  Units in kg/kg
    // Outputs:
    //   errors:  A tally of any errors that this test detects.
    //========================================================
    const Spack T_mid(T_mid_in);
    const Spack qv(qv_in);
    // Given T and qv, determine T_virt
    const Spack T_virt = physics::get_virtual_temperature(T_mid,qv,Smask(true));
    // Determine the expected virtual temperature
    Real expected_T_virt = T_mid_in*(qv_in+ep_2)/(ep_2*(1.0+qv_in));

    if (std::abs(T_virt[0]-expected_T_virt)>tol) {
      printf("get_virtual_temperature test: abs(T_virt-expected_T_virt)=%e is larger than the tol=%e\n",std::abs(T_virt[0]-expected_T_virt),tol);
      errors++;
    }

  } // virtual_temperature_test
//-----------------------------------------------------------------------------------------------//
  static void run()
  {
    constexpr int num_levs = 100; // Number of levels to use for tests.

    // Compute random values for dse test
    using view_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Real>;
    view_1d temp("temp",num_levs),
            height("height",num_levs),
            surface_height("surface_height",num_levs),
            qv("qv",num_levs);

    std::random_device rd;
    using rngAlg = std::mt19937_64;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    rngAlg engine(seed);
    using RPDF = std::uniform_real_distribution<Real>;
    RPDF pdf(1e-3,1e3);

    ekat::genRandArray(temp,engine,pdf);
    ekat::genRandArray(height,engine,pdf);
    ekat::genRandArray(surface_height,engine,pdf);
    ekat::genRandArray(qv,engine,pdf);

    int nerr = 0;
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("test_universal_physics", policy, KOKKOS_LAMBDA(const MemberType&, int& errors) {

      errors += 0;

      static constexpr Scalar p0     = C::P0;
      static constexpr Scalar rd     = C::RD;
      static constexpr Scalar inv_cp = C::INV_CP;
      static constexpr Scalar tmelt  = C::Tmelt;

      // Create dummy level data for testing:
      Real pres_top = 200.;
      Real dp       = (p0-pres_top)/(num_levs-1);

      // Run tests
      for (int k=0;k<num_levs;++k)
      {
        // Exner test
        //   - Pressure at level k is just k*dp, so each pressure level is dp in thickness.
        Real pres = p0 + k*dp;
        Real p_val  = pow( pres/p0, rd*inv_cp);
        exner_tests(pres,p_val,errors);
        // T and TH conversion TEST
        Real T_atm = tmelt + k;
        T_th_conversion_test(T_atm,pres,errors);
        // DZ Test
        //   - Determine the z-layer bottom and top as being the (k+1) thick.
        //     Allow for varying interface heights by setting zi_bot equal to the sum(0,...,k).
        Real zi_bot = (pow(k,2)+k)/2.0;
        Real zi_top = zi_bot + (k+1);
        dz_tests(zi_top,zi_bot,errors);
        // DSE Test
        dse_tests(temp(k),height(k),surface_height(k),errors);
        // virtual temperature test
        virtual_temperature_test(temp(k),qv(k),errors);
      }

    }, nerr);

    Kokkos::fence();
    REQUIRE(nerr == 0);

  } // run
}; // end of TestUniversal struct

} // namespace unit_test
} // namespace physics
} // namespace scream

namespace{

TEST_CASE("physics_universal_test", "[physics_universal_test]"){
  scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUniversal::run();

 } // TEST_CASE

} // namespace
