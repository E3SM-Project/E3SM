#include "catch2/catch.hpp"

#include "physics/share/physics_functions.hpp"
#include "physics/share/common_physics_functions.hpp"
#include "physics/share/physics_universal_impl.hpp"
#include "physics/share/common_physics_impl.hpp"
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
//ASD  KOKKOS_FUNCTION static void exner_tests(const Scalar& pressure, const Scalar& expected_val, int& errors){
//ASD
//ASD    // Allow usage of universal functions
//ASD    using physics = scream::physics::Functions<Scalar, Device>;
//ASD    // Gather the machine epsilon for the error tolerance
//ASD    static constexpr Scalar eps = C::macheps;
//ASD    Real tol = 1000*eps;
//ASD
//ASD    //========================================================
//ASD    // Test calculation of exners function given a specific pressure
//ASD    //========================================================
//ASD    // This function tests the "get_exner" universal physics function.  Comparing the output against
//ASD    // an expected exner result.
//ASD    //
//ASD    // Inputs:
//ASD    //   pressure:     A pressure value, Pa
//ASD    //   expected_val: What exners function should produce for the given pressure value, unitless.
//ASD    // Outputs:
//ASD    //   errors:       A tally of any errors that this test detects.
//ASD    //========================================================
//ASD
//ASD    const Spack pres(pressure);
//ASD    const Spack exner = physics::get_exner(pres,Smask(true));
//ASD
//ASD    if (std::abs(exner[0] - expected_val)>tol) {
//ASD      printf("exner_test: abs(exner-expected_exner)=%e is larger than the tol=%e\n",std::abs(exner[0] - expected_val),tol);
//ASD      errors++;
//ASD    }
//ASD
//ASD  } // exner_tests
//ASD
//ASD//-----------------------------------------------------------------------------------------------//
//ASD  KOKKOS_FUNCTION static void T_th_conversion_test(const Scalar& T_in, const Scalar& pres_in, int& errors){
//ASD
//ASD    // Allow usage of universal functions
//ASD    using physics = scream::physics::Functions<Scalar, Device>;
//ASD    // Gather the test tolerance
//ASD    static constexpr Scalar eps = C::macheps;
//ASD    Real tol = 1000*eps;
//ASD
//ASD    //========================================================
//ASD    // Test conversion of temperature to potential temperature, and vice versa
//ASD    //========================================================
//ASD    // This function tests both the T_to_th and th_to_T universal conversion functions.
//ASD    //
//ASD    // Inputs:
//ASD    //   T_in:    An example temperature, K.
//ASD    //   pres_in: A pressure value to use for conversion, Pa.
//ASD    // Outputs:
//ASD    //   errors:  A tally of any errors that this test detects.
//ASD    //========================================================
//ASD
//ASD    const Spack T(T_in);
//ASD    const Spack pres(pres_in);
//ASD    // Retrieve exners function for this pressure level.  Note, get_exner is tested in a separate test.
//ASD    const Spack exner = physics::get_exner(pres,Smask(true));
//ASD
//ASD    // Test conversion from temperature (T) to potential temperature (th)
//ASD    // Use T_in to convert from T to th
//ASD    // Note: T to th conversion is just a simple formula of th = T/exner
//ASD    const Spack th_atm = physics::T_to_th(T,exner,Smask(true));
//ASD    Real expected_th = T_in/exner[0];
//ASD    if (std::abs(th_atm[0]-expected_th)>tol) {
//ASD      printf("T to th test: abs(th_atm-expected_th)=%e is larger than the tol=%e\n",std::abs(th_atm[0]-expected_th),tol);
//ASD      errors++;
//ASD    }
//ASD
//ASD    // Test conversion from potential temperature (th) to temperature (T)
//ASD    // Use T_in to convert from th to T, note, we use T_in again because we are just testing the formula and T_in should be sufficient for that.
//ASD    // Note: th to T conversion is just a simple formula of T = th*exner
//ASD    const Spack T_atm = physics::th_to_T(T,exner,Smask(true));
//ASD    Real expected_T = T_in*exner[0];
//ASD    if (std::abs(T_atm[0]-expected_T)>tol) {
//ASD      printf("T to th test: abs(th_atm-expected_th)=%e is larger than the tol=%e\n",std::abs(T_atm[0]-expected_T),tol);
//ASD      errors++;
//ASD    }
//ASD
//ASD    // Test that T_to_th and th_to_T are inverses of each other.
//ASD    // Test T to th as inverses of each other
//ASD    // Converting T to th and then back to T should return the same value again.
//ASD    const Spack th_temp = physics::T_to_th(T,exner,Smask(true));
//ASD    const Spack T_new   = physics::th_to_T(th_temp,exner,Smask(true));
//ASD    if (std::abs(T_new[0]-T_in)>tol) {
//ASD      printf("T to th test: abs[ T_new (%.3e) - T_in (%.3e) ]=%e is larger than the tol=%e\n",T_new[0],T_in,std::abs(T_new[0]-T_in),tol);
//ASD      errors++;
//ASD    }
//ASD    // and vice versa
//ASD    const Spack T_temp = physics::th_to_T(T,exner,Smask(true));
//ASD    const Spack th_new = physics::T_to_th(T_temp,exner,Smask(true));
//ASD    if (std::abs(th_new[0]-T_in)>tol) {
//ASD      printf("T to th test: abs[ th_new (%.3e) - T_in (%.3e) ]=%e is larger than the tol=%e\n",th_new[0],T_in,std::abs(th_new[0]-T_in),tol);
//ASD      errors++;
//ASD    }
//ASD  } // T_th_conversion_test
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
//ASD  KOKKOS_FUNCTION static void dse_tests(const Scalar& temp, const Scalar& height, const Scalar surface_height, int& errors){
//ASD
//ASD    // Allow usage of universal functions
//ASD    using physics = scream::physics::Functions<Scalar, Device>;
//ASD    // Gather the test tolerance
//ASD    static constexpr Scalar eps = C::macheps;
//ASD    Real tol = 1000*eps;
//ASD
//ASD    //========================================================
//ASD    // Test calculation of dry static energy using get_dse
//ASD    //========================================================
//ASD    // This function tests the function "get_dse".
//ASD    //
//ASD    // Inputs:
//ASD    //   temp is the atmospheric temperature. Units in K.
//ASD    //   height is the geopotential height above surface at midpoints. Units in m.
//ASD    //   surface_height is the surface geopotential height. Units in m.
//ASD    // Outputs:
//ASD    //   errors:  A tally of any errors that this test detects.
//ASD    //========================================================
//ASD
//ASD    static constexpr Scalar cp = C::CP;
//ASD    static constexpr Scalar ggr = C::gravit;
//ASD    const Spack T_mid(temp);
//ASD    const Spack z_mid(height);
//ASD    Real expected_dse = cp*temp+ggr*height+surface_height;
//ASD    const Spack dse = physics::get_dse(T_mid, z_mid, surface_height, Smask(true));
//ASD    if (std::abs(dse[0]-expected_dse)>tol) {
//ASD      printf("get_dse test: abs(dse-expected_dse)=%e is larger than the tol=%e\n",std::abs(dse[0]-expected_dse),tol);
//ASD      errors++;
//ASD    }
//ASD
//ASD  } // dse_test
//-----------------------------------------------------------------------------------------------//
//ASD  KOKKOS_FUNCTION static void virtual_temperature_test(const Scalar& T_mid_in, const Scalar& qv_in, int& errors){
//ASD
//ASD    // Allow usage of universal functions
//ASD    using physics = scream::physics::Functions<Scalar, Device>;
//ASD    // Gather the machine epsilon for the error tolerance
//ASD    static constexpr Scalar eps = C::macheps;
//ASD    static constexpr Scalar ep_2 = C::ep_2;
//ASD    Real tol = 2000*eps;
//ASD
//ASD    //========================================================
//ASD    // Test calculation of virtual temperature using get_virtual_temperature
//ASD    //========================================================
//ASD    // This function tests the function "get_virtual_temperature".
//ASD    //
//ASD    // Inputs:
//ASD    //   T_mid_in is the atmospheric temperature. Units in K.
//ASD    //   qv_in is the atmospheric water vapor mass mixing ratio.  Units in kg/kg
//ASD    // Outputs:
//ASD    //   errors:  A tally of any errors that this test detects.
//ASD    //========================================================
//ASD    const Spack T_mid(T_mid_in);
//ASD    const Spack qv(qv_in);
//ASD    // Given T and qv, determine T_virt
//ASD    const Spack T_virt = physics::get_virtual_temperature(T_mid,qv,Smask(true));
//ASD    // Determine the expected virtual temperature
//ASD    Real expected_T_virt = T_mid_in*(qv_in+ep_2)/(ep_2*(1.0+qv_in));
//ASD
//ASD    if (std::abs(T_virt[0]-expected_T_virt)>tol) {
//ASD      printf("get_virtual_temperature test: abs(T_virt-expected_T_virt)=%e is larger than the tol=%e\n",std::abs(T_virt[0]-expected_T_virt),tol);
//ASD      errors++;
//ASD    }
//ASD  
//ASD  } // virtual_temperature_test
//-----------------------------------------------------------------------------------------------//
  static void run()
  {
    using physics = scream::physics::Functions<Scalar, Device>;
    using physicscommon = scream::physics::PhysicsFunctions<Device>;

    int num_levs = 100; // Number of levels to use for tests.

    static constexpr Scalar p0     = C::P0;
    static constexpr Scalar Rd     = C::RD;
    static constexpr Scalar inv_cp = C::INV_CP;
    static constexpr Scalar tmelt  = C::Tmelt;
    static constexpr Scalar ggr    = C::gravit;
    static constexpr Scalar eps = C::macheps;
    Real tol = 1000*eps;

    // Compute random values for dse test
    using view_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Real>;
    view_1d temperature("temperature",num_levs),
            height("height",num_levs),
            surface_height("surface_height",1),
            qv("qv",num_levs),
            pressure("pressure",num_levs),
            psuedo_density("psuedo_density",num_levs),
            dz_for_testing("dz_for_testing",num_levs);
    view_1d exner("exner",num_levs),
            theta("theta",num_levs),
            T_mid("T_mid",num_levs),
            T_virtual("T_virtual",num_levs),
            dse("dse",num_levs),
            dz("dz",num_levs),
            z_int("z_int",num_levs+1);

    std::random_device rdev;
    using rngAlg = std::mt19937_64;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rdev() : catchRngSeed;
    rngAlg engine(seed);
    using RPDF = std::uniform_real_distribution<Real>;
    RPDF pdf_qv(1e-3,1e3),
         pdf_dp(1.0,100.0),
         pdf_pres(0.0,p0),
         pdf_temp(200.0,400.0),
         pdf_height(0.0,1e5),
         pdf_surface(100.0,400.0);


    ekat::genRandArray(temperature,engine,pdf_temp);
    ekat::genRandArray(height,engine,pdf_height);
    ekat::genRandArray(surface_height,engine,pdf_surface);
    ekat::genRandArray(qv,engine,pdf_qv);
    ekat::genRandArray(pressure,engine,pdf_pres);
    ekat::genRandArray(psuedo_density,engine,pdf_dp);

    int nerr = 0;
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("test_universal_physics", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors += 0;

      
      // Create dummy level data for testing:
      Real pres_top = 200.;
      Real dp       = (p0-pres_top)/(num_levs-1);

      // Run tests
      // Exner property tests:
      // get_exner(p0) should return 1.0
      // get_exner(0.0) should return 0.0
      // get_exner(2*p0) should return 2**(Rd/cp)
      // get_exner(pressure) should work.
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_exner(p0)==1.0,"Error in get_exner(p0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_exner(0.0)==0.0,"Error in get_exner(0.0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_exner(2*p0)==pow(2.0,Rd*inv_cp),"Error in get_exner(2*p0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_exner(4.0)/physicscommon::get_exner(2.0)==pow(2.0,Rd*inv_cp),"Error in get_exner(4)/get_exner(2) property test");
      physicscommon::get_exner(team,pressure,exner);
      // Potential temperature property tests
      // theta=T when p=p0
      // theta(T=0) = 0
      // T(theta=0) = 0
      // T(theta(T0)) = T0
      // theta(T(theta0)) = theta0
      // get_theta_from_T(T,pressure) should work
      // get_T_from_theta(theta,pressure) should work
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_theta_from_T(100.0,p0)==100.0,"Error in get_theta_from_T(100,p0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_theta_from_T(0.0,1.0)==0.0,"Error in get_theta_from_T(0.0,1.0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_T_from_theta(0.0,1.0)==0.0,"Error in get_T_from_theta(0.0,1.0) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_T_from_theta(physicscommon::get_theta_from_T(100.0,1.0),1.0)==100.0,"Error in T->theta->T test"); 
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_theta_from_T(physicscommon::get_T_from_theta(100.0,1.0),1.0)==100.0,"Error in theta->T->theta test");
      physicscommon::get_theta_from_T(team,temperature,pressure,theta);
      physicscommon::get_T_from_theta(team,theta,pressure,T_mid);
      // Virtual temperature property tests
      // T_virt(T=0) = 0.0
      // T_virt(T=T0,qv=0) = T0
      // get_virtual_temperature(temperature,qv) should work
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_virtual_temperature(0.0,1e-6)==0.0,"Error in get_virtual_temperature(0,qv) property test");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_virtual_temperature(100.0,0.0)==100.0,"Error in get_virtual_temperature(T,0) property test");
      physicscommon::get_virtual_temperature(team,temperature,qv,T_virtual);
      // DSE property tests
      // dse(T=0.0, z=0.0) = surf_geopotential
      // dse(T=1/cp, z=1/gravity) = surf_potential+2
      // get_dse should work
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dse(0.0,0.0,10.0)==10.0,"");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dse(inv_cp,1.0/ggr,0.0)==2.0,"");
      physicscommon::get_dse(team,temperature,height,surface_height(0),dse);
      // DZ tests
      // get_dz should work
      // dz(psuedo_density=0) = 0
      // dz(T=0) = 0
      // dz(psuedo_density=p0,p_mid=p0,T=1.0,qv=0) = Rd/ggr
      // dz(psuedo_density=ggr,p_mid=Rd,T=T0,qv=0) = T0
      physicscommon::get_dz(team,psuedo_density,pressure,temperature,qv,dz);
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dz(0.0,p0,100.0,1e-5)==0.0,"");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dz(100.0,p0,0.0,1e-5)==0.0,"");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dz(p0,p0,1.0,0.0)==Rd/ggr,"");
      EKAT_KERNEL_REQUIRE_MSG(physicscommon::get_dz(ggr,Rd,100.0,0.0)==100.0,"");
      // z_int property tests
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs), [&] (const int k)
      {
        dz_for_testing(k) = k+1;
      });
      physicscommon::get_z_int(team,dz_for_testing,z_int);

//ASD      for (int k=0;k<num_levs;++k)
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs), [&] (const int k)
      {
        // Exner test
//TEMP        EKAT_KERNEL_REQUIRE_MSG(exner(k)==physicscommon::get_exner(pressure(k)),"Error is column-wise get_exner at level = " 
//TEMP             + std::to_string(k) + ": exner(k)=" + std::to_string(exner(k)) 
//TEMP             + " vs. get_exner(P=" + std::to_string(pressure(k)) + ")=" + std::to_string(physicscommon::get_exner(pressure(k))));
//ASD        //   - Pressure at level k is just k*dp, so each pressure level is dp in thickness.
//ASD        Real pres = p0 + k*dp;
//ASD        Real p_val  = pow( pres/p0, Rd*inv_cp);
//ASD        exner_tests(pres,p_val,errors);
        // T and TH conversion TEST
        EKAT_KERNEL_REQUIRE_MSG(theta(k)==physicscommon::get_theta_from_T(temperature(k),pressure(k)),"");
        EKAT_KERNEL_REQUIRE_MSG(T_mid(k)==physicscommon::get_T_from_theta(theta(k),pressure(k)),"");
//TEMP        EKAT_KERNEL_REQUIRE_MSG(std::abs(T_mid(k)-temperature(k))<tol,"Error in checking temperature->potential->temperature at lev = " 
//TEMP                      + std::to_string(k) + ", (T_1,T_0,,T_1-T_0P) = (" + std::to_string(T_mid(k)) + ", " + std::to_string(temperature(k)) 
//TEMP                      + ", " + std::to_string(std::abs(T_mid(k)-temperature(k))) + ", " + std::to_string(pressure(k)) + ")");
//ASD        Real T_atm = tmelt + k;
//ASD        T_th_conversion_test(T_atm,pres,errors);
        // virtual temperature test
//ASD        virtual_temperature_test(temperature(k),qv(k),errors);
        EKAT_KERNEL_REQUIRE_MSG(T_virtual(k)==physicscommon::get_virtual_temperature(temperature(k),qv(k)),"");
        // --  T_virtual must always be larger than T because moist air is less dense than dry air at the same temperature and pressure
//TEMP        EKAT_KERNEL_REQUIRE_MSG(T_virtual(k)>temperature(k),"Error in virtual temperature tests at lev = "
//TEMP                      + std::to_string(k) + ", T_virtual=" + std::to_string(T_virtual(k)) + " !> T_mid=" + std::to_string(temperature(k))); 
        // DZ Test
        //   - Determine the z-layer bottom and top as being the (k+1) thick.
        //     Allow for varying interface heights by setting zi_bot equal to the sum(0,...,k).
        Real zi_bot = (pow(k,2)+k)/2.0;
        Real zi_top = zi_bot + (k+1);
        dz_tests(zi_top,zi_bot,errors);
        // DSE Test
//ASD        dse_tests(temperature(k),height(k),surface_height(k),errors);
      });  // Kokkos parallel_for k=1,num_levs

    }, nerr); // Kokkos parallel_reduce "test_universal_physics"

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
