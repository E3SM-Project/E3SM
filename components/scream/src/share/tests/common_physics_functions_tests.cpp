#include "catch2/catch.hpp"

#include "share/util/scream_common_physics_functions.hpp"
#include "share/util/scream_common_physics_impl.hpp"
#include "physics/share/tests/physics_unit_tests_common.hpp"

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
  template<typename ScalarT, typename DeviceT>
  static void run_impl()
  {
    using physicscommon  = scream::PhysicsFunctions<Device>;
    using view_1d        = ekat::KokkosTypes<DefaultDevice>::view_1d<ScalarT>;
    using scalar_view_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Scalar>;

    static constexpr Scalar p0     = C::P0;
    static constexpr Scalar Rd     = C::RD;
    static constexpr Scalar inv_cp = C::INV_CP;
    static constexpr Scalar tmelt  = C::Tmelt;
    static constexpr Scalar ggr    = C::gravit;
    static constexpr Scalar test_tol = C::macheps*1e3;

    const int num_levs = 100; // Number of levels to use for tests.
    const Int num_pack = ekat::npack<ScalarT>(num_levs);
    const Int num_pack_int = ekat::npack<ScalarT>(num_levs+1);

    // Compute random values for tests
    view_1d        temperature("temperature",num_pack),
                   height("height",num_pack),
                   qv("qv",num_pack),
                   pressure("pressure",num_pack),
                   pseudo_density("pseudo_density",num_pack),
                   dz_for_testing("dz_for_testing",num_pack);
    scalar_view_1d surface_height("surface_height",1);
    view_1d        exner("exner",num_pack),
                   theta("theta",num_pack),
                   T_mid_from_pot("T_mid_from_pot",num_pack),
                   T_virtual("T_virtual",num_pack),
                   T_mid_from_virt("T_mid_from_virt",num_pack),
                   dse("dse",num_pack),
                   dz("dz",num_pack),
                   z_int("z_int",num_pack_int);
    // Construct random input data
    std::random_device rdev;
    using rngAlg = std::mt19937_64;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rdev() : catchRngSeed;
    // Print seed to screen to trace tests that fail.
    std::cout << "  - seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);
    using RPDF = std::uniform_real_distribution<Scalar>;
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
    ekat::genRandArray(pseudo_density,engine,pdf_dp);

    // Construct a simple set of `dz` values for testing the z_int function
    auto dz_for_testing_host = Kokkos::create_mirror_view(ekat::scalarize(dz_for_testing));
    for (int k = 0;k<num_levs;k++)
    {
      dz_for_testing_host[k] = num_levs-k;
    }
    Kokkos::deep_copy(dz_for_testing,dz_for_testing_host);

    // Run tests using single values
    Scalar t_result, T0, z0, ztest, ptest, dp0, qv0;
    // Exner property tests:
    // exner_function(p0) should return 1.0
    // exner_function(0.0) should return 0.0
    // exner_function(2*p0) should return 2**(Rd/cp)
    ptest = p0;
    REQUIRE(physicscommon::exner_function(ptest)==1.0);
    ptest = 0.0;
    REQUIRE(physicscommon::exner_function(ptest)==0.0);
    ptest = 4.0; t_result = pow(2.0,Rd*inv_cp);
    REQUIRE(std::abs(physicscommon::exner_function(ptest)/physicscommon::exner_function(ptest/2)-t_result)<test_tol);
    // Potential temperature property tests
    // theta=T when p=p0
    // theta(T=0) = 0
    // T(theta=0) = 0
    // T(theta(T0)) = T0
    // theta(T(theta0)) = theta0
    T0 = 100.0;
    REQUIRE(physicscommon::calculate_theta_from_T(T0,p0)==T0);
    REQUIRE(physicscommon::calculate_theta_from_T(0.0,1.0)==0.0);
    REQUIRE(physicscommon::calculate_T_from_theta(0.0,1.0)==0.0);
    REQUIRE(physicscommon::calculate_T_from_theta(physicscommon::calculate_theta_from_T(100.0,1.0),1.0)==100.0); 
    REQUIRE(physicscommon::calculate_theta_from_T(physicscommon::calculate_T_from_theta(100.0,1.0),1.0)==100.0);
    // Virtual temperature property tests
    // T_virt(T=0) = 0.0
    // T_virt(T=T0,qv=0) = T0
    // T(T_virt=0) = 0.0
    // T(T_virt=T0,qv=0) = T0
    // T_virt(T=T0) = T0
    // T(T_virt=T0) = T0
    REQUIRE(physicscommon::calculate_virtual_temperature(0.0,1e-6)==0.0);
    REQUIRE(physicscommon::calculate_virtual_temperature(100.0,0.0)==100.0);
    REQUIRE(physicscommon::calculate_temperature_from_virtual_temperature(0.0,1e-6)==0.0);
    REQUIRE(physicscommon::calculate_temperature_from_virtual_temperature(100.0,0.0)==100.0);
    REQUIRE(physicscommon::calculate_virtual_temperature(physicscommon::calculate_temperature_from_virtual_temperature(100.0,1.0),1.0)==100.0); 
    REQUIRE(physicscommon::calculate_temperature_from_virtual_temperature(physicscommon::calculate_virtual_temperature(100.0,1.0),1.0)==100.0);
    // DSE property tests
    // calculate_dse(T=0.0, z=0.0) = surf_geopotential
    // calculate_dse(T=1/cp, z=1/gravity) = surf_potential+2
    T0=0.0; ztest=0.0; z0=10.0;
    REQUIRE(physicscommon::calculate_dse(T0,ztest,z0)==10.0);
    T0=inv_cp; ztest=1.0/ggr; z0=0.0;
    REQUIRE(physicscommon::calculate_dse(T0,ztest,z0)==z0+2.0);
    // DZ tests
    // calculate_dz(pseudo_density=0) = 0
    // calculate_dz(T=0) = 0
    // calculate_dz(pseudo_density=p0,p_mid=p0,T=1.0,qv=0) = Rd/ggr
    // calculate_dz(pseudo_density=ggr,p_mid=Rd,T=T0,qv=0) = T0
    dp0=0.0; ptest=p0; T0=100.0; qv0=1e-5;
    REQUIRE(physicscommon::calculate_dz(dp0,ptest,T0,qv0)==0.0);
    dp0=100.0; ptest=p0; T0=0.0; qv0=1e-5;
    REQUIRE(physicscommon::calculate_dz(dp0,ptest,T0,qv0)==0.0);
    dp0=p0; ptest=p0; T0=1.0; qv0=0.0;
    REQUIRE(physicscommon::calculate_dz(dp0,ptest,T0,qv0)==Rd/ggr);
    dp0=ggr; ptest=Rd; T0=100.0; qv0=0.0;
    REQUIRE(physicscommon::calculate_dz(dp0,ptest,T0,qv0)==T0);
    // Run tests on full views
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_for("test_universal_physics", policy, KOKKOS_LAMBDA(const MemberType& team) {
      // Constants used in testing (necessary to ensure SP build works)
      // Exner property tests:
      physicscommon::exner_function(team,pressure,exner);
      // Potential temperature property tests
      // calculate_theta_from_T(T,pressure) should work
      physicscommon::calculate_theta_from_T(team,temperature,pressure,theta);
      physicscommon::calculate_T_from_theta(team,theta,pressure,T_mid_from_pot);
      // Virtual temperature property tests
      // calculate_virtual_temperature(temperature,qv) should work 
      // calculate_temperature_from_virtual_temperature should work 
      physicscommon::calculate_virtual_temperature(team,temperature,qv,T_virtual);
      physicscommon::calculate_temperature_from_virtual_temperature(team,T_virtual,qv,T_mid_from_virt);
      // DSE property tests
      // calculate_dse should work
      physicscommon::calculate_dse(team,temperature,height,surface_height(0),dse);
      // DZ tests
      // calculate_dz should work
      physicscommon::calculate_dz(team,pseudo_density,pressure,temperature,qv,dz);
      // z_int property tests
      // calculate_z_int should work
      physicscommon::calculate_z_int(team,num_levs,dz_for_testing,z_int);
    }); // Kokkos parallel_for "test_universal_physics"
    Kokkos::fence();

    // Check the properties of the full view outputi
    // - map to host view
    auto exner_host           = Kokkos::create_mirror_view(exner);
    auto theta_host           = Kokkos::create_mirror_view(theta);
    auto T_mid_from_pot_host  = Kokkos::create_mirror_view(T_mid_from_pot);
    auto T_virtual_host       = Kokkos::create_mirror_view(T_virtual);
    auto T_mid_from_virt_host = Kokkos::create_mirror_view(T_mid_from_virt);
    auto dse_host             = Kokkos::create_mirror_view(dse);
    auto z_int_host           = Kokkos::create_mirror_view(z_int);
    auto dz_host              = Kokkos::create_mirror_view(dz);
    auto temperature_host     = Kokkos::create_mirror_view(temperature);
    auto pressure_host        = Kokkos::create_mirror_view(pressure);
    auto qv_host              = Kokkos::create_mirror_view(qv);
    Kokkos::deep_copy(exner_host           , exner);
    Kokkos::deep_copy(theta_host           , theta);
    Kokkos::deep_copy(T_mid_from_pot_host  , T_mid_from_pot);
    Kokkos::deep_copy(T_virtual_host       , T_virtual);
    Kokkos::deep_copy(T_mid_from_virt_host , T_mid_from_virt);
    Kokkos::deep_copy(dse_host             , dse);
    Kokkos::deep_copy(z_int_host           , z_int);
    Kokkos::deep_copy(dz_host              , dz);
    Kokkos::deep_copy(temperature_host     , temperature);
    Kokkos::deep_copy(pressure_host        , pressure);
    Kokkos::deep_copy(qv_host              , qv);
    // Now check for all levels
    for(int k=0;k<num_pack;k++)
    {
      const auto range_pack = ekat::range<ScalarT>(k*ScalarT::n);
      const auto range_mask = range_pack < num_levs;
      // Make sure all columnwise results don't contain any obvious errors:
      // exner
      REQUIRE( !((isnan(exner_host(k)) && range_mask).any()) );
      REQUIRE( !(((exner_host(k)<0) && range_mask).any()) );
      // potential temperature
      REQUIRE( !((isnan(theta_host(k)) && range_mask).any()) );
      REQUIRE( !(((theta_host(k)<0) && range_mask).any()) );
      {
        auto theta_from_t = physicscommon::calculate_theta_from_T(temperature_host(k),pressure_host(k));
        auto test_theta_from_t = !(theta_host(k)==theta_from_t) && range_mask;
        REQUIRE( !(test_theta_from_t.any()) );
      }
      // temperature from potential temperature
      REQUIRE( !((isnan(T_mid_from_pot_host(k)) && range_mask).any()) );
      REQUIRE( !(((T_mid_from_pot_host(k)<0) && range_mask).any()) );
      {
        auto t_from_theta = physicscommon::calculate_T_from_theta(theta_host(k),pressure_host(k));
        auto test_t_from_theta = !(T_mid_from_pot_host(k)==t_from_theta) && range_mask;
        REQUIRE( !(test_t_from_theta.any()) );
      }
      // virtual temperature
      REQUIRE( !((isnan(T_virtual_host(k)) && range_mask).any()) );
      REQUIRE( !(((T_virtual_host(k)<0) && range_mask).any()) );
      {
        auto virt_temp = physicscommon::calculate_virtual_temperature(temperature_host(k),qv_host(k));
        auto test_virt_temp = !( T_virtual_host(k)==virt_temp ) && range_mask;
        REQUIRE( !(test_virt_temp.any()) );
      }
      // temperature from virtual temperature
      REQUIRE( !((isnan(T_mid_from_virt_host(k)) && range_mask).any()) );
      REQUIRE( !(((T_mid_from_virt_host(k)<0) && range_mask).any()) );
      // DSE
      REQUIRE( !((isnan(dse_host(k)) && range_mask).any()) );
      REQUIRE( !(((dse_host(k)<0) && range_mask).any()) );
      // dz
      REQUIRE( !((isnan(dz_host(k)) && range_mask).any()) );
      REQUIRE( !(((dz_host(k)<=0) && range_mask).any()) );
      // z_int
//TO-FIX      for (int i=0;i<Spack::n;i++) {
//TO-FIX        Int lev = k*Spack::n + i;  // Determine the level at this pack/vec location.
//TO-FIX        const auto lev_bwd = num_levs-lev;
//TO-FIX        Int ipack = lev_bwd / Spack::n;
//TO-FIX        Int ivec  = lev_bwd % Spack::n;
//TO-FIX        REQUIRE( z_int_host(ipack)[ivec]==lev*(lev+1)/2);
//TO-FIX      }
      REQUIRE( !((isnan(z_int_host(k)) && range_mask).any()) );
      REQUIRE( !(((z_int_host(k)<0) && range_mask).any()) );
    }
  } // run_impl()

//-----------------------------------------------------------------------------------------------//
  static void run()
  {
    // Run tests for different pack size types
    // Note, we print to screen the test type so it can be easily matched with the seed used for
    // randomization within the test.
    printf("Testing Scalar: packsize=1\n");
    run_impl<ekat::Pack<Scalar,1>,Device>();
    printf("Testing SMALL_PACK: packsize=%d\n",SCREAM_SMALL_PACK_SIZE);
    run_impl<ekat::Pack<Scalar,SCREAM_SMALL_PACK_SIZE>,Device>();
    printf("Testing PACK: packsize=%d\n",SCREAM_PACK_SIZE);
    run_impl<ekat::Pack<Scalar,SCREAM_PACK_SIZE>,Device>();
  } // run
}; // end of TestUniversal struct

} // namespace unit_test
} // namespace physics
} // namespace scream

namespace{

TEST_CASE("common_physics_functions_test", "[common_physics_functions_test]"){
  scream::physics::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUniversal::run();

 } // TEST_CASE

} // namespace
