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

    using view_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Real>;
    // Compute random values for tests
    view_1d temperature("temperature",num_levs),
            height("height",num_levs),
            surface_height("surface_height",1),
            qv("qv",num_levs),
            pressure("pressure",num_levs),
            psuedo_density("psuedo_density",num_levs),
            dz_for_testing("dz_for_testing",num_levs);
    // Allocate memory for test outputs
    view_1d exner("exner",num_levs),
            theta("theta",num_levs),
            T_mid("T_mid",num_levs),
            T_virtual("T_virtual",num_levs),
            dse("dse",num_levs),
            dz("dz",num_levs),
            z_int("z_int",num_levs+1);

    // Construct random input data
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

    // Construct a simple set of `dz` values for testing the z_int function
    auto dz_for_testing_host = Kokkos::create_mirror(dz_for_testing);
    for (int k = 0;k<num_levs;k++)
    {
      dz_for_testing_host(k) = num_levs-k;
    }
    Kokkos::deep_copy(dz_for_testing,dz_for_testing_host);

    // Run tests
    int nerr = 0;
    TeamPolicy policy(ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1));
    Kokkos::parallel_reduce("test_universal_physics", policy, KOKKOS_LAMBDA(const MemberType& team, int& errors) {

      errors += 0;

      // Run tests
      // Exner property tests:
      // get_exner(p0) should return 1.0
      // get_exner(0.0) should return 0.0
      // get_exner(2*p0) should return 2**(Rd/cp)
      // get_exner(pressure) should work.
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(p0)==1.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(0.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(2*p0)==pow(2.0,Rd*inv_cp));
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(4.0)/physicscommon::get_exner(2.0)==pow(2.0,Rd*inv_cp));
      physicscommon::get_exner(team,pressure,exner);
      // Potential temperature property tests
      // theta=T when p=p0
      // theta(T=0) = 0
      // T(theta=0) = 0
      // T(theta(T0)) = T0
      // theta(T(theta0)) = theta0
      // get_theta_from_T(T,pressure) should work
      // get_T_from_theta(theta,pressure) should work
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(100.0,p0)==100.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(0.0,1.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_T_from_theta(0.0,1.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_T_from_theta(physicscommon::get_theta_from_T(100.0,1.0),1.0)==100.0); 
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(physicscommon::get_T_from_theta(100.0,1.0),1.0)==100.0);
      physicscommon::get_theta_from_T(team,temperature,pressure,theta);
      physicscommon::get_T_from_theta(team,theta,pressure,T_mid);
      // Virtual temperature property tests
      // T_virt(T=0) = 0.0
      // T_virt(T=T0,qv=0) = T0
      // get_virtual_temperature(temperature,qv) should work
      EKAT_KERNEL_REQUIRE(physicscommon::get_virtual_temperature(0.0,1e-6)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_virtual_temperature(100.0,0.0)==100.0);
      physicscommon::get_virtual_temperature(team,temperature,qv,T_virtual);
      // DSE property tests
      // dse(T=0.0, z=0.0) = surf_geopotential
      // dse(T=1/cp, z=1/gravity) = surf_potential+2
      // get_dse should work
      EKAT_KERNEL_REQUIRE(physicscommon::get_dse(0.0,0.0,10.0)==10.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dse(inv_cp,1.0/ggr,0.0)==2.0);
      physicscommon::get_dse(team,temperature,height,surface_height(0),dse);
      // DZ tests
      // get_dz should work
      // dz(psuedo_density=0) = 0
      // dz(T=0) = 0
      // dz(psuedo_density=p0,p_mid=p0,T=1.0,qv=0) = Rd/ggr
      // dz(psuedo_density=ggr,p_mid=Rd,T=T0,qv=0) = T0
      physicscommon::get_dz(team,psuedo_density,pressure,temperature,qv,dz);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(0.0,p0,100.0,1e-5)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(100.0,p0,0.0,1e-5)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(p0,p0,1.0,0.0)==Rd/ggr);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(ggr,Rd,100.0,0.0)==100.0);
      // z_int property tests
      physicscommon::get_z_int(team,dz_for_testing,z_int);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs), [&] (const int k)
      {
        // T and TH conversion TEST
        EKAT_KERNEL_REQUIRE(theta(k)==physicscommon::get_theta_from_T(temperature(k),pressure(k)));
        EKAT_KERNEL_REQUIRE(T_mid(k)==physicscommon::get_T_from_theta(theta(k),pressure(k)));
        // virtual temperature test
        EKAT_KERNEL_REQUIRE(T_virtual(k)==physicscommon::get_virtual_temperature(temperature(k),qv(k)));
      });  // Kokkos parallel_for k=1,num_levs

    }, nerr); // Kokkos parallel_reduce "test_universal_physics"

    Kokkos::fence();
    auto z_int_host = Kokkos::create_mirror(z_int);
    Kokkos::deep_copy(z_int_host,z_int);
    //Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs+1), [&] (const int k)
    for(int k=0;k<num_levs;k++)
    {
      const auto k_bwd = num_levs-k;
      std::string z_int_err_msg = "Error at level k,k_bwd=" + std::to_string(k) + "," + std::to_string(k_bwd) + "; z_int(k_bwd) = " + std::to_string(z_int_host(k_bwd));
      EKAT_KERNEL_REQUIRE_MSG(z_int_host(k_bwd)==k*(k+1)/2,z_int_err_msg.c_str());
    }

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
