#include "catch2/catch.hpp"

#include "share/util/scream_common_physics_functions.hpp"
#include "physics/share/physics_universal_impl.hpp"
#include "share/util/scream_common_physics_impl.hpp"
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
    using physicscommon = scream::PhysicsFunctions<Device>;

    using Spack = ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>;

    int num_levs = 100; // Number of levels to use for tests.

    static constexpr Scalar p0     = C::P0;
    static constexpr Scalar Rd     = C::RD;
    static constexpr Scalar inv_cp = C::INV_CP;
    static constexpr Scalar tmelt  = C::Tmelt;
    static constexpr Scalar ggr    = C::gravit;

    using view_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Real>;
    using sview_1d = ekat::KokkosTypes<DefaultDevice>::view_1d<Spack>;
    const Int num_pack = ekat::npack<Spack>(num_levs);
    const Int num_pack_int = ekat::npack<Spack>(num_levs+1);
    // Compute random values for tests
    sview_1d temperature_packed("temperature",num_pack),
             height_packed("height",num_pack),
             qv_packed("qv",num_pack),
             pressure_packed("pressure",num_pack),
             psuedo_density_packed("psuedo_density",num_pack),
             dz_for_testing_packed("dz_for_testing",num_pack);
    view_1d surface_height("surface_height",1);
    // Allocate memory for test outputs
    view_1d exner("exner",num_levs),
            theta("theta",num_levs),
            T_mid_from_pot("T_mid_from_pot",num_levs),
            T_virtual("T_virtual",num_levs),
            T_mid_from_virt("T_mid_from_virt",num_levs),
            dse("dse",num_levs),
            dz("dz",num_levs),
            z_int("z_int",num_levs+1);
    sview_1d exner_packed("exner",num_pack),
             theta_packed("theta",num_pack),
             T_mid_from_pot_packed("T_mid_from_pot",num_pack),
             T_virtual_packed("T_virtual",num_pack),
             T_mid_from_virt_packed("T_mid_from_virt",num_pack),
             dse_packed("dse",num_pack),
             dz_packed("dz",num_pack),
             z_int_packed("z_int",num_pack_int);

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

    ekat::genRandArray(temperature_packed,engine,pdf_temp);
    ekat::genRandArray(height_packed,engine,pdf_height);
    ekat::genRandArray(surface_height,engine,pdf_surface);
    ekat::genRandArray(qv_packed,engine,pdf_qv);
    ekat::genRandArray(pressure_packed,engine,pdf_pres);
    ekat::genRandArray(psuedo_density_packed,engine,pdf_dp);

    // Construct a simple set of `dz` values for testing the z_int function
    auto dz_for_testing_host = Kokkos::create_mirror(dz_for_testing_packed);
    for (int k = 0;k<num_levs;k++)
    {
      int ipack = k / Spack::n;
      int ivec  = k % Spack::n;
      dz_for_testing_host(ipack)[ivec] = num_levs-k;
    }
    Kokkos::deep_copy(dz_for_testing_packed,dz_for_testing_host);

    view_1d temperature(reinterpret_cast<Real*>(temperature_packed.data()),num_levs),
            height(reinterpret_cast<Real*>(height_packed.data()),num_levs),
            qv(reinterpret_cast<Real*>(qv_packed.data()),num_levs),
            pressure(reinterpret_cast<Real*>(pressure_packed.data()),num_levs),
            psuedo_density(reinterpret_cast<Real*>(psuedo_density_packed.data()),num_levs),
            dz_for_testing(reinterpret_cast<Real*>(dz_for_testing_packed.data()),num_levs);

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
      // get_exner(pressure) should work for Real and for Spack.
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(p0)==1.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(0.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(2*p0)==pow(2.0,Rd*inv_cp));
      EKAT_KERNEL_REQUIRE(physicscommon::get_exner(4.0)/physicscommon::get_exner(2.0)==pow(2.0,Rd*inv_cp));
      physicscommon::get_exner(team,pressure,exner);
      physicscommon::get_exner(team,pressure_packed,exner_packed);
      // Potential temperature property tests
      // theta=T when p=p0
      // theta(T=0) = 0
      // T(theta=0) = 0
      // T(theta(T0)) = T0
      // theta(T(theta0)) = theta0
      // get_theta_from_T(T,pressure) should work
      // get_T_from_theta(theta,pressure) should work for Real and for Spack
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(100.0,p0)==100.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(0.0,1.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_T_from_theta(0.0,1.0)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_T_from_theta(physicscommon::get_theta_from_T(100.0,1.0),1.0)==100.0); 
      EKAT_KERNEL_REQUIRE(physicscommon::get_theta_from_T(physicscommon::get_T_from_theta(100.0,1.0),1.0)==100.0);
      physicscommon::get_theta_from_T(team,temperature,pressure,theta);
      physicscommon::get_T_from_theta(team,theta,pressure,T_mid_from_pot);
      physicscommon::get_theta_from_T(team,temperature_packed,pressure_packed,theta_packed);
      physicscommon::get_T_from_theta(team,theta_packed,pressure_packed,T_mid_from_pot_packed);
      // Virtual temperature property tests
      // T_virt(T=0) = 0.0
      // T_virt(T=T0,qv=0) = T0
      // T(T_virt=0) = 0.0
      // T(T_virt=T0,qv=0) = T0
      // T_virt(T=T0) = T0
      // T(T_virt=T0) = T0
      // get_virtual_temperature(temperature,qv) should work for Real and for Spack
      // get_temperature_from_virtual_temperature should work for Real and for Spack
      EKAT_KERNEL_REQUIRE(physicscommon::get_virtual_temperature(0.0,1e-6)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_virtual_temperature(100.0,0.0)==100.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_temperature_from_virtual_temperature(0.0,1e-6)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_temperature_from_virtual_temperature(100.0,0.0)==100.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_virtual_temperature(physicscommon::get_temperature_from_virtual_temperature(100.0,1.0),1.0)==100.0); 
      EKAT_KERNEL_REQUIRE(physicscommon::get_temperature_from_virtual_temperature(physicscommon::get_virtual_temperature(100.0,1.0),1.0)==100.0);
      physicscommon::get_virtual_temperature(team,temperature,qv,T_virtual);
      physicscommon::get_temperature_from_virtual_temperature(team,T_virtual,qv,T_mid_from_virt);
      physicscommon::get_virtual_temperature(team,temperature_packed,qv_packed,T_virtual_packed);
      physicscommon::get_temperature_from_virtual_temperature(team,T_virtual_packed,qv_packed,T_mid_from_virt_packed);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs), [&] (const int k)
      {
        // T and TH conversion TEST
        EKAT_KERNEL_REQUIRE(theta(k)==physicscommon::get_theta_from_T(temperature(k),pressure(k)));
        EKAT_KERNEL_REQUIRE(T_mid_from_pot(k)==physicscommon::get_T_from_theta(theta(k),pressure(k)));
        // virtual temperature test
        EKAT_KERNEL_REQUIRE(T_virtual(k)==physicscommon::get_virtual_temperature(temperature(k),qv(k)));
      });  // Kokkos parallel_for k=1,num_levs
      // DSE property tests
      // dse(T=0.0, z=0.0) = surf_geopotential
      // dse(T=1/cp, z=1/gravity) = surf_potential+2
      // get_dse should work for Real and for Spack
      EKAT_KERNEL_REQUIRE(physicscommon::get_dse(0.0,0.0,10.0)==10.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dse(inv_cp,1.0/ggr,0.0)==2.0);
      physicscommon::get_dse(team,temperature,height,surface_height(0),dse);
      physicscommon::get_dse(team,temperature_packed,height_packed,surface_height(0),dse_packed);
      // DZ tests
      // dz(psuedo_density=0) = 0
      // dz(T=0) = 0
      // dz(psuedo_density=p0,p_mid=p0,T=1.0,qv=0) = Rd/ggr
      // dz(psuedo_density=ggr,p_mid=Rd,T=T0,qv=0) = T0
      // get_dz should work for Real and for Spack
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(0.0,p0,100.0,1e-5)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(100.0,p0,0.0,1e-5)==0.0);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(p0,p0,1.0,0.0)==Rd/ggr);
      EKAT_KERNEL_REQUIRE(physicscommon::get_dz(ggr,Rd,100.0,0.0)==100.0);
      physicscommon::get_dz(team,psuedo_density,pressure,temperature,qv,dz);
      physicscommon::get_dz(team,psuedo_density_packed,pressure_packed,temperature_packed,qv_packed,dz_packed);
      // z_int property tests
      // get_z_int should work for Real and for Spack
      physicscommon::get_z_int(team,num_levs,dz_for_testing,z_int);
      physicscommon::get_z_int(team,num_levs,dz_for_testing_packed,z_int_packed);

    }, nerr); // Kokkos parallel_reduce "test_universal_physics"

    Kokkos::fence();
    auto exner_host           = Kokkos::create_mirror(exner);
    auto theta_host           = Kokkos::create_mirror(theta);
    auto T_mid_from_pot_host  = Kokkos::create_mirror(T_mid_from_pot);
    auto T_virtual_host       = Kokkos::create_mirror(T_virtual);
    auto T_mid_from_virt_host = Kokkos::create_mirror(T_mid_from_virt);
    auto dse_host             = Kokkos::create_mirror(dse);
    auto z_int_host           = Kokkos::create_mirror(z_int);
    auto dz_host              = Kokkos::create_mirror(dz);
    auto exner_pack_host           = Kokkos::create_mirror(exner_packed);
    auto theta_pack_host           = Kokkos::create_mirror(theta_packed);
    auto T_mid_from_pot_pack_host  = Kokkos::create_mirror(T_mid_from_pot_packed);
    auto T_virtual_pack_host       = Kokkos::create_mirror(T_virtual_packed);
    auto T_mid_from_virt_pack_host = Kokkos::create_mirror(T_mid_from_virt_packed);
    auto dse_pack_host             = Kokkos::create_mirror(dse_packed);
    auto z_int_pack_host           = Kokkos::create_mirror(z_int_packed);
    auto dz_pack_host              = Kokkos::create_mirror(dz_packed);
    Kokkos::deep_copy(exner_host           , exner);
    Kokkos::deep_copy(theta_host           , theta);
    Kokkos::deep_copy(T_mid_from_pot_host  , T_mid_from_pot);
    Kokkos::deep_copy(T_virtual_host       , T_virtual);
    Kokkos::deep_copy(T_mid_from_virt_host , T_mid_from_virt);
    Kokkos::deep_copy(dse_host             , dse);
    Kokkos::deep_copy(z_int_host           , z_int);
    Kokkos::deep_copy(dz_host              , dz);
    Kokkos::deep_copy(exner_pack_host           , exner_packed);
    Kokkos::deep_copy(theta_pack_host           , theta_packed);
    Kokkos::deep_copy(T_mid_from_pot_pack_host  , T_mid_from_pot_packed);
    Kokkos::deep_copy(T_virtual_pack_host       , T_virtual_packed);
    Kokkos::deep_copy(T_mid_from_virt_pack_host , T_mid_from_virt_packed);
    Kokkos::deep_copy(dse_pack_host             , dse_packed);
    Kokkos::deep_copy(z_int_pack_host           , z_int_packed);
    Kokkos::deep_copy(dz_pack_host              , dz_packed);
    //Kokkos::parallel_for(Kokkos::TeamThreadRange(team,num_levs+1), [&] (const int k)
    for(int k=0;k<num_levs;k++)
    {
     int ipack = k / Spack::n;
     int ivec  = k % Spack::n;
      // Make sure all columnwise results don't contain any obvious errors:
      // exner
      EKAT_REQUIRE(exner_host(k)==exner_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(exner_host(k)));
      EKAT_REQUIRE(!(exner_host(k)<0));
      // potential temperature
      EKAT_REQUIRE(theta_host(k)==theta_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(theta_host(k)));
      EKAT_REQUIRE(!(theta_host(k)<0));
      // temperature from potential temperature
      EKAT_REQUIRE(T_mid_from_pot_host(k)==T_mid_from_pot_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(T_mid_from_pot_host(k)));
      EKAT_REQUIRE(!(T_mid_from_pot_host(k)<0));
      // virtual temperature
      EKAT_REQUIRE(T_virtual_host(k)==T_virtual_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(T_virtual_host(k)));
      EKAT_REQUIRE(!(T_virtual_host(k)<0));
      // temperature from virtual temperature
      EKAT_REQUIRE(T_mid_from_virt_host(k)==T_mid_from_virt_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(T_mid_from_virt_host(k)));
      EKAT_REQUIRE(!(T_mid_from_virt_host(k)<0));
      // DSE
      EKAT_REQUIRE(dse_host(k)==dse_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(dse_host(k)));
      EKAT_REQUIRE(!(dse_host(k)<0));
      // dz
      EKAT_REQUIRE(dz_host(k)==dz_pack_host(ipack)[ivec]);
      EKAT_REQUIRE(!isnan(dz_host(k)));
      EKAT_REQUIRE(!(dz_host(k)<=0));
      // z_int
      EKAT_REQUIRE_MSG(z_int_host(k)==z_int_pack_host(ipack)[ivec],"ERROR: " + std::to_string(z_int_host(k)) + " != " + std::to_string(z_int_pack_host(ipack)[ivec]) + ", k = " + std::to_string(k));
      const auto k_bwd = num_levs-k;
      EKAT_REQUIRE(z_int_host(k_bwd)==k*(k+1)/2);
      EKAT_REQUIRE(!isnan(z_int_host(k)));
      EKAT_REQUIRE(!(z_int_host(k)<0));
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
