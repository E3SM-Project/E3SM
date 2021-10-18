#include "catch2/catch.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_common_physics_functions.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace{

template<typename ScalarT, int NumLevels>
struct ChecksHelpers {

  static bool is_non_negative (const ScalarT& s, const int k) {
    return not ( k<NumLevels && (s<0 || isnan(s)) );
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return lhs==rhs;
  }
  static bool approx_equal (const ScalarT lhs, const ScalarT rhs,
                            const int k, const ScalarT tol) {
    using std::abs;
    return not ( k<NumLevels && abs(lhs-rhs)>=tol );
  }
  static bool approx_equal (const ScalarT computed, const ScalarT expected, const ScalarT tol) {
    using std::abs;
    return abs(computed-expected)/abs(expected) < tol;
  }
};

template<typename T, int N, int NumLevels>
struct ChecksHelpers<ekat::Pack<T,N>,NumLevels> {
  using ScalarT = ekat::Pack<T,N>;

  static bool is_non_negative (const ScalarT& s, const int k) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return ( range_mask && (s<0 || isnan(s) ) ).none();
  }
  static bool equal (const ScalarT& lhs, const ScalarT& rhs) {
    return (lhs==rhs).all();
  }
  static bool approx_equal (const ScalarT& lhs, const ScalarT& rhs,
                            const int k, const T tol) {
    const auto range = ekat::range<ScalarT>(k*N);
    const auto range_mask = range < NumLevels;
    return (range_mask && abs(lhs-rhs)>=tol).none();
  }
  static bool approx_equal (const ScalarT& computed, const ScalarT& expected, const T tol) {
    return (abs(computed-expected)/abs(expected) < tol).all();
  }
};

// Helper function. Create Mirror View and Deep-Copy (CMVDC)
template<typename ViewT>
auto cmvdc (const ViewT& v_d) -> typename ViewT::HostMirror {
  auto v_h = Kokkos::create_mirror_view(v_d);
  Kokkos::deep_copy(v_h,v_d);
  return v_h;
}

//-----------------------------------------------------------------------------------------------//
template<typename ScalarT, typename DeviceT>
void run(std::mt19937_64& engine)
{
  using STraits    = ekat::ScalarTraits<ScalarT>;
  using RealType   = typename STraits::scalar_type;

  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<RealType>;

  using KT         = ekat::KokkosTypes<DeviceT>;
  using ExecSpace  = typename KT::ExeSpace;
  using TeamPolicy = typename KT::TeamPolicy;
  using MemberType = typename KT::MemberType;
  using view_1d    = typename KT::template view_1d<ScalarT>;
  using rview_1d   = typename KT::template view_1d<RealType>;

  static constexpr auto Rd       = PC::RD;
  static constexpr auto inv_cp   = PC::INV_CP;
  static constexpr auto g        = PC::gravit;
  static constexpr auto test_tol = PC::macheps*1e3;

  constexpr int pack_size = sizeof(ScalarT) / sizeof(RealType);
  using pack_info = ekat::PackInfo<pack_size>;

  constexpr int num_levs = 100; // Number of levels to use for tests.
  const     int num_mid_packs = pack_info::num_packs(num_levs);
  const     int num_int_packs = pack_info::num_packs(num_levs+1);

  using Check = ChecksHelpers<ScalarT,num_levs>;

  // Input (randomized) views
  view_1d temperature("temperature",num_mid_packs),
          height("height",num_mid_packs),
          qv("qv",num_mid_packs),
          pressure("pressure",num_mid_packs),
          pseudo_density("pseudo_density",num_mid_packs),
          dz_for_testing("dz_for_testing",num_mid_packs),
          mmr_for_testing("mass_mixing_ratio_for_testing",num_mid_packs),
          wetmmr_for_testing("wet_mass_mixing_ratio_for_testing",num_mid_packs);
  // Output views
  view_1d exner("exner",num_mid_packs),
          theta("theta",num_mid_packs),
          T_from_Theta("T_from_Theta",num_mid_packs),
          Tv("T_virtual",num_mid_packs),
          T_from_Tv("T_from_T_virtual",num_mid_packs),
          dse("dse",num_mid_packs),
          dz("dz",num_mid_packs),
          z_int("z_int",num_int_packs),
          vmr("volume_mixing_ratio",num_mid_packs),
          mmr("mass_mixing_ratio",num_mid_packs),
          wetmmr("wet_mass_mixing_ratio",num_mid_packs),
          drymmr("dry_mass_mixing_ratio",num_mid_packs),
          density("density",num_mid_packs);

  auto dview_as_real = [&] (const view_1d& v) -> rview_1d {
    return rview_1d(reinterpret_cast<RealType*>(v.data()),v.size()*pack_size);
  };
  auto hview_as_real = [&] (const typename view_1d::HostMirror& v) -> typename rview_1d::HostMirror {
    return typename rview_1d::HostMirror(reinterpret_cast<RealType*>(v.data()),v.size()*pack_size);
  };

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_qv(1e-6,1e-3),
       pdf_dp(1.0,100.0),
       pdf_pres(0.0,PC::P0),
       pdf_temp(200.0,400.0),
       pdf_height(0.0,1e5),
       pdf_dz(1.0,1e5),
       pdf_surface(100.0,400.0),
       pdf_mmr(0,0.99);

  //contruct random integers
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf_rand_int(1,100);

  ekat::genRandArray(dview_as_real(temperature),     engine,pdf_temp);
  ekat::genRandArray(dview_as_real(height),          engine,pdf_height);
  ekat::genRandArray(dview_as_real(qv),              engine,pdf_qv);
  ekat::genRandArray(dview_as_real(pressure),        engine,pdf_pres);
  ekat::genRandArray(dview_as_real(pseudo_density),  engine,pdf_dp);
  ekat::genRandArray(dview_as_real(mmr_for_testing), engine,pdf_mmr);
  ekat::genRandArray(dview_as_real(wetmmr_for_testing), engine,pdf_mmr);

  // Construct a simple set of `dz` values for testing the z_int function
  auto dz_for_testing_host = Kokkos::create_mirror_view(dz_for_testing);
  auto dz_real_host = hview_as_real(dz_for_testing_host);
  for (int k = 0; k<num_mid_packs*pack_size; ++k) {
    dz_real_host[k] = num_levs-k;
  }
  Kokkos::deep_copy(dz_for_testing,dz_for_testing_host);

  // ---- Single-scalar functions tests ---- //

  const ScalarT p0   = PC::P0;
  const ScalarT zero = 0.0;
  const ScalarT one  = 1.0;

  ScalarT p, T0, theta0, tmp, qv0, dp0, mmr0, vmr0, wetmmr0, dz0, Tv0, rho0, rand_int0;
  RealType surf_height;

  // calculate density property tests:
  //  - calculate_density(pseudo_density=zero) should return 0.0
  //  - density using ideal gas law should match the calculated density.
  dz0  = pdf_dz(engine);
  p    = pdf_pres(engine);
  dp0  = pdf_dp(engine);
  T0   = pdf_temp(engine);
  qv0  = pdf_qv(engine);
  Tv0  = PF::calculate_virtual_temperature(T0,qv0);
  dz0  = PF::calculate_dz(dp0,p,T0,qv0);
  rho0 = p / Tv0 / Rd;  // Ideal gas law
  rand_int0 = pdf_rand_int(engine);//random integers

  REQUIRE( Check::equal(PF::calculate_density(zero,dz0),zero) );
  REQUIRE( Check::approx_equal(PF::calculate_density(dp0,dz0),rho0,test_tol) );

  // Exner property tests:
  //  - exner_function(p0) should return 1.0
  //  - exner_function(0.0) should return 0.0
  //  - exner_function(2*x) should return 2**(Rd/cp)*exner_function(x)
  REQUIRE( Check::equal(PF::exner_function(p0),one) );
  REQUIRE( Check::equal(PF::exner_function(zero),zero) );

  p = pdf_pres(engine);
  const auto exner1 = PF::exner_function(p);
  const auto exner2 = PF::exner_function(p/2);
  const auto factor = pow(2.0,Rd*inv_cp);
  REQUIRE( Check::approx_equal(exner1,factor*exner2,test_tol) );

  // Potential temperature property tests:
  //  - theta=T when p=p0
  //  - theta(T=0) = 0
  //  - T(theta=0) = 0
  //  - T(theta(T0)) = T0  (in exact arithmetic)
  //  - theta(T(theta0)) = theta0 (in exact arithmetic)
  //  - T_from_theta ant theta_from_T are one the inverse of the other (up to roundoff errors)
  T0     = pdf_temp(engine);
  theta0 = pdf_temp(engine);
  p      = pdf_pres(engine);
  REQUIRE( Check::equal(PF::calculate_theta_from_T(T0,p0),T0) );
  REQUIRE( Check::equal(PF::calculate_theta_from_T(zero,p),zero) );
  REQUIRE( Check::equal(PF::calculate_T_from_theta(theta0,p0),theta0) );
  REQUIRE( Check::equal(PF::calculate_T_from_theta(zero,p),zero) );
  REQUIRE( Check::approx_equal(PF::calculate_T_from_theta(PF::calculate_theta_from_T(T0,p),p),T0,test_tol) );
  REQUIRE( Check::approx_equal(PF::calculate_theta_from_T(PF::calculate_T_from_theta(theta0,p),p),theta0,test_tol) );

  // Virtual temperature property tests:
  //  - T_virt(T=0) = 0.0
  //  - T_virt(T=T0,qv=0) = T0
  //  - T(T_virt=0) = 0.0
  //  - T(T_virt=T0,qv=0) = T0
  //  - T_virt(T=T0) = T0
  //  - T(T_virt=T0) = T0
  //  - T_from_Tvirt and Tvirt_from_T are one the inverse of the other (up to roundoff errors)
  qv0 = pdf_qv(engine);
  T0  = pdf_temp(engine);
  REQUIRE( Check::equal(PF::calculate_virtual_temperature(zero,qv0),zero) );
  REQUIRE( Check::equal(PF::calculate_virtual_temperature(T0,zero),T0) );
  REQUIRE( Check::equal(PF::calculate_temperature_from_virtual_temperature(zero,qv0),zero) );
  REQUIRE( Check::equal(PF::calculate_temperature_from_virtual_temperature(T0,zero),T0) );

  tmp = PF::calculate_temperature_from_virtual_temperature(T0,qv0);
  REQUIRE( Check::approx_equal(PF::calculate_virtual_temperature(tmp,qv0),T0,test_tol) );
  tmp = PF::calculate_virtual_temperature(T0,qv0);
  REQUIRE( Check::approx_equal(PF::calculate_temperature_from_virtual_temperature(tmp,qv0),T0,test_tol) );

  // DSE property tests:
  //  - calculate_dse(T=0.0, z=0.0) = surf_height
  //  - calculate_dse(T=1/cp, z=1/gravity) = surf_height+2
  surf_height = pdf_surface(engine);
  REQUIRE( Check::equal(PF::calculate_dse(zero,zero,surf_height),ScalarT(surf_height)) );
  REQUIRE( Check::equal(PF::calculate_dse(ScalarT(inv_cp),ScalarT(1/g),surf_height),ScalarT(surf_height+2.0)) );

  // WETMMR to DRYMMR (and vice versa) property tests
  wetmmr0 = pdf_mmr(engine);// get initial inputs for wetmmr_from_drymmr and drymmr_from_wetmmr functions
  qv0  = pdf_qv(engine);  // This is an input for mmr_tests, so it won't be modified by mmr tests
  // mmr_test1: For zero drymmr, wetmmr should be zero
  // mmr_test2: For zero wetmmr, drymmr should be zero
  // mmr_test3: Compute drymmr from wetmmr0 and then use the result to compute wetmmr, which should be approximately
  //            equal to wetmmr0

  // mmr_test4: [to test mathematical properties of the function]
  //            Compute drymmr from wetmmr0 and then use the result to
  //            compute wetmmr. "qv" should be in the form of 2^k+1 (1 is added as we have [1-qv] in the numerator or
  //            denominator of this function), so that the results of wet->dry->wet are exactly
  //            the same as the initial input mmr(wetmmr0). This is a mathematical property test, so unphysical invalid qv
  //            values are also acceptable
  //
  //            *WARNING* This test might fail if a check is added to this function to accept only valid values for qv0

  REQUIRE( Check::equal(PF::calculate_wetmmr_from_drymmr(zero,qv0),zero) ); //mmr_test1
  REQUIRE( Check::equal(PF::calculate_drymmr_from_wetmmr(zero,qv0),zero) ); //mmr_test2

  //mmr_test3
  tmp = PF::calculate_drymmr_from_wetmmr(wetmmr0,qv0);//get drymmr from wetmmr0
  tmp = PF::calculate_wetmmr_from_drymmr(tmp, qv0);//convert it back to wetmmr0
  REQUIRE( Check::approx_equal(tmp,wetmmr0,test_tol) );// wetmmr0 should be equal to tmp

  //mmr_test4
  ScalarT qv0_2k_m1 = pow(2,rand_int0)+1; // 2^rand_int+1
  tmp = PF::calculate_drymmr_from_wetmmr(wetmmr0,qv0_2k_m1);//get drymmr from wetmmr0
  tmp = PF::calculate_wetmmr_from_drymmr(tmp, qv0_2k_m1);//convert it back to wetmmr0
  REQUIRE( Check::equal(tmp,wetmmr0) );// wetmmr0 should be exactly equal to tmp

  // DZ property tests:
  //  - calculate_dz(pseudo_density=0) = 0
  //  - calculate_dz(T=0) = 0
  //  - calculate_dz(pseudo_density=p0,p_mid=p0,T=1.0,qv=0) = Rd/g
  //  - calculate_dz(pseudo_density=g,p_mid=Rd,T=T0,qv=0) = T0
  p = pdf_pres(engine);
  T0 = pdf_temp(engine);
  qv0 = pdf_qv(engine);
  dp0 = pdf_dp(engine);
  REQUIRE( Check::equal(PF::calculate_dz(zero,p,T0,qv0),zero) );
  REQUIRE( Check::equal(PF::calculate_dz(dp0,p,zero,qv0),zero) );
  REQUIRE( Check::equal(PF::calculate_dz(ScalarT(p0),ScalarT(p0),one,zero),ScalarT(Rd/g)) );
  REQUIRE( Check::approx_equal(PF::calculate_dz(ScalarT(g),ScalarT(Rd),T0,zero),T0,test_tol) );

  // MMR and VMR property tests:
  //  - calculate_vmr_from_mmr(mmr=0) = 0
  //  - calculate_mmr_from_vmr(vmr=0) = 0
  //  - calculate_vmr_from_mmr(calculate_mmr_from_vmr(gas_name="h2o",vmr0)) = vmr
  //  - calculate_vmr_from_mmr(calculate_mmr_from_vmr(gas_name="o2",vmr0)) = vmr, test that changing gas name changes the result.
  //  - calculate_mmr_from_vmr(calculate_vmr_from_mmr(gas_name="h2o",mmr0)) = vmr
  //  - calculate_mmr_from_vmr(calculate_vmr_from_mmr(gas_name="o2",mmr0)) != vmr, test that changing gas name changes the result.
  mmr0 = pdf_mmr(engine);
  vmr0 = pdf_mmr(engine);
  qv0 = pdf_qv(engine);
  const auto h2o_mol = PC::get_gas_mol_weight("h2o");
  const auto o2_mol  = PC::get_gas_mol_weight("o2");
  REQUIRE( Check::equal(PF::calculate_vmr_from_mmr(h2o_mol,qv0,zero),zero) );
  REQUIRE( Check::equal(PF::calculate_mmr_from_vmr(h2o_mol,qv0,zero),zero) );
  tmp = PF::calculate_vmr_from_mmr(h2o_mol,qv0,mmr0);
  REQUIRE( Check::approx_equal(PF::calculate_mmr_from_vmr(h2o_mol,qv0,tmp),mmr0,test_tol) );
  REQUIRE( !Check::approx_equal(PF::calculate_mmr_from_vmr(o2_mol,qv0,tmp),mmr0,test_tol) );
  tmp = PF::calculate_mmr_from_vmr(h2o_mol,qv0,vmr0);
  REQUIRE( Check::approx_equal(PF::calculate_vmr_from_mmr(h2o_mol,qv0,tmp),vmr0,test_tol) );
  REQUIRE( !Check::approx_equal(PF::calculate_vmr_from_mmr(o2_mol,qv0,tmp),vmr0,test_tol) );

  // --------- Run tests on full columns of data ----------- //
  TeamPolicy policy(ekat::ExeSpaceUtils<ExecSpace>::get_default_team_policy(1, 1));
  Kokkos::parallel_for("test_universal_physics", policy, KOKKOS_LAMBDA(const MemberType& team) {

    // Compute density(dp,dz)
    PF::calculate_density(team,pseudo_density,dz_for_testing,density);

    // Compute exner(p)
    PF::exner_function(team,pressure,exner);

    // Compute theta(T), and T(theta(T))
    PF::calculate_theta_from_T(team,temperature,pressure,theta);
    PF::calculate_T_from_theta(team,theta,pressure,T_from_Theta);

    // Compute Tv(T,qv), and T(Tv(T,qv),qv)
    PF::calculate_virtual_temperature(team,temperature,qv,Tv);
    PF::calculate_temperature_from_virtual_temperature(team,Tv,qv,T_from_Tv);

    // Compute dse(T,z,z_surf)
    PF::calculate_dse(team,temperature,height,surf_height,dse);

    // Compute dz(dp,p,T,qv)
    PF::calculate_dz(team,pseudo_density,pressure,temperature,qv,dz);

    // Compute z_int(dz,z_surf)
    PF::calculate_z_int(team,num_levs,dz_for_testing,surf_height,z_int);

    // Compute vmr from mmr and vice versa
    PF::calculate_vmr_from_mmr(team,h2o_mol,qv,mmr_for_testing,vmr);
    PF::calculate_mmr_from_vmr(team,h2o_mol,qv,vmr,mmr);

    // Compute drymmr from wetmmr
    PF::calculate_drymmr_from_wetmmr(team,wetmmr_for_testing,qv,drymmr);

    // Convert drymmr computed above to wetmmr
    PF::calculate_wetmmr_from_drymmr(team,drymmr,qv,wetmmr);

  }); // Kokkos parallel_for "test_universal_physics"
  Kokkos::fence();

  // Deep copy to host, and check the properties of the full view output
  auto temperature_host     = cmvdc(temperature);
  auto theta_host           = cmvdc(theta);
  auto pressure_host        = cmvdc(pressure);
  auto qv_host              = cmvdc(qv);

  auto density_host         = cmvdc(density);
  auto exner_host           = cmvdc(exner);
  auto T_from_Theta_host    = cmvdc(T_from_Theta);
  auto Tv_host              = cmvdc(Tv);
  auto T_from_Tv_host       = cmvdc(T_from_Tv);
  auto dse_host             = cmvdc(dse);
  auto z_int_host           = cmvdc(z_int);
  auto dz_host              = cmvdc(dz);
  auto vmr_host             = cmvdc(vmr);
  auto mmr_host             = cmvdc(mmr);
  auto mmr_for_testing_host = cmvdc(mmr_for_testing);
  auto wetmmr_host          = cmvdc(wetmmr);
  auto drymmr_host          = cmvdc(drymmr);
  auto wetmmr_for_testing_host = cmvdc(wetmmr_for_testing);

  for (int k=0; k<num_mid_packs; ++k) {

    // Make sure all results don't contain invalid numbers
    REQUIRE( Check::is_non_negative(density_host(k),k) );
    REQUIRE( Check::is_non_negative(exner_host(k),k) );
    REQUIRE( Check::is_non_negative(theta_host(k),k) );
    REQUIRE( Check::is_non_negative(T_from_Theta_host(k),k) );
    REQUIRE( Check::is_non_negative(T_from_Tv_host(k),k) );
    REQUIRE( Check::is_non_negative(Tv_host(k),k) );
    REQUIRE( Check::is_non_negative(dz_host(k),k) );
    REQUIRE( Check::is_non_negative(z_int_host(k),k) );
    REQUIRE( Check::is_non_negative(dse_host(k),k) );
    REQUIRE( Check::is_non_negative(dse_host(k),k) );
    REQUIRE( Check::is_non_negative(vmr_host(k),k) );
    REQUIRE( Check::is_non_negative(mmr_host(k),k) );
    REQUIRE( Check::is_non_negative(wetmmr_host(k),k) );
    REQUIRE( Check::is_non_negative(drymmr_host(k),k) );

    // Check T(Theta(T))==T (up to roundoff tolerance)
    REQUIRE ( Check::approx_equal(T_from_Theta_host(k),temperature_host(k),k,test_tol) );

    // Check T(T_virtual(T))==T (up to roundoff tolerance)
    REQUIRE ( Check::approx_equal(T_from_Tv_host(k),temperature_host(k),k,test_tol) );

    // Check vmr(mmr(vmr))==mmr (up to roundoff tolerance)
    REQUIRE ( Check::approx_equal(mmr_host(k),mmr_for_testing_host(k),test_tol) );

    // Check wetmmr(drymmr(wetmmr))==wetmmr (up to roundoff tolerance)
    REQUIRE ( Check::approx_equal(wetmmr_host(k),wetmmr_for_testing_host(k),test_tol) );
  }
} // run()

TEST_CASE("common_physics_functions_test", "[common_physics_functions_test]"){
  // Run tests for both Real and Pack, and for (potentially) different pack sizes
  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  std::random_device rdev;
  using rngAlg = std::mt19937_64;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rdev() : catchRngSeed;
  // Print seed to screen to trace tests that fail.
  std::cout << " -> Random number generator seed:: " << seed << "\n";
  if (catchRngSeed==0) {
    std::cout << "    Note: catch rng seed was 0 (default). We interpret that as a request to pick a random seed.\n"
                 "    To reproduce a previous run, use --rng-seed N to provide the rng seed.\n\n";
  }
  std::cout << " -> Nnumber of randomized runs: " << num_runs << "\n\n";
  rngAlg engine(seed);

  printf(" -> Testing Real scalar type...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Real,Device>(engine);
  }
  printf("ok!\n");

  printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_SMALL_PACK_SIZE);
  for (int irun=0; irun<num_runs; ++irun) {
    run<ekat::Pack<Real,SCREAM_SMALL_PACK_SIZE>,Device>(engine);
  }
  printf("ok!\n");

  if (SCREAM_PACK_SIZE!=SCREAM_SMALL_PACK_SIZE) {
    printf(" -> Testing Pack<Real,%d> scalar type...",SCREAM_PACK_SIZE);
    for (int irun=0; irun<num_runs; ++irun) {
      run<ekat::Pack<Real,SCREAM_PACK_SIZE>,Device>(engine);
    }
    printf("ok!\n");
  }

  printf("\n");

} // TEST_CASE

} // namespace
