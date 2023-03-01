#include "catch2/catch.hpp"

#include "physics/share/physics_constants.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_common_physics_functions.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_test_utils.hpp"

#include <iomanip>

namespace{

template<typename ScalarT, int NumLevels>
struct ChecksHelpers {

  static bool is_non_negative (const ScalarT& s, const int k) {
    return not ( k<NumLevels && (s<0 || std::isnan(s)) );
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

template<typename DeviceT>
void run_scalar_valued_fns(std::mt19937_64& engine)
{
  /*
  Most of the common physics functions are templated to operate on scalars or on packs of vertical indices. 
  The functions tested here don't include any vertical dimension (e.g. they handle variables only defined 
  at the surface), so are only defined for scalar reals. This is fundamentally different than the other 
  functions, so these functions get their own run test.
  */	

  using RealType   = scream::Real;
  using PF         = scream::PhysicsFunctions<DeviceT>;
  using PC         = scream::physics::Constants<RealType>;
  using Check = ChecksHelpers<RealType,1>; //1 is for number of levels.

  static constexpr auto pi       = PC::Pi;
  static constexpr auto Rd       = PC::RD;
  static constexpr auto g        = PC::gravit;
  static constexpr auto test_tol = PC::macheps*1e3;
  static constexpr auto coeff_1  = PC::earth_ellipsoid1;
  static constexpr auto coeff_2  = PC::earth_ellipsoid2;
  static constexpr auto coeff_3  = PC::earth_ellipsoid3;

  // Construct random input data
  using RPDF = std::uniform_real_distribution<RealType>;
  RPDF pdf_lat(-pi/2.0,pi/2.0),
       pdf_area(1e-8,pi*pi);

  //calculate_surface_air_T property tests:
  // If z_mid_bot==0, output should equal T_mid_bot
  // If z_mid_bot>0, output should be warmer that T_mid_bot
  // If z_mid_bot<0, output should be colder than T_mid_bot
  // if z_mid_bot is 1 km, output should be exactly 6.5 K warmer
  REQUIRE( Check::approx_equal(PF::calculate_surface_air_T(300,0),300,test_tol) );
  REQUIRE( Check::approx_equal(PF::calculate_surface_air_T(300,1000),306.5,test_tol) );
  REQUIRE( Check::approx_equal(PF::calculate_surface_air_T(250.3,-10.2),250.2337,test_tol) );
  
  // lapse_T_for_psl property tests:
  // If T_ground = 0, T_ground_tmp should be 255/2 and lapse should be 0.0065 Really cold case.
  // If T_ground = 300 K with phi_ground>0 m2/s2, lapse = 0 and T_ground_tmp=0.5*(290.5+T_ground). Really hot case.
  // If T_ground = 290 K and phi_ground=10000 m2/s2, T_ground_tmp=T_ground and lapse
  //    is such that T_sl=T_ground+lapse*phi_ground/gravit is within roundoff of 290.5 K. Marginally hot case.
  // If T_ground = 280 K and phi_ground=100 m2/s2, T_ground_tmp=T_ground and lapse=6.5 K/km (typical conditions)
  RealType lapse;
  RealType T_ground_tmp;
  RealType T_ground=0;
  RealType phi_ground=100;
  PF::lapse_T_for_psl(T_ground, phi_ground, lapse, T_ground_tmp );
  REQUIRE( Check::approx_equal(T_ground_tmp,0.5*255.,test_tol) );
  REQUIRE( Check::equal(lapse,0.0065) );
  T_ground=300;
  PF::lapse_T_for_psl(T_ground, phi_ground, lapse, T_ground_tmp );
  REQUIRE( Check::approx_equal(T_ground_tmp,0.5*(290.5+T_ground),test_tol) );
  REQUIRE( Check::equal(lapse,0) );
  T_ground=290;
  phi_ground=10000;
  PF::lapse_T_for_psl(T_ground, phi_ground, lapse, T_ground_tmp );
  REQUIRE( Check::equal(T_ground_tmp, T_ground) );
  REQUIRE( Check::approx_equal( T_ground_tmp+lapse*phi_ground/g, 290.5 , test_tol) );
  T_ground=280;
  phi_ground=100;
  PF::lapse_T_for_psl(T_ground, phi_ground, lapse, T_ground_tmp );
  REQUIRE( Check::equal(T_ground_tmp, T_ground) );
  REQUIRE( Check::equal(lapse,0.0065) );
  
  // PSL property tests:
  // PSL==surface pressure if surface height = 0 m
  // computed value close to exact solution when lapse rate is zero (i.e. very warm conditions)
  // computed value close to exact solution when lapse rate is 6.5 K/km (typical conditions)
  // PSL lower than surface pressure whenever surface height > 0 m
  // PSL greater than surface pressure whenever surface height < 0 m
  RealType p_ground=100000; //can't use p0 because that is a ScalarT which could be a pack
  REQUIRE( Check::equal(PF::calculate_psl(310 , p_ground, 0 ), p_ground) );
  REQUIRE( Check::equal(PF::calculate_psl( 2 , 2, 0 ), 2) );
  
  T_ground=300;
  T_ground_tmp=0.5*(290.5+T_ground);
  RealType psl_exact = p_ground*std::exp(phi_ground/(Rd*T_ground_tmp));
  RealType psl=PF::calculate_psl( T_ground , p_ground, phi_ground );
  REQUIRE( Check::approx_equal( psl_exact, psl, test_tol) );
  REQUIRE( psl>p_ground);
  
  T_ground=280;
  phi_ground=-100;
  lapse=0.0065;
  psl_exact = p_ground*std::pow( 1. + lapse/g*phi_ground/T_ground, g/(Rd*lapse));
  psl=PF::calculate_psl( T_ground , p_ground, phi_ground );
  REQUIRE( Check::approx_equal( psl_exact, psl, test_tol) );
  REQUIRE( psl<p_ground );

  // Get dx from grid cell area property tests:
  RealType area, lat;
  area = 0.0; lat = 1.0;
  REQUIRE( Check::equal(PF::calculate_dx_from_area(area,lat),0.0) );
  area = (pi/180.0)*(pi/180.0);
  // Note, it is assumed that lat is in degrees
  lat = 0.0;
  REQUIRE( Check::equal(PF::calculate_dx_from_area(area,lat), coeff_1-coeff_2+coeff_3) );
  lat = 90.0;
  REQUIRE( Check::equal(PF::calculate_dx_from_area(area,lat), coeff_1+coeff_2+coeff_3) );
  lat = 45.0;
  REQUIRE( Check::approx_equal(PF::calculate_dx_from_area(area,lat), coeff_1 - coeff_3, test_tol) );
  lat = 22.5;
  REQUIRE( Check::approx_equal(PF::calculate_dx_from_area(area,lat), coeff_1-std::sqrt(2.0)/2.0*coeff_2, test_tol) );
  lat     = pdf_lat(engine);
  area    = pdf_area(engine);
  REQUIRE( Check::equal(PF::calculate_dx_from_area(area,lat),PF::calculate_dx_from_area(area,-lat)) );

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
  static constexpr auto cp       = PC::CP;
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
          qv_dry("qv_dry",num_mid_packs),
          qv_wet("qv_wet",num_mid_packs),
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
          T_from_dse("T_from_dse",num_mid_packs),
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
       pdf_mmr(0,0.99),
       pdf_dse(1e5,5e5);

  //contruct random integers
  using IPDF = std::uniform_int_distribution<int>;
  IPDF pdf_rand_int(1,100);

  ekat::genRandArray(dview_as_real(temperature),     engine,pdf_temp);
  ekat::genRandArray(dview_as_real(height),          engine,pdf_height);
  ekat::genRandArray(dview_as_real(qv),              engine,pdf_qv);
  ekat::genRandArray(dview_as_real(qv_wet),          engine,pdf_qv);
  ekat::genRandArray(dview_as_real(qv_dry),          engine,pdf_qv);
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

  ScalarT p, T0, theta0, tmp, qv0, dp0, mmr0, vmr0, dz0, Tv0, rho0, z0, dse0;
  ScalarT wetmmr0, drymmr0, qv_dry0, qv_wet0;
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
  //  - calculate_temperature_from_dse(dse=cp, z=cp/gravity) = -1.0*surf_height/cp (up to roundoff errors)
  //  - calculate_dse and calculate_temperature_from_dse are one the inverse of the other (up to roundoff errors)
  surf_height = pdf_surface(engine);
  z0          = pdf_height(engine);
  dse0        = pdf_dse(engine);
  T0          = pdf_temp(engine);
  REQUIRE( Check::equal(PF::calculate_dse(zero,zero,surf_height),ScalarT(surf_height)) );
  REQUIRE( Check::equal(PF::calculate_dse(ScalarT(inv_cp),ScalarT(1/g),surf_height),ScalarT(surf_height+2.0)) );
  REQUIRE( Check::approx_equal(PF::calculate_temperature_from_dse(ScalarT(cp),ScalarT(cp/g),surf_height),
                               ScalarT(-1.0*surf_height/cp), test_tol) );
  REQUIRE( Check::approx_equal(PF::calculate_dse(PF::calculate_temperature_from_dse(dse0,z0,surf_height),
                               z0,surf_height),dse0,test_tol) );
  REQUIRE( Check::approx_equal(PF::calculate_temperature_from_dse(PF::calculate_dse(T0,z0,surf_height),
                               z0,surf_height),T0,test_tol) );

  // WETMMR to DRYMMR (and vice versa) property tests
  wetmmr0 = pdf_mmr(engine);// get initial wet mmr
  drymmr0 = pdf_mmr(engine);// get initial dry mmr
  qv_wet0 = pdf_qv(engine); // get initial qv in wet mmr
  qv_dry0 = pdf_qv(engine); // get initial qv in dry mmr
  // mmr_test1: For zero drymmr, wetmmr should be zero
  // mmr_test2: For zero wetmmr, drymmr should be zero
  // mmr_test3: Compute drymmr from wetmmr0 and then use the result to compute wetmmr, which should be approximately
  //            equal to wetmmr0. NOTE: calculate_wetmmr_from_drymmr takes qv in dry mmr as an argument and
  //            calculate_drymmr_from_wetmmr takes qv in wet mmr as an argument. Therefore, we need to convert
  //            qv from wet mmr to dry mmr as well as part of these conversions

  REQUIRE( Check::equal(PF::calculate_wetmmr_from_drymmr(zero,qv_dry0),zero) ); //mmr_test1
  REQUIRE( Check::equal(PF::calculate_drymmr_from_wetmmr(zero,qv_wet0),zero) ); //mmr_test2

  //mmr_test3
  drymmr0 = PF::calculate_drymmr_from_wetmmr(wetmmr0,qv_wet0);//get drymmr from wetmmr0 using qv_wet0
  // Now convert qv wet mmr to qv dry mmr as qv dry mmr is an input for the "calculate_wetmmr_from_drymmr" function
  qv_dry0 = PF::calculate_drymmr_from_wetmmr(qv_wet0, qv_wet0);
  tmp     = PF::calculate_wetmmr_from_drymmr(drymmr0, qv_dry0);//convert it back to wetmmr0
  REQUIRE( Check::approx_equal(tmp,wetmmr0,test_tol) );// wetmmr0 should be equal to tmp

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

    // Compute temperature from dse (dse,z,z_surf)
    PF::calculate_temperature_from_dse(team,dse,height,surf_height,T_from_dse);

    // Compute dz(dp,p,T,qv)
    PF::calculate_dz(team,pseudo_density,pressure,temperature,qv,dz);

    // Compute z_int(dz,z_surf)
    PF::calculate_z_int(team,num_levs,dz_for_testing,surf_height,z_int);

    // Compute vmr from mmr and vice versa
    PF::calculate_vmr_from_mmr(team,h2o_mol,qv,mmr_for_testing,vmr);
    PF::calculate_mmr_from_vmr(team,h2o_mol,qv,vmr,mmr);

    // Compute drymmr from wetmmr
    PF::calculate_drymmr_from_wetmmr(team,wetmmr_for_testing,qv_wet,drymmr);
    
    //convert qv_wet to qv_dry
    PF::calculate_drymmr_from_wetmmr(team,qv_wet,qv_wet,qv_dry);

    // Convert drymmr computed above to wetmmr
    PF::calculate_wetmmr_from_drymmr(team,drymmr,qv_dry,wetmmr);

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
  auto T_from_dse_host      = cmvdc(T_from_dse);
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
    REQUIRE( Check::is_non_negative(T_from_dse_host(k),k) );
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

//===============================================================================
TEST_CASE("common_physics_functions_test", "[common_physics_functions_test]"){

  using scream::Real;
  using Device = scream::DefaultDevice;

  constexpr int num_runs = 5;

  auto engine = scream::setup_random_test();

// Run tests of vertically dimensioned-functions for both Real and Pack,
// and for (potentially) different pack sizes

  printf(" -> Number of randomized runs: %d\n\n", num_runs);

  printf(" -> Testing Real scalar type...");
  for (int irun=0; irun<num_runs; ++irun) {
    run<Real,Device>(engine);
    run_scalar_valued_fns<Device>(engine);
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
