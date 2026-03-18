#include "gw_functions.hpp"
#include "eamxx_gw_process_interface.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

#include <array>

namespace scream
{

/*------------------------------------------------------------------------------------------------*/
GWDrag::GWDrag(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params) {
  // Nothing to do here
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::create_requests() {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  constexpr int pack_size = Pack::n;

  m_grid = m_grids_manager->get_grid("physics");

  const auto& grid_name = m_grid->name();
  const auto layout     = m_grid->get_3d_scalar_layout(true);

  // retrieve local grid parameters
  m_ncol = m_grid->get_num_local_dofs();
  m_nlev = m_grid->get_num_vertical_levels();
  m_npgw = m_params.get<int>("pgwv") * 2 + 1;

  m_lat = m_grid->get_geometry_data("lat");

  const auto nondim = Units::nondimensional();
  const auto m2     = pow(m,2);
  const auto s2     = pow(s,2);
  const auto K2     = pow(K,2);

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();        // 2D variables
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);    // 3D variables at mid-levels
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);   // 3D variables at interfaces
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);  // horiz_wind field

  // Input variables
  add_field<Required>("p_mid",                scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("p_int",                scalar3d_int, Pa,     grid_name, pack_size);
  add_field<Required>("pseudo_density",       scalar3d_mid, Pa,     grid_name, pack_size);
  add_field<Required>("phis",                 scalar2d    , m2/s2,  grid_name);
  add_field<Required>("omega",                scalar3d_mid, Pa/s,   grid_name, pack_size);
  add_field<Required>("landfrac",             scalar2d    , nondim, grid_name);
  add_field<Required>("sgh",                  scalar2d    , nondim, grid_name);

  // Input/Output variables
  add_field <Updated>("T_mid",                scalar3d_mid, K,      grid_name, pack_size);
  add_tracer<Updated>("qv",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qc",                   m_grid,       kg/kg,             pack_size);
  add_tracer<Updated>("qi",                   m_grid,       kg/kg,             pack_size);
  add_field <Updated>("horiz_winds",          vector3d_mid, m/s,    grid_name, pack_size);

  // // Output variables
  // add_field <Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");
  // add_field <Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name, "ACCUMULATED");

  // Diagnostic Outputs
  add_field<Computed>("gw_activity",          scalar2d,     nondim, grid_name);
  add_field<Computed>("gw_T_mid_tend",        scalar3d_mid, K/s,    grid_name, pack_size);
  add_field<Computed>("gw_qv_tend",           scalar3d_mid, kg/kg/s,grid_name, pack_size);
  add_field<Computed>("gw_u_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
  add_field<Computed>("gw_v_tend",            scalar3d_mid, m/s/s,  grid_name, pack_size);
}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::initialize_impl (const RunType) {

  // defaults for gw_common_init()
  bool do_molec_diff_default = false; // Flag for molecular diffusion
  int nbot_molec_default = 0;         // bottom level for molecular diffusion
  int ktop_default = 0;               // Top level for gravity waves.
  Real kwv_default = 6.28e-5;         // Effective horizontal wave number (100 km wavelength)

  // calculate interface reference pressures
  const auto hyai = m_grid->get_geometry_data("hyai").get_view<const Real*>();
  const auto hybi = m_grid->get_geometry_data("hybi").get_view<const Real*>();
  Kokkos::View<Real*, Kokkos::HostSpace> pref_int("pref_int", hyai.size());
  Kokkos::parallel_for("calculate_pref_int", 
    Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, hyai.size()), 
    KOKKOS_LAMBDA (const int k) {
      pref_int(k) = PC::P0.value * hyai(k) + PC::P0.value * hybi(k);
  });
  Kokkos::fence();

  GWF::gw_common_init( m_params, 
                       m_nlev,
                       pref_int,
                       do_molec_diff_default,
                       nbot_molec_default,
                       ktop_default,
                       kwv_default);

  std::string gw_drag_file = m_params.get<std::string>("gw_drag_file");
  scorpio::register_file(gw_drag_file,scorpio::FileMode::Read);
  const int PS_dim_size = scorpio::get_dimlen(gw_drag_file, "PS"); // Phase Speed [m/s]
  const int MW_dim_size = scorpio::get_dimlen(gw_drag_file, "MW"); // Mean Wind in Heating [m/s]
  const int HD_dim_size = scorpio::get_dimlen(gw_drag_file, "HD"); // Heating Depth [km]
  Kokkos::View<Real***, Kokkos::HostSpace> mfcc_host("mfcc_host", HD_dim_size, MW_dim_size, PS_dim_size);
  scorpio::read_var(gw_drag_file,"mfcc",mfcc_host.data());
  scorpio::release_file(gw_drag_file);

  GWF::gw_convect_init( m_params, mfcc_host );

  GWF::gw_front_init( m_params, pref_int );

  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),       m_grid,100.0,400.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"), m_grid,-200.0, 200.0,false);

}

/*------------------------------------------------------------------------------------------------*/
void GWDrag::run_impl (const double dt) {
  using PC = scream::physics::Constants<Real>;
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  const auto team_policy = TPF::get_default_team_policy(m_ncol, nlev_mid_packs);
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(m_ncol, nlev_mid_packs);
  // Use one workspace with the biggest size or use two, one for pver, one for pver*2*pgwv?
  WSM wsm( (m_nlev+1)*m_npgw, 9, team_policy);
  //----------------------------------------------------------------------------
  // get fields

  // variables not updated by GWD
  const auto& phis        = get_field_in("phis")          .get_view<const Real*>();

  const auto& p_mid       = get_field_in("p_mid")         .get_view<const Pack**>();
  const auto& p_int       = get_field_in("p_int")         .get_view<const Pack**>();
  const auto& p_del       = get_field_in("pseudo_density").get_view<const Pack**>();
  const auto& omega       = get_field_in("omega")         .get_view<const Pack**>();
  const auto& landfrac    = get_field_in("landfrac")      .get_view<const Real*>();
  const auto& sgh         = get_field_in("sgh")           .get_view<const Real*>();

  // variables updated by GWD
  const auto& T_mid       = get_field_out("T_mid")        .get_view<Pack**>();
  const auto& qv          = get_field_out("qv")           .get_view<Pack**>();
  const auto& qc          = get_field_out("qc")           .get_view<Pack**>();
  const auto& qi          = get_field_out("qi")           .get_view<Pack**>();
  const auto& hwinds_fld  = get_field_out("horiz_winds");
  const auto& uwind       = hwinds_fld.get_component(0)   .get_view<Pack**>();
  const auto& vwind       = hwinds_fld.get_component(1)   .get_view<Pack**>();

  auto m_lat_v = m_lat.get_view<const Real*>();
  //----------------------------------------------------------------------------
  // create local temporaries to avoid "Implicit capture" warning
  const auto loc_phis  = phis;
  const auto loc_p_mid = p_mid;
  const auto loc_p_int = p_int;
  const auto loc_p_del = p_del;

  auto loc_T_mid       = T_mid;
  auto loc_qv          = qv;
  auto loc_qc          = qc;
  auto loc_qi          = qi;
  auto loc_uwind       = uwind;
  auto loc_vwind       = vwind;
  auto loc_landfrac    = landfrac;

  //----------------------------------------------------------------------------
  // NaN check after gw_oro_src
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const Int i = team.league_rank();
    const auto uwind_i      = ekat::scalarize(ekat::subview(loc_uwind, i));
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      int nan_count = 0, inf_count = 0;
      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(uwind_i(k))) nan_count++; if (Kokkos::isinf(uwind_i(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) { printf("[run_impl  IN ] col %d: uwind  NaN=%d Inf=%d\n", i, nan_count, inf_count); }
    });
  });
  //----------------------------------------------------------------------------
  // // calculate altitude on interfaces (z_int) and mid-points (z_mid)

  // // create temporaries to avoid "Implicit capture" warning
  // const auto loc_p_mid = p_mid;
  // const auto loc_p_del = p_del;
  // const auto loc_T_mid = T_mid;
  // const auto loc_qv    = qv;
  // auto loc_z_mid = m_buffer.z_mid;
  // auto loc_z_del = m_buffer.z_del;
  // auto loc_z_int = m_buffer.z_int;
  // auto loc_nlev = m_nlev;

  // Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
  //   const int i = team.league_rank();
  //   const auto p_mid_i = ekat::subview(loc_p_mid, i);
  //   const auto p_del_i = ekat::subview(loc_p_del, i);
  //   const auto T_mid_i = ekat::subview(loc_T_mid, i);
  //   const auto qv_i    = ekat::subview(loc_qv,    i);
  //   auto z_mid_i = ekat::subview(loc_z_mid, i);
  //   auto z_del_i = ekat::subview(loc_z_del, i);
  //   auto z_int_i = ekat::subview(loc_z_int, i);
  //   auto z_surf = 0.0; // z_mid & z_int are altitude above the surface
  //   PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
  //   team.team_barrier();
  //   PF::calculate_z_int(team, loc_nlev, z_del_i, z_surf, z_int_i);
  //   team.team_barrier();
  //   PF::calculate_z_mid(team, loc_nlev, z_int_i, z_mid_i);
  //   team.team_barrier();
  // });
  //----------------------------------------------------------------------------

  auto loc_z_mid       = m_buffer.z_mid;
  auto loc_z_del       = m_buffer.z_del;
  auto loc_z_int       = m_buffer.z_int;
  auto loc_p_del_rcp   = m_buffer.p_del_rcp;
  auto loc_p_int_log   = m_buffer.p_int_log;
  auto loc_T_int       = m_buffer.T_int;
  auto loc_N_mid       = m_buffer.N_mid;
  auto loc_N_int       = m_buffer.N_int;
  auto loc_rho_int     = m_buffer.rho_int;
  auto loc_q_combined  = m_buffer.q_combined;
  auto loc_tau         = m_buffer.tau;
  auto loc_ubm         = m_buffer.ubm;
  auto loc_ubi         = m_buffer.ubi;
  auto loc_c           = m_buffer.c;
  auto loc_kvtt        = m_buffer.kvtt;
  auto loc_dse         = m_buffer.dse;
  auto loc_utgw        = m_buffer.utgw;
  auto loc_vtgw        = m_buffer.vtgw;
  auto loc_ttgw        = m_buffer.ttgw;
  auto loc_qtgw        = m_buffer.qtgw;
  auto loc_gw_tend_u   = m_buffer.gw_tend_u;
  auto loc_gw_tend_v   = m_buffer.gw_tend_v;
  auto loc_gw_tend_t   = m_buffer.gw_tend_t;
  auto loc_gw_tend_q   = m_buffer.gw_tend_q;
  auto loc_taucd       = m_buffer.taucd;
  auto loc_egwdffi     = m_buffer.egwdffi;
  auto loc_gwut        = m_buffer.gwut;
  auto loc_dttdf       = m_buffer.dttdf;
  auto loc_dttke       = m_buffer.dttke;
  //----------------------------------------------------------------------------
  // populate q_combined
  Kokkos::parallel_for(KT::RangePolicy(0, m_ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid_packs;
    const int k = idx%nlev_mid_packs;
    loc_q_combined(i, k, 0) = loc_qv(i,k);
    loc_q_combined(i, k, 1) = loc_qc(i,k);
    loc_q_combined(i, k, 2) = loc_qi(i,k);
  });
  //----------------------------------------------------------------------------
  // calculate altitude on interfaces (z_int) and mid-points (z_mid)
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto p_mid_i = ekat::subview(loc_p_mid, i);
    const auto p_del_i = ekat::subview(loc_p_del, i);
    const auto T_mid_i = ekat::subview(loc_T_mid, i);
    const auto qv_i = ekat::subview(loc_qv,    i);
    auto z_mid_i = ekat::subview(loc_z_mid, i);
    auto z_del_i = ekat::subview(loc_z_del, i);
    auto z_int_i = ekat::subview(loc_z_int, i);
    auto z_surf = 0.0; // z_mid & z_int are altitude above the surface
    PF::calculate_dz(team, p_del_i, p_mid_i, T_mid_i, qv_i, z_del_i);
    team.team_barrier();
    PF::calculate_z_int(team, m_nlev, z_del_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, m_nlev, z_int_i, z_mid_i);
    team.team_barrier();
  });
  //----------------------------------------------------------------------------
  // miscellaneous calculations
  Kokkos::parallel_for(KT::RangePolicy(0, m_ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid_packs;
    const int k = idx%nlev_mid_packs;
    loc_dse(i,k) = PF::calculate_dse(loc_T_mid(i,k),loc_z_mid(i,k),loc_phis(i));
    loc_p_del_rcp(i,k) = 1.0/loc_p_del(i,k);
  });
  Kokkos::parallel_for(KT::RangePolicy(0, m_ncol*nlev_int_packs), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_int_packs;
    const int k = idx%nlev_int_packs;
    loc_p_int_log(i,k) = ekat::log(loc_p_int(i,k));
  });
  //----------------------------------------------------------------------------
  // initialize output tendencies
  Kokkos::deep_copy(loc_gw_tend_u,0.0);
  Kokkos::deep_copy(loc_gw_tend_v,0.0);
  Kokkos::deep_copy(loc_gw_tend_t,0.0);
  Kokkos::deep_copy(loc_gw_tend_q,0.0);
  // initialize intermediate per-source tendency and diffusivity buffers;
  // gw_drag_prof may not write all levels (e.g., above the source level), so
  // uninitialized memory would corrupt the accumulation into gw_tend_*
  Kokkos::deep_copy(loc_tau,  0.0);
  Kokkos::deep_copy(loc_utgw, 0.0);
  Kokkos::deep_copy(loc_vtgw, 0.0);
  Kokkos::deep_copy(loc_ttgw, 0.0);
  // Kokkos::deep_copy(loc_gwut, 0.0);
  Kokkos::deep_copy(loc_qtgw, 0.0);
  Kokkos::deep_copy(loc_kvtt, 0.0);

  //----------------------------------------------------------------------------
  // Compute profiles of background state
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const Int i = team.league_rank();
    // Get single-column subviews of all inputs
    // Note: scalarize includes pack padding, so we need to slice to actual sizes
    const auto T_mid_s   = ekat::scalarize(ekat::subview(loc_T_mid,  i));
    const auto p_mid_s   = ekat::scalarize(ekat::subview(loc_p_mid,  i));
    const auto p_int_s   = ekat::scalarize(ekat::subview(loc_p_int,  i));
    const auto T_int_s   = ekat::scalarize(ekat::subview(loc_T_int,  i));
    const auto N_mid_s   = ekat::scalarize(ekat::subview(loc_N_mid,  i));
    const auto N_int_s   = ekat::scalarize(ekat::subview(loc_N_int,  i));
    const auto rho_int_s = ekat::scalarize(ekat::subview(loc_rho_int,i));

    // Slice to exclude padding: mid levels = [0, m_nlev), int levels = [0, m_nlev+1)
    const auto T_mid_i   = ekat::subview(T_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto p_mid_i   = ekat::subview(p_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto p_int_i   = ekat::subview(p_int_s,   Kokkos::pair<int, int>{0, m_nlev+1});
    // const auto T_int_i   = ekat::subview(T_int_s,   Kokkos::pair<int, int>{0, m_nlev}); // note: interface data w/ m_nlev dimension
    const auto T_int_i   = ekat::subview(T_int_s,   Kokkos::pair<int, int>{0, m_nlev+1});
    const auto N_mid_i   = ekat::subview(N_mid_s,   Kokkos::pair<int, int>{0, m_nlev});
    const auto N_int_i   = ekat::subview(N_int_s,   Kokkos::pair<int, int>{0, m_nlev+1});
    const auto rho_int_i = ekat::subview(rho_int_s, Kokkos::pair<int, int>{0, m_nlev+1});

    GWF::gw_prof(team, m_nlev, PC::Cpair.value,
                 T_mid_i, p_mid_i, p_int_i,
                 rho_int_i, T_int_i, N_mid_i, N_int_i);

    // // NaN check after gw_prof
    // Kokkos::single(Kokkos::PerTeam(team), [&]() {
    //   int nan_count = 0, inf_count = 0;

    //   // --- inputs ---
    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_i(k))) nan_count++; if (Kokkos::isinf(T_mid_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) { printf("[gw_prof  IN ] col %d: T_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(p_mid_i(k))) nan_count++; if (Kokkos::isinf(p_mid_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) { printf("[gw_prof  IN ] col %d: p_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(p_int_i(k))) nan_count++; if (Kokkos::isinf(p_int_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) { printf("[gw_prof  IN ] col %d: p_int  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

    //   // --- outputs ---
    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(T_int_i(k))) nan_count++; if (Kokkos::isinf(T_int_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) {
    //     printf("[gw_prof  OUT] col %d: T_int  NaN=%d Inf=%d\n", i, nan_count, inf_count);
    //     // for (int k = 0; k <= m_nlev; k++) {
    //     //   if (Kokkos::isnan(T_int_i(k))) { printf("  T_int[%d]=NaN\n", k); }
    //     // }
    //   }
    //   // Print boundary values to detect OOB corruption
    //   printf("[gw_prof  DBG] col %d: T_int[%d]=%e T_int[%d]=%e\n",
    //          i, m_nlev-1, (double)T_int_i(m_nlev-1), m_nlev, (double)T_int_i(m_nlev));

    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(rho_int_i(k))) nan_count++; if (Kokkos::isinf(rho_int_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) {
    //     printf("[gw_prof  OUT] col %d: rho_int NaN=%d Inf=%d\n", i, nan_count, inf_count);
    //     // for (int k = 0; k <= m_nlev; k++) {
    //     //   if (Kokkos::isnan(rho_int_i(k))) { printf("  rho_int[%d]=NaN\n", k); }
    //     // }
    //   }

    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(N_mid_i(k))) nan_count++; if (Kokkos::isinf(N_mid_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) {
    //     printf("[gw_prof  OUT] col %d: N_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count);
    //     // for (int k = 0; k < m_nlev; k++) {
    //     //   if (Kokkos::isnan(N_mid_i(k))) { printf("  N_mid[%d]=NaN\n", k); }
    //     // }
    //   }

    //   nan_count = 0; inf_count = 0;
    //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(N_int_i(k))) nan_count++; if (Kokkos::isinf(N_int_i(k))) inf_count++; }
    //   if (nan_count > 0 || inf_count > 0) {
    //     printf("[gw_prof  OUT] col %d: N_int  NaN=%d Inf=%d\n", i, nan_count, inf_count);
    //     // for (int k = 0; k <= m_nlev; k++) {
    //     //   if (Kokkos::isnan(N_int_i(k))) { printf("  N_int[%d]=NaN\n", k); }
    //     // }
    //   }
    // });

  });

  


  //----------------------------------------------------------------------------
  // Calculate local molecular diffusivity
  // NOTE - if we need moelcular diffusion then this is where we would calculate
  // the local diffusivity. However, the molecular viscosity coefficient increases
  // exponentially with altitude, so this process become increasingly important
  // in thermal conduction, species transport, and momentum damping at higher
  // altitudes. The critical altitude is the homopause/turbopause, where the molecular
  // diffusion/viscosity becomes comparable with turbulence diffusion, and is
  // dominant above. So for model tops above 100 km molecular diffusion/viscosity
  // is needed (i.e. WACCM), but currently this is not a priority for EAMxx.
  // if (do_molec_diff) { ??? }
  //----------------------------------------------------------------------------
  // Convective gravity waves (Beres scheme)
  if (GWF::s_common_init.use_gw_convect) {

    // // Determine wave sources
    // GWF::gw_beres_src();

    // // Solve for the drag profile
    // GWF::gw_drag_prof();

    // add the diffusion coefficients
    // do k = 0, pver
    //   egwdffi_tot(:,k) = egwdffi_tot(:,k) + egwdffi(:,k)
    // end do

    // Store constituents tendencies
    // do m=1, pcnst
    //    do k = 1, pver
    //       ptend%q(:ncol,k,m) = qtgw(:,k,m)
    //    end do
    // end do

    // ! add the momentum tendencies to the output tendency arrays
    // do k = 1, pver
    //    ptend%u(:ncol,k) = utgw(:,k)
    //    ptend%v(:ncol,k) = vtgw(:,k)
    //    ptend%s(:ncol,k) = ttgw(:,k)
    // end do

    // // Momentum & energy conservation
    // GWF::momentum_energy_conservation();

  }

  // Frontally generated gravity waves
  if (GWF::s_common_init.use_gw_frontal) {
    // GWF::gw_cm_src();
    // GWF::gw_drag_prof();
    // GWF::momentum_energy_conservation();
  }

  // Orographic stationary gravity waves
  if (GWF::s_common_init.use_gw_orographic) {

    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();

      const Real landfrac_i = loc_landfrac(i);

      // Get single-column subviews of all inputs
      const auto uwind_i      = ekat::scalarize(ekat::subview(loc_uwind, i));
      const auto vwind_i      = ekat::scalarize(ekat::subview(loc_vwind, i));
      const auto T_mid_i      = ekat::scalarize(ekat::subview(loc_T_mid, i));
      const auto T_int_i      = ekat::scalarize(ekat::subview(loc_T_int, i));
      const auto p_mid_i      = ekat::scalarize(ekat::subview(loc_p_mid, i));
      const auto p_int_i      = ekat::scalarize(ekat::subview(loc_p_int, i));
      const auto p_del_i      = ekat::scalarize(ekat::subview(loc_p_del, i));
      const auto p_del_rcp_i  = ekat::scalarize(ekat::subview(loc_p_del_rcp, i));
      const auto p_int_log_i  = ekat::scalarize(ekat::subview(loc_p_int_log, i));
      const auto z_mid_i      = ekat::scalarize(ekat::subview(loc_z_mid, i));
      const auto N_mid_i      = ekat::scalarize(ekat::subview(loc_N_mid, i));
      const auto N_int_i      = ekat::scalarize(ekat::subview(loc_N_int, i));
      const auto rho_int_i    = ekat::scalarize(ekat::subview(loc_rho_int,i));
      const auto tau_i        = ekat::scalarize(ekat::subview(loc_tau, i));
      const auto ubm_i        = ekat::scalarize(ekat::subview(loc_ubm, i));
      const auto ubi_i        = ekat::scalarize(ekat::subview(loc_ubi, i));
      const auto c_i          = ekat::scalarize(ekat::subview(loc_c, i));
      const auto kvtt_i       = ekat::scalarize(ekat::subview(loc_kvtt, i));
      const auto dse_i        = ekat::scalarize(ekat::subview(loc_dse, i));
      const auto utgw_i       = ekat::scalarize(ekat::subview(loc_utgw, i));
      const auto vtgw_i       = ekat::scalarize(ekat::subview(loc_vtgw, i));
      const auto ttgw_i       = ekat::scalarize(ekat::subview(loc_ttgw, i));
      const auto qtgw_i       = ekat::scalarize(ekat::subview(loc_qtgw, i));
      const auto gw_tend_u_i  = ekat::scalarize(ekat::subview(loc_gw_tend_u, i));
      const auto gw_tend_v_i  = ekat::scalarize(ekat::subview(loc_gw_tend_v, i));
      const auto gw_tend_t_i  = ekat::scalarize(ekat::subview(loc_gw_tend_t, i));
      const auto gw_tend_q_i  = ekat::scalarize(ekat::subview(loc_gw_tend_q, i));
      const auto taucd_i      = ekat::scalarize(ekat::subview(loc_taucd, i));
      const auto egwdffi_i    = ekat::scalarize(ekat::subview(loc_egwdffi, i));
      const auto gwut_i       = ekat::scalarize(ekat::subview(loc_gwut, i));
      const auto dttdf_i      = ekat::scalarize(ekat::subview(loc_dttdf, i));
      const auto dttke_i      = ekat::scalarize(ekat::subview(loc_dttke, i));

      const auto q_comb_s = ekat::scalarize(ekat::subview(loc_q_combined, i));
      const auto q_2d = Kokkos::subview(q_comb_s,Kokkos::pair<int,int>{0, m_nlev},Kokkos::ALL);

      Int src_lev; // level index of gravity wave source
      Int tnd_lev; // lowest level index where tendencies are allowed
      Real xv;      // zonal unit vector of source wind
      Real yv;      // meridional unit vector of source wind

      // Determine the orographic wave source
      GWF::gw_oro_src(team, GWF::s_common_init, m_nlev, GWF::s_common_init.pgwv,
                      uwind_i, vwind_i, T_mid_i, sgh(i),
                      p_mid_i, p_int_i, p_del_i, z_mid_i, N_mid_i,
                      src_lev, tnd_lev,
                      tau_i, ubm_i, ubi_i, xv, yv, c_i );

      // Kokkos::fence();

      // // NaN check after gw_oro_src
      // Kokkos::single(Kokkos::PerTeam(team), [&]() {
      //   int nan_count = 0, inf_count = 0;
      //   printf("[oro_src  DBG] col %d: src_lev=%d tnd_lev=%d xv=%e yv=%e\n",
      //          i, (int)src_lev, (int)tnd_lev, (double)xv, (double)yv);

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(uwind_i(k))) nan_count++; if (Kokkos::isinf(uwind_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: uwind  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(vwind_i(k))) nan_count++; if (Kokkos::isinf(vwind_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: vwind  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_i(k))) nan_count++; if (Kokkos::isinf(T_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: T_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(N_mid_i(k))) nan_count++; if (Kokkos::isinf(N_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: N_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   // tau: shape (m_npgw, lev) - interface levels
      //   nan_count = 0; inf_count = 0;
      //   for (int pg = 0; pg < m_npgw; pg++) {
      //     for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(tau_i(pg,k))) nan_count++; if (Kokkos::isinf(tau_i(pg,k))) inf_count++; }
      //   }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[oro_src  OUT] col %d: tau    NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     for (int pg = 0; pg < m_npgw; pg++) {
      //       for (int k = 0; k <= m_nlev; k++) {
      //         if (Kokkos::isnan(tau_i(pg,k))) { printf("  tau[pg=%d,lev=%d]=NaN\n", pg, k); }
      //       }
      //     }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(ubm_i(k))) nan_count++; if (Kokkos::isinf(ubm_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[oro_src  OUT] col %d: ubm    NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(ubm_i(k))) { printf("  ubm[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(ubi_i(k))) nan_count++; if (Kokkos::isinf(ubi_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[oro_src  OUT] col %d: ubi    NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k <= m_nlev; k++) {
      //     //   if (Kokkos::isnan(ubi_i(k))) { printf("  ubi[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(T_int_i(k))) nan_count++; if (Kokkos::isinf(T_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[oro_src  IN ] col %d: T_int  NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k <= m_nlev; k++) {
      //     //   if (Kokkos::isnan(T_int_i(k))) { printf("  T_int[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(p_int_i(k))) nan_count++; if (Kokkos::isinf(p_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: p_int  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(p_del_i(k))) nan_count++; if (Kokkos::isinf(p_del_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: p_del  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(N_int_i(k))) nan_count++; if (Kokkos::isinf(N_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: N_int  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(rho_int_i(k))) nan_count++; if (Kokkos::isinf(rho_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[oro_src  IN ] col %d: rho_int NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      // });

      // Solve for the drag profile with orographic sources
      Int max_lev = tnd_lev; // ???

      // // NaN check: gw_drag_prof inputs
      // Kokkos::single(Kokkos::PerTeam(team), [&]() {
      //   int nan_count = 0, inf_count = 0;
      //   printf("[drag_prof IN ] col %d: src_lev=%d max_lev=%d tnd_lev=%d\n",
      //          i, (int)src_lev, (int)max_lev, (int)tnd_lev);

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_i(k))) nan_count++; if (Kokkos::isinf(T_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: T_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(T_int_i(k))) nan_count++; if (Kokkos::isinf(T_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof IN ] col %d: T_int  NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k <= m_nlev; k++) {
      //     //   if (Kokkos::isnan(T_int_i(k))) { printf("  T_int[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(N_mid_i(k))) nan_count++; if (Kokkos::isinf(N_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof IN ] col %d: N_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(N_mid_i(k))) { printf("  N_mid[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(N_int_i(k))) nan_count++; if (Kokkos::isinf(N_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: N_int  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(rho_int_i(k))) nan_count++; if (Kokkos::isinf(rho_int_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: rho_int NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(ubi_i(k))) nan_count++; if (Kokkos::isinf(ubi_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: ubi    NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(ubm_i(k))) nan_count++; if (Kokkos::isinf(ubm_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: ubm    NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   // tau: shape (m_npgw, lev)
      //   nan_count = 0; inf_count = 0;
      //   for (int pg = 0; pg < m_npgw; pg++) {
      //     for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(tau_i(pg,k))) nan_count++; if (Kokkos::isinf(tau_i(pg,k))) inf_count++; }
      //   }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: tau    NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(dse_i(k))) nan_count++; if (Kokkos::isinf(dse_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof IN ] col %d: dse    NaN=%d Inf=%d\n", i, nan_count, inf_count); }
      // });

      // std::cout << "WHDEBUG 1" << std::endl;
      // std::cout.flush();
      // Kokkos::fence();
      // std::cout << "WHDEBUG 2" << std::endl;
      // std::cout.flush();
      GWF::gw_drag_prof(team, wsm.get_workspace(team), GWF::s_common_init,
                        m_nlev, GWF::s_common_init.pgwv, src_lev, max_lev, tnd_lev,
                        GWF::s_common_init.do_taper, dt, m_lat_v(i),
                        T_mid_i, T_int_i, p_mid_i, p_int_i,
                        p_del_i, p_del_rcp_i, p_int_log_i, rho_int_i,
                        N_mid_i, N_int_i, ubm_i, ubi_i, xv, yv,
                        GWF::s_common_init.gw_orographic_eff,
                        c_i, kvtt_i, q_2d, dse_i, tau_i,
                        utgw_i, vtgw_i, ttgw_i, qtgw_i,
                        taucd_i, egwdffi_i, gwut_i, dttdf_i, dttke_i);
      // Kokkos::fence();
      // std::cout << "WHDEBUG 3" << std::endl;
      // std::cout.flush();

      // // NaN check: gw_drag_prof outputs
      // Kokkos::single(Kokkos::PerTeam(team), [&]() {
      //   int nan_count = 0, inf_count = 0;

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(utgw_i(k))) nan_count++; if (Kokkos::isinf(utgw_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof OUT] col %d: utgw   NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(utgw_i(k))) { printf("  utgw[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(vtgw_i(k))) nan_count++; if (Kokkos::isinf(vtgw_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof OUT] col %d: vtgw   NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(vtgw_i(k))) { printf("  vtgw[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(ttgw_i(k))) nan_count++; if (Kokkos::isinf(ttgw_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof OUT] col %d: ttgw   NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(ttgw_i(k))) { printf("  ttgw[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(dttdf_i(k))) nan_count++; if (Kokkos::isinf(dttdf_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof OUT] col %d: dttdf  NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k < m_nlev; k++) {
      //     //   if (Kokkos::isnan(dttdf_i(k))) { printf("  dttdf[%d]=NaN\n", k); }
      //     // }
      //   }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(dttke_i(k))) nan_count++; if (Kokkos::isinf(dttke_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: dttke  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(kvtt_i(k))) nan_count++; if (Kokkos::isinf(kvtt_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: kvtt   NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(egwdffi_i(k))) nan_count++; if (Kokkos::isinf(egwdffi_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: egwdffi NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   // taucd: shape (lev, 4) - interface levels
      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k <= m_nlev; k++) {
      //     for (int m = 0; m < 4; m++) { if (Kokkos::isnan(taucd_i(k,m))) nan_count++; if (Kokkos::isinf(taucd_i(k,m))) inf_count++; }
      //   }
      //   if (nan_count > 0 || inf_count > 0) {
      //     printf("[drag_prof OUT] col %d: taucd  NaN=%d Inf=%d\n", i, nan_count, inf_count);
      //     // for (int k = 0; k <= m_nlev; k++) {
      //     //   for (int m = 0; m < 4; m++) {
      //     //     if (Kokkos::isnan(taucd_i(k,m))) { printf("  taucd[lev=%d,m=%d]=NaN\n", k, m); }
      //     //   }
      //     // }
      //   }

      //   // tau after drag_prof
      //   nan_count = 0; inf_count = 0;
      //   for (int pg = 0; pg < m_npgw; pg++) {
      //     for (int k = 0; k <= m_nlev; k++) { if (Kokkos::isnan(tau_i(pg,k))) nan_count++; if (Kokkos::isinf(tau_i(pg,k))) inf_count++; }
      //   }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: tau    NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   // N_mid after drag_prof (should be unchanged but check for corruption)
      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(N_mid_i(k))) nan_count++; if (Kokkos::isinf(N_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: N_mid  NaN=%d Inf=%d (CORRUPTED)\n", i, nan_count, inf_count); }

      //   // T_mid after drag_prof (should be unchanged)
      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_i(k))) nan_count++; if (Kokkos::isinf(T_mid_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[drag_prof OUT] col %d: T_mid  NaN=%d Inf=%d (CORRUPTED)\n", i, nan_count, inf_count); }

      // });

      // add tendencies to aggregate output tendencies
      team.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, m_nlev), [&] (const int k) {
        gw_tend_u_i(k)   += utgw_i(k)   * landfrac_i;
        gw_tend_v_i(k)   += vtgw_i(k)   * landfrac_i;
        gw_tend_t_i(k)   += ttgw_i(k)   * landfrac_i;
        gw_tend_q_i(k,0) += qtgw_i(k,0) * landfrac_i;
        gw_tend_q_i(k,1) += qtgw_i(k,1) * landfrac_i;
        gw_tend_q_i(k,2) += qtgw_i(k,2) * landfrac_i;
      });

      // // NaN check after tendency accumulation
      // team.team_barrier();
      // Kokkos::single(Kokkos::PerTeam(team), [&]() {
      //   int nan_count = 0, inf_count = 0;

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(gw_tend_u_i(k))) nan_count++; if (Kokkos::isinf(gw_tend_u_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[tend_accum  ] col %d: gw_tend_u NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(gw_tend_v_i(k))) nan_count++; if (Kokkos::isinf(gw_tend_v_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[tend_accum  ] col %d: gw_tend_v NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      //   nan_count = 0; inf_count = 0;
      //   for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(gw_tend_t_i(k))) nan_count++; if (Kokkos::isinf(gw_tend_t_i(k))) inf_count++; }
      //   if (nan_count > 0 || inf_count > 0) { printf("[tend_accum  ] col %d: gw_tend_t NaN=%d Inf=%d\n", i, nan_count, inf_count); }
      // });

    });

    //----------------------------------------------------------------------------
    // GW energy fixer
    Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
      const Int i = team.league_rank();

      const auto p_del_i = ekat::scalarize(ekat::subview(loc_p_del,      i));
      const auto p_int_i = ekat::scalarize(ekat::subview(loc_p_int,      i));
      const auto uwind_i = ekat::scalarize(ekat::subview(loc_uwind,      i));
      const auto vwind_i = ekat::scalarize(ekat::subview(loc_vwind,      i));
      const auto gw_tend_u_i = ekat::scalarize(ekat::subview(loc_gw_tend_u, i));
      const auto gw_tend_v_i = ekat::scalarize(ekat::subview(loc_gw_tend_v, i));
      const auto gw_tend_t_i = ekat::scalarize(ekat::subview(loc_gw_tend_t, i));

      Real dE = 0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, 0, m_nlev), [&] (const int k, Real& lsum) {
        lsum -= p_del_i(k) * ( gw_tend_u_i(k) * (uwind_i(k) + gw_tend_u_i(k) * GWF::GWC::half * dt)
                              +gw_tend_v_i(k) * (vwind_i(k) + gw_tend_v_i(k) * GWF::GWC::half * dt)
                              +gw_tend_t_i(k));
      }, Kokkos::Sum<Real>(dE));

      dE /= (p_int_i(m_nlev) - p_int_i(0));

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, m_nlev), [&] (const int k) {
        gw_tend_t_i(k) += dE;
      });
    });

  } // use_gw_orographic

  //----------------------------------------------------------------------------
  // Convert the tendencies for the dry constituents to dry air basis

  // The fortran version used the conversion below - but since the constituent
  // tendencies only use the first 3 tracers (vapor, cloud liquid, cloud ice)
  // we can skip this step, and just note in the usage changes in the future

  // do m = 1, pcnst
  //    if (cnst_type(m).eq.'dry') then
  //       do k = 1, pver
  //          do i = 1, ncol
  //             ptend%q(i,k,m) = ptend%q(i,k,m)*state1%pdel(i,k)/state1%pdeldry(i,k)
  //          end do
  //       end do
  //    end if
  // end do

  //----------------------------------------------------------------------------
  // NaN check before applying tendencies to prognostic fields
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const Int i = team.league_rank();
    const auto uwind_s = ekat::scalarize(ekat::subview(loc_uwind, i));
    const auto vwind_s = ekat::scalarize(ekat::subview(loc_vwind, i));
    const auto T_mid_s = ekat::scalarize(ekat::subview(loc_T_mid, i));
    const auto tend_u_s = ekat::scalarize(ekat::subview(loc_gw_tend_u, i));
    const auto tend_v_s = ekat::scalarize(ekat::subview(loc_gw_tend_v, i));
    const auto tend_t_s = ekat::scalarize(ekat::subview(loc_gw_tend_t, i));
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      int nan_count = 0, inf_count = 0;

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(uwind_s(k))) nan_count++; if (Kokkos::isinf(uwind_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) {
        printf("[pre_update  ] col %d: uwind  NaN=%d Inf=%d\n", i, nan_count, inf_count);
        // for (int k = 0; k < m_nlev; k++) {
        //   if (Kokkos::isnan(uwind_s(k))) { printf("  uwind[%d]=NaN\n", k); }
        // }
      }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(vwind_s(k))) nan_count++; if (Kokkos::isinf(vwind_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) {
        printf("[pre_update  ] col %d: vwind  NaN=%d Inf=%d\n", i, nan_count, inf_count);
        // for (int k = 0; k < m_nlev; k++) {
        //   if (Kokkos::isnan(vwind_s(k))) { printf("  vwind[%d]=NaN\n", k); }
        // }
      }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_s(k))) nan_count++; if (Kokkos::isinf(T_mid_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) { printf("[pre_update  ] col %d: T_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(tend_u_s(k))) nan_count++; if (Kokkos::isinf(tend_u_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) { printf("[pre_update  ] col %d: tend_u NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(tend_v_s(k))) nan_count++; if (Kokkos::isinf(tend_v_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) { printf("[pre_update  ] col %d: tend_v NaN=%d Inf=%d\n", i, nan_count, inf_count); }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(tend_t_s(k))) nan_count++; if (Kokkos::isinf(tend_t_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) { printf("[pre_update  ] col %d: tend_t NaN=%d Inf=%d\n", i, nan_count, inf_count); }
    });
  });

  //----------------------------------------------------------------------------
  // update prognostic fields
  Kokkos::parallel_for(KT::RangePolicy(0, m_ncol*nlev_mid_packs), KOKKOS_LAMBDA (const int idx) {
    const int i = idx/nlev_mid_packs;
    const int k = idx%nlev_mid_packs;
    T_mid(i,k) += loc_gw_tend_t(i,k) / PC::Cpair.value * dt;
    uwind(i,k) += loc_gw_tend_u(i,k) * dt;
    vwind(i,k) += loc_gw_tend_v(i,k) * dt;
    qv(i,k) += loc_gw_tend_q(i,k,0) * dt;
    qc(i,k) += loc_gw_tend_q(i,k,1) * dt;
    qi(i,k) += loc_gw_tend_q(i,k,2) * dt;
  });

  //----------------------------------------------------------------------------
  // NaN check after applying tendencies to prognostic fields
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const KT::MemberType& team) {
    const Int i = team.league_rank();
    const auto uwind_s = ekat::scalarize(ekat::subview(uwind, i));
    const auto vwind_s = ekat::scalarize(ekat::subview(vwind, i));
    const auto T_mid_s = ekat::scalarize(ekat::subview(T_mid, i));
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      int nan_count = 0, inf_count = 0;

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(uwind_s(k))) nan_count++; if (Kokkos::isinf(uwind_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) {
        printf("[post_update ] col %d: uwind  NaN=%d Inf=%d\n", i, nan_count, inf_count);
        // for (int k = 0; k < m_nlev; k++) {
        //   if (Kokkos::isnan(uwind_s(k))) { printf("  uwind[%d]=NaN\n", k); }
        // }
      }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(vwind_s(k))) nan_count++; if (Kokkos::isinf(vwind_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) {
        printf("[post_update ] col %d: vwind  NaN=%d Inf=%d\n", i, nan_count, inf_count);
        // for (int k = 0; k < m_nlev; k++) {
        //   if (Kokkos::isnan(vwind_s(k))) { printf("  vwind[%d]=NaN\n", k); }
        // }
      }

      nan_count = 0; inf_count = 0;
      for (int k = 0; k < m_nlev; k++) { if (Kokkos::isnan(T_mid_s(k))) nan_count++; if (Kokkos::isinf(T_mid_s(k))) inf_count++; }
      if (nan_count > 0 || inf_count > 0) {
        printf("[post_update ] col %d: T_mid  NaN=%d Inf=%d\n", i, nan_count, inf_count);
        // for (int k = 0; k < m_nlev; k++) {
        //   if (Kokkos::isnan(T_mid_s(k))) { printf("  T_mid[%d]=NaN\n", k); }
        // }
      }
    });
  });

}
/*------------------------------------------------------------------------------------------------*/
size_t GWDrag::requested_buffer_size_in_bytes() const
{
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  constexpr int pcnst = Buffer::pcnst;
  size_t gw_buffer_size = 0;

  gw_buffer_size += Buffer::num_3d_int_views*m_ncol*m_npgw*nlev_int_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_3d_pcnst_views*m_ncol*nlev_mid_packs*pcnst*sizeof(Pack);
  gw_buffer_size += Buffer::num_3d_cd_int_views*m_ncol*nlev_int_packs*4*sizeof(Pack);
  gw_buffer_size += Buffer::num_3d_pgw_views*m_ncol*nlev_mid_packs*m_npgw*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_mid_views*m_ncol*nlev_mid_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_int_views*m_ncol*nlev_int_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_pgw_views*m_ncol*m_npgw*sizeof(Pack);

  return gw_buffer_size;
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  //----------------------------------------------------------------------------
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  constexpr int pcnst = Buffer::pcnst;
  //----------------------------------------------------------------------------
  uview_3d* buffer_3d_int_view_ptrs[Buffer::num_3d_int_views] = {
    &m_buffer.tau
  };
  for (int i=0; i<Buffer::num_3d_int_views; ++i) {
    *buffer_3d_int_view_ptrs[i] = uview_3d(mem, m_ncol, m_npgw, nlev_int_packs);
    mem += buffer_3d_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_3d* buffer_3d_pcnst_view_ptrs[Buffer::num_3d_pcnst_views] = {
    &m_buffer.q_combined,
    &m_buffer.qtgw,
    &m_buffer.gw_tend_q
  };
  for (int i=0; i<Buffer::num_3d_pcnst_views; ++i) {
    *buffer_3d_pcnst_view_ptrs[i] = uview_3d(mem, m_ncol, nlev_mid_packs, pcnst);
    mem += buffer_3d_pcnst_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_3d* buffer_3d_cd_int_view_ptrs[Buffer::num_3d_cd_int_views] = {
    &m_buffer.taucd
  };
  for (int i=0; i<Buffer::num_3d_cd_int_views; ++i) {
    *buffer_3d_cd_int_view_ptrs[i] = uview_3d(mem, m_ncol, nlev_int_packs, 4);
    mem += buffer_3d_cd_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_3d* buffer_3d_pgw_view_ptrs[Buffer::num_3d_pgw_views] = {
    &m_buffer.gwut
  };
  for (int i=0; i<Buffer::num_3d_pgw_views; ++i) {
    *buffer_3d_pgw_view_ptrs[i] = uview_3d(mem, m_ncol, nlev_mid_packs, m_npgw);
    mem += buffer_3d_pgw_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_mid_view_ptrs[Buffer::num_2d_mid_views] = {
    &m_buffer.z_del,
    &m_buffer.z_mid,
    &m_buffer.N_mid,
    &m_buffer.ubm,
    &m_buffer.dse,
    &m_buffer.utgw,
    &m_buffer.vtgw,
    &m_buffer.ttgw,
    &m_buffer.gw_tend_u,
    &m_buffer.gw_tend_v,
    &m_buffer.gw_tend_t,
    &m_buffer.dttdf,
    &m_buffer.dttke,
    &m_buffer.p_del_rcp
  };
  for (int i=0; i<Buffer::num_2d_mid_views; ++i) {
    *buffer_2d_mid_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_mid_packs);
    mem += buffer_2d_mid_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_int_view_ptrs[Buffer::num_2d_int_views] = {
    &m_buffer.z_int,
    &m_buffer.T_int,
    &m_buffer.N_int,
    &m_buffer.rho_int,
    &m_buffer.kvtt,
    &m_buffer.ubi,
    &m_buffer.egwdffi,
    &m_buffer.p_int_log
  };
  for (int i=0; i<Buffer::num_2d_int_views; ++i) {
    *buffer_2d_int_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_int_packs);
    mem += buffer_2d_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_2d_pgw_view_ptrs[Buffer::num_2d_pgw_views] = {
    &m_buffer.c
  };
  for (int i=0; i<Buffer::num_2d_pgw_views; ++i) {
    *buffer_2d_pgw_view_ptrs[i] = uview_2d(mem, m_ncol, m_npgw);
    mem += buffer_2d_pgw_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for GWDrag.");
}
/*------------------------------------------------------------------------------------------------*/
size_t GWDrag::requested_buffer_size_in_bytes() const
{
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  size_t gw_buffer_size = 0;

  gw_buffer_size += Buffer::num_2d_mid_views*m_ncol*nlev_mid_packs*sizeof(Pack);
  gw_buffer_size += Buffer::num_2d_int_views*m_ncol*nlev_int_packs*sizeof(Pack);

  return gw_buffer_size;
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::init_buffers(const ATMBufferManager &buffer_manager)
{
  auto buffer_chk = ( buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes() );
  EKAT_REQUIRE_MSG(buffer_chk,"Error! Buffers size not sufficient.\n");
  //----------------------------------------------------------------------------
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());
  const int nlev_mid_packs = ekat::npack<Pack>(m_nlev);
  const int nlev_int_packs = ekat::npack<Pack>(m_nlev+1);
  //----------------------------------------------------------------------------
  uview_2d* buffer_mid_view_ptrs[Buffer::num_2d_mid_views] = {
    &m_buffer.z_del,
    &m_buffer.z_mid
  };
  for (int i=0; i<Buffer::num_2d_mid_views; ++i) {
    *buffer_mid_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_mid_packs);
    mem += buffer_mid_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  uview_2d* buffer_int_view_ptrs[Buffer::num_2d_int_views] = {
    &m_buffer.z_int
  };
  for (int i=0; i<Buffer::num_2d_int_views; ++i) {
    *buffer_int_view_ptrs[i] = uview_2d(mem, m_ncol, nlev_int_packs);
    mem += buffer_int_view_ptrs[i]->size();
  }
  //----------------------------------------------------------------------------
  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for TurbulentMountainStress.");
}
/*------------------------------------------------------------------------------------------------*/
void GWDrag::finalize_impl ()
{
  GWF::gw_finalize();
}
/*------------------------------------------------------------------------------------------------*/
} // namespace scream
