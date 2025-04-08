#include "physics/iop_forcing/eamxx_iop_forcing_process_interface.hpp"

#include "share/field/field_utils.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

namespace scream
{
// =========================================================================================
void IOPForcing::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();

  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);

  constexpr int pack_size = Pack::n;

  add_field<Required>("ps", scalar2d, Pa, grid_name);

  add_field<Updated>("horiz_winds", vector3d_mid, m/s, grid_name, pack_size);
  add_field<Updated>("T_mid", scalar3d_mid, K, grid_name, pack_size);

  add_tracer<Updated>("qv", m_grid, kg/kg, pack_size);
  add_group<Updated>("tracers", grid_name, pack_size, MonolithicAlloc::Required);

  // Sanity check that iop data manager is setup by driver
  EKAT_REQUIRE_MSG(m_iop_data_manager,
    "Error! IOPDataManager not setup by driver, but IOPForcing"
    "being used as an ATM process.\n");

  // Create helper fields for finding horizontal means
  auto level_only_scalar_layout = scalar3d_mid.clone().strip_dim(0);
  auto level_only_vector_layout = vector3d_mid.clone().strip_dim(0);
  const auto iop_nudge_tq = m_iop_data_manager->get_params().get<bool>("iop_nudge_tq");
  const auto iop_nudge_uv = m_iop_data_manager->get_params().get<bool>("iop_nudge_uv");
  if (iop_nudge_tq or iop_nudge_uv) {
    create_helper_field("horiz_mean_weights", scalar2d, grid_name, pack_size);
  }
  if (iop_nudge_tq) {
    create_helper_field("qv_mean", level_only_scalar_layout, grid_name, pack_size);
    create_helper_field("t_mean",  level_only_scalar_layout, grid_name, pack_size);
  }
  if (iop_nudge_uv) {
    create_helper_field("horiz_winds_mean", level_only_vector_layout, grid_name, pack_size);
  }
}
// =========================================================================================
void IOPForcing::
set_computed_group_impl (const FieldGroup& group)
{
  EKAT_REQUIRE_MSG(group.m_info->size() >= 1,
                   "Error! IOPForcing requires at least qv as tracer input.\n");

  const auto& name = group.m_info->m_group_name;

  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! IOPForcing was not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(group.m_info->m_monolithic_allocation,
      "Error! IOPForcing expects a monolithic allocation for tracers.\n");

  m_num_tracers = group.m_info->size();
}
// =========================================================================================
size_t IOPForcing::requested_buffer_size_in_bytes() const
{
  // Number of bytes needed by the WorkspaceManager passed to shoc_main
  const int nlevi_packs  = ekat::npack<Pack>(m_num_levs+1);
  const auto policy      = ESU::get_default_team_policy(m_num_cols, nlevi_packs);
  const size_t wsm_bytes = WorkspaceMgr::get_total_bytes_needed(nlevi_packs, 7+m_num_tracers, policy);

  return wsm_bytes;
}
// =========================================================================================
void IOPForcing::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  const int nlevi_packs = ekat::npack<Pack>(m_num_levs+1);
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());

  // WSM data
  m_buffer.wsm_data = mem;

  const auto policy       = ESU::get_default_team_policy(m_num_cols, nlevi_packs);
  const size_t wsm_npacks = WorkspaceMgr::get_total_bytes_needed(nlevi_packs, 7+m_num_tracers, policy)/sizeof(Pack);
  mem += wsm_npacks;

  size_t used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for IOPForcing.\n");
}
// =========================================================================================
void IOPForcing::create_helper_field (const std::string& name,
                                      const FieldLayout& layout,
                                      const std::string& grid_name,
                                      const int          ps)
{
  using namespace ekat::units;
  FieldIdentifier id(name,layout,Units::nondimensional(),grid_name);

  // Create the field. Init with NaN's, so we spot instances of uninited memory usage
  Field f(id);
  f.get_header().get_alloc_properties().request_allocation(ps);
  f.allocate_view();
  f.deep_copy(ekat::ScalarTraits<Real>::invalid());

  m_helper_fields[name] = f;
}
// =========================================================================================
void IOPForcing::initialize_impl (const RunType run_type)
{
  // Set field property checks for the fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"),m_grid,-400.0,400.0,false);
  // For qv, ensure it doesn't get negative, by allowing repair of any neg value.
  // TODO: use a repairable lb that clips only "small" negative values
  add_postcondition_check<Interval>(get_field_out("qv"),m_grid,0,0.2,true);

  // Setup WSM for internal local variables
  const auto nlevi_packs = ekat::npack<Pack>(m_num_levs+1);
  const auto policy = ESU::get_default_team_policy(m_num_cols, nlevi_packs);
  m_workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 7+m_num_tracers, policy);

  // Compute field for horizontal contraction weights (1/num_global_dofs)
  const auto iop_nudge_tq = m_iop_data_manager->get_params().get<bool>("iop_nudge_tq");
  const auto iop_nudge_uv = m_iop_data_manager->get_params().get<bool>("iop_nudge_uv");
  const Real one_over_num_dofs = 1.0/m_grid->get_num_global_dofs();
  if (iop_nudge_tq or iop_nudge_uv) m_helper_fields.at("horiz_mean_weights").deep_copy(one_over_num_dofs);
}
// =========================================================================================
KOKKOS_FUNCTION
void IOPForcing::
advance_iop_subsidence(const MemberType& team,
                       const int nlevs,
                       const Real dt,
                       const Real ps,
                       const view_1d<const Pack>& ref_p_mid,
                       const view_1d<const Pack>& ref_p_int,
                       const view_1d<const Pack>& ref_p_del,
                       const view_1d<const Pack>& omega,
                       const Workspace& workspace,
                       const view_1d<Pack>& u,
                       const view_1d<Pack>& v,
                       const view_1d<Pack>& T,
                       const view_2d<Pack>& Q)
{
  constexpr Real Rair = C::Rair;
  constexpr Real Cpair = C::Cpair;

  const auto n_q_tracers = Q.extent_int(0);
  const auto nlev_packs = ekat::npack<Pack>(nlevs);

  // Get some temporary views from WS
  uview_1d<Pack> omega_int, delta_u, delta_v, delta_T, tmp;
  workspace.take_many_contiguous_unsafe<4>({"omega_int", "delta_u", "delta_v", "delta_T"},
                                           {&omega_int,  &delta_u,  &delta_v,  &delta_T});
  const auto delta_Q_slot = workspace.take_macro_block("delta_Q", n_q_tracers);
  uview_2d<Pack> delta_Q(delta_Q_slot.data(), n_q_tracers, nlev_packs);

  auto s_ref_p_mid = ekat::scalarize(ref_p_mid);
  auto s_omega = ekat::scalarize(omega);
  auto s_delta_u = ekat::scalarize(delta_u);
  auto s_delta_v = ekat::scalarize(delta_v);
  auto s_delta_T = ekat::scalarize(delta_T);
  auto s_delta_Q = ekat::scalarize(delta_Q);
  auto s_omega_int = ekat::scalarize(omega_int);

  // Compute omega on the interface grid by using a weighted average in pressure
  const int pack_begin = 1/Pack::n, pack_end = (nlevs-1)/Pack::n;
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, pack_begin, pack_end+1), [&] (const int k){
    auto range_pack = ekat::range<IntPack>(k*Pack::n);
    range_pack.set(range_pack<1, 1);
    Pack ref_p_mid_k, ref_p_mid_km1, omega_k, omega_km1;
    ekat::index_and_shift<-1>(s_ref_p_mid, range_pack, ref_p_mid_k, ref_p_mid_km1);
    ekat::index_and_shift<-1>(s_omega, range_pack, omega_k, omega_km1);

    const auto weight = (ref_p_int(k) - ref_p_mid_km1)/(ref_p_mid_k - ref_p_mid_km1);
    omega_int(k).set(range_pack>=1 and range_pack<=nlevs-1,
                      weight*omega_k + (1-weight)*omega_km1);
  });
  omega_int(0)[0] = 0;
  omega_int(nlevs/Pack::n)[nlevs%Pack::n] = 0;

  // Compute delta views for u, v, T, and Q (e.g., u(k+1) - u(k), k=0,...,nlevs-2)
  ColOps::compute_midpoint_delta(team, nlevs-1, u, delta_u);
  ColOps::compute_midpoint_delta(team, nlevs-1, v, delta_v);
  ColOps::compute_midpoint_delta(team, nlevs-1, T, delta_T);
  for (int iq=0; iq<n_q_tracers; ++iq) {
    auto tracer       = Kokkos::subview(Q,       iq, Kokkos::ALL());
    auto delta_tracer = Kokkos::subview(delta_Q, iq, Kokkos::ALL());
    ColOps::compute_midpoint_delta(team, nlevs-1, tracer, delta_tracer);
  }
  team.team_barrier();

  // Compute updated temperature, horizontal winds, and tracers
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
    auto range_pack = ekat::range<IntPack>(k*Pack::n);
    const auto at_top = range_pack==0;
    const auto not_at_top = not at_top;
    const auto at_bot = range_pack==nlevs-1;
    const auto not_at_bot = not at_bot;
    const bool any_at_top = at_top.any();
    const bool any_at_bot = at_bot.any();

    // Get delta(k-1) packs. The range pack should not
    // contain index 0 (so that we don't attempt to access
    // k=-1 index) or index > nlevs-2 (since delta_* views
    // are size nlevs-1).
    auto range_pack_for_m1_shift = range_pack;
    range_pack_for_m1_shift.set(range_pack<1, 1);
    range_pack_for_m1_shift.set(range_pack>nlevs-2, nlevs-2);
    Pack delta_u_k, delta_u_km1,
         delta_v_k, delta_v_km1,
         delta_T_k, delta_T_km1;
    ekat::index_and_shift<-1>(s_delta_u, range_pack_for_m1_shift, delta_u_k, delta_u_km1);
    ekat::index_and_shift<-1>(s_delta_v, range_pack_for_m1_shift, delta_v_k, delta_v_km1);
    ekat::index_and_shift<-1>(s_delta_T, range_pack_for_m1_shift, delta_T_k, delta_T_km1);

    // At the top and bottom of the model, set the end points for
    // delta_*_k and delta_*_km1 to be the first and last entries
    // of delta_*, respectively.
    if (any_at_top) {
      delta_u_k.set(at_top, s_delta_u(0));
      delta_v_k.set(at_top, s_delta_v(0));
      delta_T_k.set(at_top, s_delta_T(0));
    }
    if (any_at_bot) {
      delta_u_km1.set(at_bot, s_delta_u(nlevs-2));
      delta_v_km1.set(at_bot, s_delta_v(nlevs-2));
      delta_T_km1.set(at_bot, s_delta_T(nlevs-2));
    }

    // Get omega_int(k+1) pack. The range pack should not
    // contain index > nlevs-1 (since omega_int is size nlevs+1).
    auto range_pack_for_p1_shift = range_pack;
    range_pack_for_p1_shift.set(range_pack>nlevs-1, nlevs-1);
    Pack omega_int_k, omega_int_kp1;
    ekat::index_and_shift<1>(s_omega_int, range_pack, omega_int_k, omega_int_kp1);

    const auto fac = (dt/2)/ref_p_del(k);

    // Update u
    u(k).update(not_at_bot, fac*omega_int_kp1*delta_u_k, -1, 1);
    u(k).update(not_at_top, fac*omega_int_k*delta_u_km1, -1, 1);

    // Update v
    v(k).update(not_at_bot, fac*omega_int_kp1*delta_v_k, -1, 1);
    v(k).update(not_at_top, fac*omega_int_k*delta_v_km1, -1, 1);

    // Before updating T, first scale using thermal
    // expansion term due to LS vertical advection
    T(k) *= 1 + (dt*Rair/Cpair)*omega(k)/ref_p_mid(k);

    // Update T
    T(k).update(not_at_bot, fac*omega_int_kp1*delta_T_k, -1, 1);
    T(k).update(not_at_top, fac*omega_int_k*delta_T_km1, -1, 1);

    // Update Q
    Pack delta_tracer_k, delta_tracer_km1;
    for (int iq=0; iq<n_q_tracers; ++iq) {
      auto s_delta_tracer = Kokkos::subview(s_delta_Q, iq, Kokkos::ALL());
      ekat::index_and_shift<-1>(s_delta_tracer, range_pack_for_m1_shift, delta_tracer_k, delta_tracer_km1);
      if (any_at_top) delta_tracer_k.set(at_top, s_delta_tracer(0));
      if (any_at_bot) delta_tracer_km1.set(at_bot, s_delta_tracer(nlevs-2));

      Q(iq, k).update(not_at_bot, fac*omega_int_kp1*delta_tracer_k, -1, 1);
      Q(iq, k).update(not_at_top, fac*omega_int_k*delta_tracer_km1, -1, 1);
    }
  });

  // Release WS views
  workspace.release_macro_block(delta_Q_slot, n_q_tracers);
  workspace.release_many_contiguous<4>({&omega_int,  &delta_u,  &delta_v,  &delta_T});
}
// =========================================================================================
KOKKOS_FUNCTION
void IOPForcing::
advance_iop_forcing(const MemberType& team,
                         const int nlevs,
                         const Real dt,
                         const view_1d<const Pack>& divT,
                         const view_1d<const Pack>& divq,
                         const view_1d<Pack>& T,
                         const view_1d<Pack>& qv)
{
  const auto nlev_packs = ekat::npack<Pack>(nlevs);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
    T(k).update(divT(k), dt, 1.0);
    qv(k).update(divq(k), dt, 1.0);
  });
}
// =========================================================================================
KOKKOS_FUNCTION
void IOPForcing::
iop_apply_coriolis(const MemberType& team,
                   const int nlevs,
                   const Real dt,
                   const Real lat,
                   const view_1d<const Pack>& u_ls,
                   const view_1d<const Pack>& v_ls,
                   const view_1d<Pack>& u,
                   const view_1d<Pack>& v)
{
  constexpr Real pi = C::Pi;
  constexpr Real earth_rotation = C::omega;

  // Compute coriolis force
  const auto fcor = 2*earth_rotation*std::sin(lat*pi/180);

  const auto nlev_packs = ekat::npack<Pack>(nlevs);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&] (const int k) {
    const auto u_cor = v(k) - v_ls(k);
    const auto v_cor = u(k) - u_ls(k);
    u(k).update(u_cor, dt*fcor, 1.0);
    v(k).update(v_cor, -dt*fcor, 1.0);
  });
}
// =========================================================================================
void IOPForcing::run_impl (const double dt)
{
  // Pack dimensions
  const auto nlev_packs  = ekat::npack<Pack>(m_num_levs);

  // Hybrid coord values
  const auto ps0 = C::P0;
  const auto hyam = m_grid->get_geometry_data("hyam").get_view<const Real*>();
  const auto hybm = m_grid->get_geometry_data("hybm").get_view<const Real*>();
  const auto hyai = m_grid->get_geometry_data("hyai").get_view<const Real*>();
  const auto hybi = m_grid->get_geometry_data("hybi").get_view<const Real*>();

  // Get FM fields
  const auto ps = get_field_in("ps").get_view<const Real*>();
  const auto horiz_winds = get_field_out("horiz_winds").get_view<Pack***>();
  const auto T_mid = get_field_out("T_mid").get_view<Pack**>();
  const auto qv = get_field_out("qv").get_view<Pack**>();
  const auto Q = get_group_out("tracers").m_monolithic_field->get_view<Pack***>();

  // Load data from IOP files, if necessary
  // TODO: this is using the TS from the beg of the step. Should it use end_of_step_ts() instead?
  m_iop_data_manager->read_iop_file_data(start_of_step_ts());

  // Define local IOP param values
  const auto iop_dosubsidence     = m_iop_data_manager->get_params().get<bool>("iop_dosubsidence");
  const auto iop_coriolis         = m_iop_data_manager->get_params().get<bool>("iop_coriolis");
  const auto iop_nudge_tq         = m_iop_data_manager->get_params().get<bool>("iop_nudge_tq");
  const auto iop_nudge_uv         = m_iop_data_manager->get_params().get<bool>("iop_nudge_uv");
  const auto use_large_scale_wind = m_iop_data_manager->get_params().get<bool>("use_large_scale_wind");
  const auto use_3d_forcing       = m_iop_data_manager->get_params().get<bool>("use_3d_forcing");
  const auto target_lat           = m_iop_data_manager->get_params().get<Real>("target_latitude");
  const auto iop_nudge_tscale     = m_iop_data_manager->get_params().get<Real>("iop_nudge_tscale");
  const auto iop_nudge_tq_low     = m_iop_data_manager->get_params().get<Real>("iop_nudge_tq_low");
  const auto iop_nudge_tq_high    = m_iop_data_manager->get_params().get<Real>("iop_nudge_tq_high");

  // Define local IOP field views
  const Real ps_iop = m_iop_data_manager->get_iop_field("Ps").get_view<const Real, Host>()();
  view_1d<const Pack> omega, divT, divq, u_ls, v_ls, qv_iop, t_iop, u_iop, v_iop;
  divT = use_3d_forcing ? m_iop_data_manager->get_iop_field("divT3d").get_view<const Pack*>()
                        : m_iop_data_manager->get_iop_field("divT").get_view<const Pack*>();
  divq = use_3d_forcing ? m_iop_data_manager->get_iop_field("divq3d").get_view<const Pack*>()
                        : m_iop_data_manager->get_iop_field("divq").get_view<const Pack*>();
  if (iop_dosubsidence) {
    omega = m_iop_data_manager->get_iop_field("omega").get_view<const Pack*>();
  }
  if (iop_coriolis) {
    u_ls = m_iop_data_manager->get_iop_field("u_ls").get_view<const Pack*>();
    v_ls = m_iop_data_manager->get_iop_field("v_ls").get_view<const Pack*>();
  }
  if (iop_nudge_tq) {
    qv_iop = m_iop_data_manager->get_iop_field("q").get_view<const Pack*>();
    t_iop  = m_iop_data_manager->get_iop_field("T").get_view<const Pack*>();
  }
  if (iop_nudge_uv) {
    u_iop = use_large_scale_wind ? m_iop_data_manager->get_iop_field("u_ls").get_view<const Pack*>()
                                 : m_iop_data_manager->get_iop_field("u").get_view<const Pack*>();
    v_iop  = use_large_scale_wind ? m_iop_data_manager->get_iop_field("v_ls").get_view<const Pack*>()
                                  : m_iop_data_manager->get_iop_field("v").get_view<const Pack*>();
  }

  // Team policy and workspace manager for eamxx
  const auto policy_iop = ESU::get_default_team_policy(m_num_cols, nlev_packs);

  // Reset internal WSM variables.
  m_workspace_mgr.reset_internals();

  // Avoid implicit capture of this
  auto wsm = m_workspace_mgr;
  auto num_levs = m_num_levs;

  // Apply IOP forcing
  Kokkos::parallel_for("apply_iop_forcing", policy_iop, KOKKOS_LAMBDA (const MemberType& team) {
    const int icol  =  team.league_rank();

    auto ps_i = ps(icol);
    auto u_i = Kokkos::subview(horiz_winds, icol, 0, Kokkos::ALL());
    auto v_i = Kokkos::subview(horiz_winds, icol, 1, Kokkos::ALL());
    auto T_mid_i = ekat::subview(T_mid, icol);
    auto qv_i = ekat::subview(qv, icol);
    auto Q_i = Kokkos::subview(Q, icol, Kokkos::ALL(), Kokkos::ALL());

    auto ws = wsm.get_workspace(team);
    uview_1d<Pack> ref_p_mid, ref_p_int, ref_p_del;
    ws.take_many_contiguous_unsafe<3>({"ref_p_mid", "ref_p_int", "ref_p_del"},
                                      {&ref_p_mid,  &ref_p_int,  &ref_p_del});

    // Compute reference pressures and layer thickness.
    // TODO: Allow geometry data to allocate packsize
    auto s_ref_p_mid = ekat::scalarize(ref_p_mid);
    auto s_ref_p_int = ekat::scalarize(ref_p_int);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_levs+1), [&](const int& k) {
      s_ref_p_int(k) = hyai(k)*ps0 + hybi(k)*ps_i;
      if (k < num_levs) {
        s_ref_p_mid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      }
    });
    team.team_barrier();
    ColOps::compute_midpoint_delta(team, num_levs, ref_p_int, ref_p_del);
    team.team_barrier();

    if (iop_dosubsidence) {
    // Compute subsidence due to large-scale forcing
      advance_iop_subsidence(team, num_levs, dt, ps_i, ref_p_mid, ref_p_int, ref_p_del, omega, ws, u_i, v_i, T_mid_i, Q_i);
    }

    // Update T and qv according to large scale forcing as specified in IOP file.
    advance_iop_forcing(team, num_levs, dt, divT, divq, T_mid_i, qv_i);

    if (iop_coriolis) {
      // Apply coriolis forcing to u and v winds
      iop_apply_coriolis(team, num_levs, dt, target_lat, u_ls, v_ls, u_i, v_i);
    }

    // Release WS views
    ws.release_many_contiguous<3>({&ref_p_mid, &ref_p_int, &ref_p_del});
  });

  // Nudge the domain based on the domain mean
  // and observed quantities of T, Q, u, and v
  if (iop_nudge_tq or iop_nudge_uv) {
    // Compute domain mean of qv, T_mid, u, and v
    view_1d<Pack> qv_mean, t_mean;
    view_2d<Pack> horiz_winds_mean;
    if (iop_nudge_tq){
      horiz_contraction<Real>(m_helper_fields.at("qv_mean"), get_field_out("qv"),
                              m_helper_fields.at("horiz_mean_weights"), &m_comm);
      qv_mean = m_helper_fields.at("qv_mean").get_view<Pack*>();

      horiz_contraction<Real>(m_helper_fields.at("t_mean"), get_field_out("T_mid"),
                              m_helper_fields.at("horiz_mean_weights"), &m_comm);
      t_mean = m_helper_fields.at("t_mean").get_view<Pack*>();
    }
    if (iop_nudge_uv){
      horiz_contraction<Real>(m_helper_fields.at("horiz_winds_mean"), get_field_out("horiz_winds"),
                              m_helper_fields.at("horiz_mean_weights"), &m_comm);
      horiz_winds_mean = m_helper_fields.at("horiz_winds_mean").get_view<Pack**>();
    }

    // Apply relaxation
    const auto rtau = std::max(dt, iop_nudge_tscale);
    Kokkos::parallel_for("apply_domain_relaxation",
                          policy_iop,
                          KOKKOS_LAMBDA (const MemberType& team) {
      const int icol = team.league_rank();

      auto ps_i = ps(icol);
      auto u_i = Kokkos::subview(horiz_winds, icol, 0, Kokkos::ALL());
      auto v_i = Kokkos::subview(horiz_winds, icol, 1, Kokkos::ALL());
      auto T_mid_i = ekat::subview(T_mid, icol);
      auto qv_i = ekat::subview(qv, icol);

      auto ws = wsm.get_workspace(team);
      uview_1d<Pack> ref_p_mid;
      ws.take_many_contiguous_unsafe<1>({"ref_p_mid"},{&ref_p_mid});

      // Compute reference pressures and layer thickness.
      // TODO: Allow geometry data to allocate packsize
      auto s_ref_p_mid = ekat::scalarize(ref_p_mid);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, num_levs), [&](const int& k) {
        s_ref_p_mid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      });
      team.team_barrier();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int& k) {
        if (iop_nudge_tq) {
          // Restrict nudging of T and qv to certain levels if requested by user
          // IOP pressure variable is in unitis of [Pa], while iop_nudge_tq_low/high
          // is in units of [hPa], thus convert iop_nudge_tq_low/high
          Mask nudge_level(false);
          int max_size = hyam.size();
          for (int lev=k*Pack::n, p = 0; p < Pack::n && lev < max_size; ++lev, ++p) {
            const auto pressure_from_iop = hyam(lev)*ps0 + hybm(lev)*ps_iop;
            nudge_level.set(p, pressure_from_iop <= iop_nudge_tq_low*100
                                and
                                pressure_from_iop >= iop_nudge_tq_high*100);
          }

          qv_i(k).update(nudge_level, qv_mean(k) - qv_iop(k), -dt/rtau, 1.0);
          T_mid_i(k).update(nudge_level, t_mean(k) - t_iop(k), -dt/rtau, 1.0);
        }
        if (iop_nudge_uv) {
          u_i(k).update(horiz_winds_mean(0, k) - u_iop(k), -dt/rtau, 1.0);
          v_i(k).update(horiz_winds_mean(1, k) - v_iop(k), -dt/rtau, 1.0);
        }
      });

      // Release WS views
      ws.release_many_contiguous<1>({&ref_p_mid});
    });
  }
}
// =========================================================================================
} // namespace scream
