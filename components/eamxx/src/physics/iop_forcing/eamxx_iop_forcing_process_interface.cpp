#include "physics/iop_forcing/eamxx_iop_forcing_process_interface.hpp"

#include "share/field/field_utils.hpp"
#include "share/property_checks/field_within_interval_check.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

#include <ekat_math_utils.hpp>

namespace scream
{
// =========================================================================================
void
IOPForcing::create_requests()
{
  using namespace ekat::units;

  m_grid                = m_grids_manager->get_grid("physics");
  const auto &grid_name = m_grid->name();

  m_num_cols = m_grid->get_num_local_dofs();      // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels(); // Number of levels per column

  // Define the different field layouts that will be used for this process
  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true, 2);

  constexpr int pack_size = Pack::n;

  add_field<Required>("ps", scalar2d, Pa, grid_name);

  add_field<Updated>("horiz_winds", vector3d_mid, m / s, grid_name, pack_size);
  add_field<Updated>("T_mid", scalar3d_mid, K, grid_name, pack_size);

  add_tracer<Updated>("qv", m_grid, kg / kg, pack_size);
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
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Number of bytes needed by the WorkspaceManager passed to shoc_main
  const int nlevi_packs  = ekat::npack<Pack>(m_num_levs+1);
  const auto policy      = TPF::get_default_team_policy(m_num_cols, nlevi_packs);
  const size_t wsm_bytes = WorkspaceMgr::get_total_bytes_needed(nlevi_packs, 7+m_num_tracers, policy);

  return wsm_bytes;
}
// =========================================================================================
void IOPForcing::init_buffers(const ATMBufferManager &buffer_manager)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Buffers size not sufficient.\n");

  const int nlevi_packs = ekat::npack<Pack>(m_num_levs+1);
  Pack* mem = reinterpret_cast<Pack*>(buffer_manager.get_memory());

  // WSM data
  m_buffer.wsm_data = mem;

  const auto policy       = TPF::get_default_team_policy(m_num_cols, nlevi_packs);
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
  f.deep_copy(ekat::invalid<Real>());

  m_helper_fields[name] = f;
}
// =========================================================================================
void IOPForcing::initialize_impl (const RunType run_type)
{
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Set field property checks for the fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"),m_grid,-400.0,400.0,false);
  // For qv, ensure it doesn't get negative, by allowing repair of any neg value.
  // TODO: use a repairable lb that clips only "small" negative values
  add_postcondition_check<Interval>(get_field_out("qv"),m_grid,0,0.2,true);

  // Setup WSM for internal local variables
  const auto nlevi_packs = ekat::npack<Pack>(m_num_levs+1);
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlevi_packs);
  m_workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 7+m_num_tracers, policy);

  // Compute field for horizontal contraction weights (1/num_global_dofs)
  const auto iop_nudge_tq = m_iop_data_manager->get_params().get<bool>("iop_nudge_tq");
  const auto iop_nudge_uv = m_iop_data_manager->get_params().get<bool>("iop_nudge_uv");
  const Real one_over_num_dofs = 1.0/m_grid->get_num_global_dofs();
  if (iop_nudge_tq or iop_nudge_uv) m_helper_fields.at("horiz_mean_weights").deep_copy(one_over_num_dofs);
}
// =========================================================================================
// =========================================================================================
// =========================================================================================
// Inline helper for scalar 1D linear interpolation
KOKKOS_FUNCTION
inline Real linear_interp_1d(const Real* x, const Real* f, const int n, const Real x_interp) {
  if (x_interp <= x[0]) return f[0];
  if (x_interp >= x[n-1]) return f[n-1];

  for (int i = 0; i < n-1; ++i) {
    if (x_interp >= x[i] && x_interp <= x[i+1]) {
      Real t = (x_interp - x[i]) / (x[i+1] - x[i]);
      return (1.0 - t) * f[i] + t * f[i+1];
    }
  }
  return f[n-1]; // fallback
}

// =========================================================================================
// Main semi-Lagrangian subsidence routine
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
  constexpr Real Rair = C::Rair.value;
  constexpr Real Cpair = C::Cpair.value;

  const int n_q_tracers = Q.extent_int(0);

  // Scalar views for pack-based inputs
  auto s_ref_p_mid = ekat::scalarize(ref_p_mid);
  auto s_omega     = ekat::scalarize(omega);
  auto s_u         = ekat::scalarize(u);
  auto s_v         = ekat::scalarize(v);
  auto s_T         = ekat::scalarize(T);
  auto s_Q         = ekat::scalarize(Q);

  const Real* x_ptr     = s_ref_p_mid.data();
  const Real* omega_ptr = s_omega.data();
  Real* u_ptr           = s_u.data();
  Real* v_ptr           = s_v.data();
  Real* T_ptr           = s_T.data();
  Real* Q_ptr           = s_Q.data();

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlevs), [&](const int k) {
    const Real p_ref   = x_ptr[k];
    const Real omega_k = omega_ptr[k];
    const Real p_dep   = p_ref - dt * omega_k;

    // Note that I know I should probably be using ekat's linear interp
    //  routine but I couldn't figure out for the life of me how to interface
    //  with this without getting an onslaught of compile issues.  HELP!?

    // Interpolate u, v, T at departure level
    u_ptr[k] = linear_interp_1d(x_ptr, u_ptr, nlevs, p_dep);
    v_ptr[k] = linear_interp_1d(x_ptr, v_ptr, nlevs, p_dep);
    T_ptr[k] = linear_interp_1d(x_ptr, T_ptr, nlevs, p_dep);

    // Add thermal expansion correction
    T_ptr[k] *= 1.0 + (dt * Rair / Cpair) * omega_k / p_ref;

    // Interpolate each tracer
    for (int m = 0; m < n_q_tracers; ++m) {
      Real* tracer_ptr = &Q_ptr[m * nlevs];
      tracer_ptr[k] = linear_interp_1d(x_ptr, tracer_ptr, nlevs, p_dep);
    }
  });
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
  constexpr Real earth_rotation = C::omega.value;

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
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  // Pack dimensions
  const auto nlev_packs  = ekat::npack<Pack>(m_num_levs);

  // Hybrid coord values
  const Real ps0 = C::P0.value;
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
  const auto policy_iop = TPF::get_default_team_policy(m_num_cols, nlev_packs);

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
      horiz_contraction(m_helper_fields.at("qv_mean"), get_field_out("qv"),
                        m_helper_fields.at("horiz_mean_weights"), true, &m_comm);
      qv_mean = m_helper_fields.at("qv_mean").get_view<Pack*>();

      horiz_contraction(m_helper_fields.at("t_mean"), get_field_out("T_mid"),
                        m_helper_fields.at("horiz_mean_weights"), true, &m_comm);
      t_mean = m_helper_fields.at("t_mean").get_view<Pack*>();
    }
    if (iop_nudge_uv){
      horiz_contraction(m_helper_fields.at("horiz_winds_mean"), get_field_out("horiz_winds"),
                        m_helper_fields.at("horiz_mean_weights"), true, &m_comm);
      horiz_winds_mean = m_helper_fields.at("horiz_winds_mean").get_view<Pack**>();
    }

    // Apply relaxation
    const auto rtau = std::max(dt, iop_nudge_tscale);
    Kokkos::parallel_for("apply_domain_relaxation",
                          policy_iop,
                          KOKKOS_LAMBDA (const MemberType& team) {
      const int icol = team.league_rank();

      auto u_i = Kokkos::subview(horiz_winds, icol, 0, Kokkos::ALL());
      auto v_i = Kokkos::subview(horiz_winds, icol, 1, Kokkos::ALL());
      auto T_mid_i = ekat::subview(T_mid, icol);
      auto qv_i = ekat::subview(qv, icol);

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
    });
  }
}
// =========================================================================================
} // namespace scream
