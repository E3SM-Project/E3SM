#include "eamxx_homme_process_interface.hpp"

// EAMxx includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/iop/intensive_observation_period.hpp"
#include "share/util/scream_column_ops.hpp"

// Homme includes
#include "Context.hpp"
#include "ColumnOps.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"

// SCREAM includes
#include "share/util/scream_common_physics_functions.hpp"

// EKAT includes
#include "ekat/ekat_workspace.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

namespace scream {

// Compute effects of large scale subsidence on T, q, u, and v.
KOKKOS_FUNCTION
void HommeDynamics::
advance_iop_subsidence(const KT::MemberType& team,
                       const int nlevs,
                       const Real dt,
                       const Real ps,
                       const view_1d<const Pack>& pmid,
                       const view_1d<const Pack>& pint,
                       const view_1d<const Pack>& pdel,
                       const view_1d<const Pack>& omega,
                       const Workspace& workspace,
                       const view_1d<Pack>& u,
                       const view_1d<Pack>& v,
                       const view_1d<Pack>& T,
                       const view_2d<Pack>& Q)
{
  using ColOps = ColumnOps<DefaultDevice, Real>;
  using C = physics::Constants<Real>;
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

  auto s_pmid = ekat::scalarize(pmid);
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
    Pack pmid_k, pmid_km1, omega_k, omega_km1;
    ekat::index_and_shift<-1>(s_pmid, range_pack, pmid_k, pmid_km1);
    ekat::index_and_shift<-1>(s_omega, range_pack, omega_k, omega_km1);

    const auto weight = (pint(k) - pmid_km1)/(pmid_k - pmid_km1);
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

    const auto fac = (dt/2)/pdel(k);

    // Update u
    u(k).update(not_at_bot, fac*omega_int_kp1*delta_u_k, -1, 1);
    u(k).update(not_at_top, fac*omega_int_k*delta_u_km1, -1, 1);

    // Update v
    v(k).update(not_at_bot, fac*omega_int_kp1*delta_v_k, -1, 1);
    v(k).update(not_at_top, fac*omega_int_k*delta_v_km1, -1, 1);

    // Before updating T, first scale using thermal
    // expansion term due to LS vertical advection
    T(k) *= 1 + (dt*Rair/Cpair)*omega(k)/pmid(k);

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

// Apply large scale forcing for temperature and water vapor as provided by the IOP file
KOKKOS_FUNCTION
void HommeDynamics::
advance_iop_forcing(const KT::MemberType& team,
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

// Provide coriolis forcing to u and v winds, using large scale winds specified in IOP forcing file.
KOKKOS_FUNCTION
void HommeDynamics::
iop_apply_coriolis(const KT::MemberType& team,
                   const int nlevs,
                   const Real dt,
                   const Real lat,
                   const view_1d<const Pack>& u_ls,
                   const view_1d<const Pack>& v_ls,
                   const view_1d<Pack>& u,
                   const view_1d<Pack>& v)
{
  using C = physics::Constants<Real>;
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

void HommeDynamics::
apply_iop_forcing(const Real dt)
{
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using PF = PhysicsFunctions<DefaultDevice>;
  using ColOps = ColumnOps<DefaultDevice, Real>;

  // Homme objects
  const auto& c = Homme::Context::singleton();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto& params = c.get<Homme::SimulationParams>();

  // Dimensions
  constexpr int NGP   = HOMMEXX_NP;
  constexpr int NLEV  = HOMMEXX_NUM_LEV;
  constexpr int NLEVI = HOMMEXX_NUM_LEV_P;
  const auto nelem    = m_dyn_grid->get_num_local_dofs()/(NGP*NGP);
  const auto total_levels = m_dyn_grid->get_num_vertical_levels();
  const auto qsize = params.qsize;

  // Sanity checks since we will be switching between ekat::Pack
  // and Homme::Scalar view types
  EKAT_ASSERT_MSG(NLEV  == ekat::npack<Pack>(total_levels),
    "Error! Dimension for vectorized Homme levels does not match level dimension "
    "of the packed views used here. Check that Pack typedef is using a pack size "
    "consistent with Homme's vector size.\n");
  EKAT_ASSERT_MSG(NLEVI == ekat::npack<Pack>(total_levels+1),
    "Error! Dimension for vectorized Homme levels does not match level dimension "
    "of the packed views used here. Check that Pack typedef is using a pack size "
    "consistent with Homme's vector size.\n");

  // Hybrid coord values
  const auto ps0 = hvcoord.ps0;
  const auto hyam = m_dyn_grid->get_geometry_data("hyam").get_view<const Real*>();
  const auto hybm = m_dyn_grid->get_geometry_data("hybm").get_view<const Real*>();
  const auto hyai = m_dyn_grid->get_geometry_data("hyai").get_view<const Real*>();
  const auto hybi = m_dyn_grid->get_geometry_data("hybi").get_view<const Real*>();

  // Homme element states
  auto ps_dyn = get_internal_field("ps_dyn").get_view<Real***>();
  auto dp3d_dyn = get_internal_field("dp3d_dyn").get_view<Pack****>();
  auto vtheta_dp_dyn = get_internal_field("vtheta_dp_dyn").get_view<Pack****>();
  auto phi_int_dyn = get_internal_field("phi_int_dyn").get_view<Pack****>();
  auto v_dyn = get_internal_field("v_dyn").get_view<Pack*****>();
  auto Q_dyn = m_helper_fields.at("Q_dyn").get_view<Pack*****>();
  auto Qdp_dyn = get_internal_field("Qdp_dyn").get_view<Pack*****>();

  // Load data from IOP files, if necessary
  m_iop->read_iop_file_data(timestamp());

  // Define local IOP param values
  const auto iop_dosubsidence = m_iop->get_params().get<bool>("iop_dosubsidence");
  const auto iop_coriolis = m_iop->get_params().get<bool>("iop_coriolis");
  const auto iop_nudge_tq = m_iop->get_params().get<bool>("iop_nudge_tq");
  const auto iop_nudge_uv = m_iop->get_params().get<bool>("iop_nudge_uv");
  const auto use_large_scale_wind = m_iop->get_params().get<bool>("use_large_scale_wind");
  const auto use_3d_forcing = m_iop->get_params().get<bool>("use_3d_forcing");
  const auto lat = m_iop->get_params().get<Real>("target_latitude");
  const auto iop_nudge_tscale = m_iop->get_params().get<Real>("iop_nudge_tscale");
  const auto iop_nudge_tq_low = m_iop->get_params().get<Real>("iop_nudge_tq_low");
  const auto iop_nudge_tq_high = m_iop->get_params().get<Real>("iop_nudge_tq_high");

  // Define local IOP field views
  const Real ps_iop = m_iop->get_iop_field("Ps").get_view<const Real, Host>()();
  view_1d<const Pack> omega, divT, divq, u_ls, v_ls, qv_iop, t_iop, u_iop, v_iop;
  divT = use_3d_forcing ? m_iop->get_iop_field("divT3d").get_view<const Pack*>()
                        : m_iop->get_iop_field("divT").get_view<const Pack*>();
  divq = use_3d_forcing ? m_iop->get_iop_field("divq3d").get_view<const Pack*>()
                        : m_iop->get_iop_field("divq").get_view<const Pack*>();
  if (iop_dosubsidence) {
    omega = m_iop->get_iop_field("omega").get_view<const Pack*>();
  }
  if (iop_coriolis) {
    u_ls = m_iop->get_iop_field("u_ls").get_view<const Pack*>();
    v_ls = m_iop->get_iop_field("v_ls").get_view<const Pack*>();
  }
  if (iop_nudge_tq) {
    qv_iop = m_iop->get_iop_field("q").get_view<const Pack*>();
    t_iop  = m_iop->get_iop_field("T").get_view<const Pack*>();
  }
  if (iop_nudge_uv) {
    u_iop = use_large_scale_wind ? m_iop->get_iop_field("u_ls").get_view<const Pack*>()
                                 : m_iop->get_iop_field("u").get_view<const Pack*>();
    v_iop  = use_large_scale_wind ? m_iop->get_iop_field("v_ls").get_view<const Pack*>()
                                  : m_iop->get_iop_field("v").get_view<const Pack*>();
  }

  // Team policy and workspace manager for eamxx
  const auto policy_iop = ESU::get_default_team_policy(nelem*NGP*NGP, NLEV);

  // TODO: Create a memory buffer for this class
  //       and add the below WSM and views
  WorkspaceMgr iop_wsm(NLEVI, 7+qsize, policy_iop);
  view_Nd<Pack, 4>
    temperature("temperature", nelem, NGP, NGP, NLEV);

  // Lambda for computing temperature
  auto compute_temperature = [&] () {
    Kokkos::parallel_for("compute_temperature_for_iop", policy_iop, KOKKOS_LAMBDA (const KT::MemberType& team) {
      const int ie  =  team.league_rank()/(NGP*NGP);
      const int igp = (team.league_rank()/NGP)%NGP;
      const int jgp =  team.league_rank()%NGP;

      // Get temp views from workspace
      auto ws = iop_wsm.get_workspace(team);
      uview_1d<Pack> pmid;
      ws.take_many_contiguous_unsafe<1>({"pmid"},{&pmid});

      auto ps_i          = ps_dyn(ie, igp, jgp);
      auto dp3d_i        = ekat::subview(dp3d_dyn, ie, igp, jgp);
      auto vtheta_dp_i   = ekat::subview(vtheta_dp_dyn, ie, igp, jgp);
      auto qv_i          = ekat::subview(Q_dyn, ie, 0, igp, jgp);
      auto temperature_i = ekat::subview(temperature, ie, igp, jgp);

      // Compute reference pressures and layer thickness.
      // TODO: Allow geometry data to allocate packsize
      auto s_pmid = ekat::scalarize(pmid);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, total_levels), [&](const int& k) {
        s_pmid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      });
      team.team_barrier();

      // Compute temperature from virtual potential temperature
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, NLEV), [&] (const int k) {
        auto T_val = vtheta_dp_i(k);
        T_val /= dp3d_i(k);
        T_val = PF::calculate_temperature_from_virtual_temperature(T_val,qv_i(k));
        temperature_i(k) = PF::calculate_T_from_theta(T_val,pmid(k));
      });

      // Release WS views
      ws.release_many_contiguous<1>({&pmid});
    });
  };

  // Preprocess some homme states to get temperature
  compute_temperature();
  Kokkos::fence();

  // Apply IOP forcing
  Kokkos::parallel_for("apply_iop_forcing", policy_iop, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank()/(NGP*NGP);
    const int igp = (team.league_rank()/NGP)%NGP;
    const int jgp =  team.league_rank()%NGP;

    // Get temp views from workspace
    auto ws = iop_wsm.get_workspace(team);
    uview_1d<Pack> pmid, pint, pdel;
    ws.take_many_contiguous_unsafe<3>({"pmid", "pint", "pdel"},
                                      {&pmid,  &pint,  &pdel});

    auto ps_i = ps_dyn(ie, igp, jgp);
    auto u_i = ekat::subview(v_dyn, ie, 0, igp, jgp);
    auto v_i = ekat::subview(v_dyn, ie, 1, igp, jgp);
    auto temperature_i = ekat::subview(temperature, ie, igp, jgp);
    auto qv_i = ekat::subview(Q_dyn, ie, 0, igp, jgp);
    auto Q_i = Kokkos::subview(Q_dyn, ie, Kokkos::ALL(), igp, jgp, Kokkos::ALL());

    // Compute reference pressures and layer thickness.
    // TODO: Allow geometry data to allocate packsize
    auto s_pmid = ekat::scalarize(pmid);
    auto s_pint = ekat::scalarize(pint);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, total_levels+1), [&](const int& k) {
      s_pint(k) = hyai(k)*ps0 + hybi(k)*ps_i;
      if (k < total_levels) {
        s_pmid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      }
    });
    team.team_barrier();
    ColOps::compute_midpoint_delta(team, total_levels, pint, pdel);
    team.team_barrier();

    if (iop_dosubsidence) {
    // Compute subsidence due to large-scale forcing
      advance_iop_subsidence(team, total_levels, dt, ps_i, pmid, pint, pdel, omega, ws, u_i, v_i, temperature_i, Q_i);
    }

    // Update T and qv according to large scale forcing as specified in IOP file.
    advance_iop_forcing(team, total_levels, dt, divT, divq, temperature_i, qv_i);

    if (iop_coriolis) {
      // Apply coriolis forcing to u and v winds
      iop_apply_coriolis(team, total_levels, dt, lat, u_ls, v_ls, u_i, v_i);
    }

    // Release WS views
    ws.release_many_contiguous<3>({&pmid, &pint, &pdel});
  });
  Kokkos::fence();

  // Postprocess homme states Qdp and vtheta_dp
  Kokkos::parallel_for("compute_qdp_and_vtheta_dp", policy_iop, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank()/(NGP*NGP);
    const int igp = (team.league_rank()/NGP)%NGP;
    const int jgp =  team.league_rank()%NGP;

    // Get temp views from workspace
    auto ws = iop_wsm.get_workspace(team);
    uview_1d<Pack> pmid, pint, pdel;
    ws.take_many_contiguous_unsafe<3>({"pmid", "pint", "pdel"},
                                      {&pmid,  &pint,  &pdel});

    auto ps_i          = ps_dyn(ie, igp, jgp);
    auto dp3d_i        = ekat::subview(dp3d_dyn, ie, igp, jgp);
    auto vtheta_dp_i   = ekat::subview(vtheta_dp_dyn, ie, igp, jgp);
    auto qv_i          = ekat::subview(Q_dyn, ie, 0, igp, jgp);
    auto Q_i           = Kokkos::subview(Q_dyn, ie, Kokkos::ALL(), igp, jgp, Kokkos::ALL());
    auto Qdp_i         = Kokkos::subview(Qdp_dyn, ie, Kokkos::ALL(), igp, jgp, Kokkos::ALL());
    auto temperature_i = ekat::subview(temperature, ie, igp, jgp);

    // Compute reference pressures and layer thickness.
    // TODO: Allow geometry data to allocate packsize
    auto s_pmid = ekat::scalarize(pmid);
    auto s_pint = ekat::scalarize(pint);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, total_levels+1), [&](const int& k) {
      s_pint(k) = hyai(k)*ps0 + hybi(k)*ps_i;
      if (k < total_levels) {
        s_pmid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      }
    });

    team.team_barrier();

    // Compute Qdp from updated Q
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, NLEV*qsize), [&] (const int k) {
      const int ilev = k/qsize;
      const int q = k%qsize;

      Qdp_i(q, ilev) = Q_i(q, ilev)*dp3d_i(ilev);
      // For BFB on restarts, Q needs to be updated after we compute Qdp
      Q_i(q, ilev) = Qdp_i(q, ilev)/dp3d_i(ilev);
    });
      team.team_barrier();

    // Convert updated temperature back to psuedo density virtual potential temperature
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, NLEV), [&] (const int k) {
        const auto th = PF::calculate_theta_from_T(temperature_i(k),pmid(k));
        vtheta_dp_i(k) = PF::calculate_virtual_temperature(th,qv_i(k))*dp3d_i(k);
    });

    // Release WS views
    ws.release_many_contiguous<3>({&pmid, &pint, &pdel});
  });

  if (iop_nudge_tq or iop_nudge_uv) {
    // Nudge the domain based on the domain mean
    // and observed quantities of T, Q, u, and

    if (iop_nudge_tq) {
      // Compute temperature
      compute_temperature();
      Kokkos::fence();
    }

    // Compute domain mean of qv, temperature, u, and v

    // TODO: add to local mem buffer
    view_1d<Pack> qv_mean, t_mean, u_mean, v_mean;
    if (iop_nudge_tq) {
      qv_mean = view_1d<Pack>("u_mean", NLEV),
      t_mean = view_1d<Pack>("v_mean", NLEV);
    }
    if (iop_nudge_uv){
      u_mean = view_1d<Pack>("u_mean", NLEV),
      v_mean = view_1d<Pack>("v_mean", NLEV);
    }

    const auto qv_mean_h = Kokkos::create_mirror_view(qv_mean);
    const auto t_mean_h  = Kokkos::create_mirror_view(t_mean);
    const auto u_mean_h  = Kokkos::create_mirror_view(u_mean);
    const auto v_mean_h  = Kokkos::create_mirror_view(v_mean);

    for (int k=0; k<total_levels; ++k) {
      if (iop_nudge_tq){
        Real& qv_mean_k = qv_mean_h(k/Pack::n)[k%Pack::n];
        Real& t_mean_k = t_mean_h(k/Pack::n)[k%Pack::n];
        Kokkos::parallel_reduce("compute_domain_means_tq",
                                nelem*NGP*NGP,
                                KOKKOS_LAMBDA (const int idx, Real& q_sum, Real& t_sum) {
          const int ie  =  idx/(NGP*NGP);
          const int igp = (idx/NGP)%NGP;
          const int jgp =  idx%NGP;

          q_sum += Q_dyn(ie, 0, igp, jgp, k/Pack::n)[k%Pack::n];
          t_sum += temperature(ie, igp, jgp, k/Pack::n)[k%Pack::n];
        },
        qv_mean_k,
        t_mean_k);

        m_comm.all_reduce(&qv_mean_k, 1, MPI_SUM);
        m_comm.all_reduce(&t_mean_k, 1, MPI_SUM);

        qv_mean_k /= m_dyn_grid->get_num_global_dofs();
        t_mean_k /= m_dyn_grid->get_num_global_dofs();
      }
      if (iop_nudge_uv){
        Real& u_mean_k = u_mean_h(k/Pack::n)[k%Pack::n];
        Real& v_mean_k = v_mean_h(k/Pack::n)[k%Pack::n];
        Kokkos::parallel_reduce("compute_domain_means_uv",
                                nelem*NGP*NGP,
                                KOKKOS_LAMBDA (const int idx, Real& u_sum, Real& v_sum) {
          const int ie  =  idx/(NGP*NGP);
          const int igp = (idx/NGP)%NGP;
          const int jgp =  idx%NGP;

          u_sum += v_dyn(ie, 0, igp, jgp, k/Pack::n)[k%Pack::n];
          v_sum += v_dyn(ie, 1, igp, jgp, k/Pack::n)[k%Pack::n];
        },
        u_mean_k,
        v_mean_k);

        m_comm.all_reduce(&u_mean_k, 1, MPI_SUM);
        m_comm.all_reduce(&v_mean_k, 1, MPI_SUM);

        u_mean_k /= m_dyn_grid->get_num_global_dofs();
        v_mean_k /= m_dyn_grid->get_num_global_dofs();
      }
    }
    Kokkos::deep_copy(qv_mean, qv_mean_h);
    Kokkos::deep_copy(t_mean,  t_mean_h);
    Kokkos::deep_copy(u_mean,  u_mean_h);
    Kokkos::deep_copy(v_mean,  v_mean_h);

    // Apply relaxation
    const auto rtau = std::max(dt, iop_nudge_tscale);
    Kokkos::parallel_for("apply_domain_relaxation",
                          policy_iop,
                          KOKKOS_LAMBDA (const KT::MemberType& team) {

      const int ie  =  team.league_rank()/(NGP*NGP);
      const int igp = (team.league_rank()/NGP)%NGP;
      const int jgp =  team.league_rank()%NGP;

      // Get temp views from workspace
      auto ws = iop_wsm.get_workspace(team);
      uview_1d<Pack> pmid;
      ws.take_many_contiguous_unsafe<1>({"pmid"},{&pmid});

      auto ps_i          = ps_dyn(ie, igp, jgp);
      auto dp3d_i        = ekat::subview(dp3d_dyn, ie, igp, jgp);
      auto vtheta_dp_i   = ekat::subview(vtheta_dp_dyn, ie, igp, jgp);
      auto qv_i          = ekat::subview(Q_dyn, ie, 0, igp, jgp);
      auto temperature_i = ekat::subview(temperature, ie, igp, jgp);
      auto u_i           = ekat::subview(v_dyn, ie, 0, igp, jgp);
      auto v_i           = ekat::subview(v_dyn, ie, 1, igp, jgp);

      // Compute reference pressures and layer thickness.
      // TODO: Allow geometry data to allocate packsize
      auto s_pmid = ekat::scalarize(pmid);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, total_levels), [&](const int& k) {
        s_pmid(k) = hyam(k)*ps0 + hybm(k)*ps_i;
      });
      team.team_barrier();

      if (iop_nudge_tq or iop_nudge_uv) {
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, NLEV), [&](const int& k) {
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
            temperature_i(k).update(nudge_level, t_mean(k) - t_iop(k), -dt/rtau, 1.0);

            // Convert updated temperature back to virtual potential temperature
            const auto th = PF::calculate_theta_from_T(temperature_i(k),pmid(k));
            vtheta_dp_i(k) = PF::calculate_virtual_temperature(th,qv_i(k))*dp3d_i(k);
          }
          if (iop_nudge_uv) {
            u_i(k).update(u_mean(k) - u_iop(k), -dt/rtau, 1.0);
            v_i(k).update(v_mean(k) - v_iop(k), -dt/rtau, 1.0);
          }
        });
      }

    // Release WS views
    ws.release_many_contiguous<1>({&pmid});
    });
  }
}

} // namespace scream
