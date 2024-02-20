#include "eamxx_homme_process_interface.hpp"

// EAMxx includes
#include "control/intensive_observation_period.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/homme_dynamics_helpers.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/util/scream_column_ops.hpp"

// Homme includes
#include "Context.hpp"
#include "ColumnOps.hpp"
#include "ElementOps.hpp"
#include "EquationOfState.hpp"
#include "HommexxEnums.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "SimulationParams.hpp"
#include "Types.hpp"

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
    const auto at_bot = range_pack==nlevs-1;
    const auto at_mid = not at_top and not at_bot;
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
    auto& u_k = u(k);
    u_k.set(at_top, u_k - fac*omega_int_kp1*delta_u_k);
    u_k.set(at_bot, u_k - fac*omega_int_k*delta_u_km1);
    u_k.set(at_mid, u_k - fac*(omega_int_kp1*delta_u_k + omega_int_k*delta_u_km1));

    // Update v
    auto& v_k = v(k);
    v_k.set(at_top, v_k - fac*omega_int_kp1*delta_v_k);
    v_k.set(at_bot, v_k - fac*omega_int_k*delta_v_km1);
    v_k.set(at_mid, v_k - fac*(omega_int_kp1*delta_v_k + omega_int_k*delta_v_km1));

    // Before updating T, first scale using thermal
    // expansion term due to LS vertical advection
    auto& T_k = T(k);
    T_k *= 1 + (dt*Rair/Cpair)*omega(k)/pmid(k);

    // Update T
    T_k.set(at_top, T_k - fac*omega_int_kp1*delta_T_k);
    T_k.set(at_bot, T_k - fac*omega_int_k*delta_T_km1);
    T_k.set(at_mid, T_k - fac*(omega_int_kp1*delta_T_k + omega_int_k*delta_T_km1));

    // Update Q
    Pack delta_tracer_k, delta_tracer_km1;
    for (int iq=0; iq<n_q_tracers; ++iq) {
      auto s_delta_tracer = Kokkos::subview(s_delta_Q, iq, Kokkos::ALL());
      ekat::index_and_shift<-1>(s_delta_tracer, range_pack_for_m1_shift, delta_tracer_k, delta_tracer_km1);
      if (any_at_top) delta_tracer_k.set(at_top, s_delta_tracer(0));
      if (any_at_bot) delta_tracer_km1.set(at_bot, s_delta_tracer(nlevs-2));

      auto& Q_k = Q(iq, k);
      Q_k.set(at_top, Q_k - fac*omega_int_kp1*delta_tracer_k);
      Q_k.set(at_bot, Q_k - fac*omega_int_k*delta_tracer_km1);
      Q_k.set(at_mid, Q_k - fac*(omega_int_kp1*delta_tracer_k + omega_int_k*delta_tracer_km1));
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
    T(k) += dt*divT(k);
    qv(k) += dt*divq(k);
  });
}

void HommeDynamics::
apply_iop_forcing(const Real dt)
{
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;

  using EOS = Homme::EquationOfState;
  using ElementOps = Homme::ElementOps;
  using KV = Homme::KernelVariables;

  using ColOps = ColumnOps<DefaultDevice, Real>;
  using C = physics::Constants<Real>;
  constexpr Real Rair = C::Rair;

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

  // Homme element states and EOS/EO classes
  auto ps_dyn = get_internal_field("ps_dyn").get_view<Real***>();
  auto dp3d_dyn = get_internal_field("dp3d_dyn").get_view<Pack****>();
  auto vtheta_dp_dyn = get_internal_field("vtheta_dp_dyn").get_view<Pack****>();
  auto phi_int_dyn = get_internal_field("phi_int_dyn").get_view<Pack****>();
  auto v_dyn = get_internal_field("v_dyn").get_view<Pack*****>();
  auto Q_dyn = m_helper_fields.at("Q_dyn").get_view<Pack*****>();
  auto Qdp_dyn = get_internal_field("Qdp_dyn").get_view<Pack*****>();

  EOS eos;
  eos.init(params.theta_hydrostatic_mode, hvcoord);

  ElementOps elem_ops;
  elem_ops.init(hvcoord);
  const bool use_moisture = (params.moisture == Homme::MoistDry::MOIST);

  // Define local IOP param values and views
  const auto iop_dosubsidence = m_iop->get_params().get<bool>("iop_dosubsidence");
  const auto use_3d_forcing = m_iop->get_params().get<bool>("use_3d_forcing");
  const auto omega = m_iop->get_iop_field("omega").get_view<const Pack*>();
  const auto divT = use_3d_forcing ? m_iop->get_iop_field("divT3d").get_view<const Pack*>()
                                   : m_iop->get_iop_field("divT").get_view<const Pack*>();
  const auto divq = use_3d_forcing ? m_iop->get_iop_field("divq3d").get_view<const Pack*>()
                                   : m_iop->get_iop_field("divq").get_view<const Pack*>();

  // Team policy and workspace manager for both homme and scream
  // related loops. We need separate policies since hommexx functions used here
  // assume they are called inside nested loops for elements and Gaussian points,
  // whereas EAMxx function we use expects a single level of parallelism
  // for elements and Guassian points.
  // TODO: scream::ColumnOps functions could take an arbitary loop boundary
  //       (TeamVectorRange, TeamThreadRange, ThreadVectorRange) so that
  //       all 3 kernel launches here could be combined.
  const auto policy_homme = ESU::get_default_team_policy(nelem, NLEV);
  const auto policy_eamxx = ESU::get_default_team_policy(nelem*NGP*NGP, NLEV);

  // TODO: Create a memory buffer for this class
  //       and add the below WSM and views
  WorkspaceMgr eamxx_wsm(NLEVI, 7+qsize, policy_eamxx);
  WorkspaceMgr homme_wsm(NLEV,  32,      policy_homme);
  view_Nd<Pack, 4>
    temperature("temperature", nelem, NGP, NGP, NLEV),
    exner("exner", nelem, NGP, NGP, NLEV);

  // Preprocess some homme states to get temperature and exner
  Kokkos::parallel_for("compute_t_and_exner", policy_homme, KOKKOS_LAMBDA (const KT::MemberType& team) {
    KV kv(team);
    const int ie  =  team.league_rank();

    // Get temp views from workspace
    auto ws = homme_wsm.get_workspace(team);
    auto pnh_slot   = ws.take_macro_block("pnh"  , NGP*NGP);
    auto rstar_slot = ws.take_macro_block("rstar", NGP*NGP);
    uview_2d<Pack>
      pnh  (reinterpret_cast<Pack*>(pnh_slot.data()),   NGP*NGP, NLEV),
      rstar(reinterpret_cast<Pack*>(rstar_slot.data()), NGP*NGP, NLEV);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NGP*NGP), [&] (const int idx) {
      const int igp = idx/NGP;
      const int jgp = idx%NGP;

      auto dp3d_i      = ekat::subview(dp3d_dyn, ie, igp, jgp);
      auto vtheta_dp_i = ekat::subview(vtheta_dp_dyn, ie, igp, jgp);
      auto phi_int_i   = ekat::subview(phi_int_dyn, ie, igp, jgp);
      auto qv_i        = ekat::subview(Q_dyn, ie, 0, igp, jgp);
      auto pnh_i         = ekat::subview(pnh, idx);
      auto rstar_i       = ekat::subview(rstar, idx);
      auto exner_i       = ekat::subview(exner, ie, igp, jgp);
      auto temperature_i = ekat::subview(temperature, ie, igp, jgp);

      // Reinterperate into views of Homme::Scalar for calling Hommexx function.
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> dp3d_scalar(reinterpret_cast<Homme::Scalar*>(dp3d_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> vtheta_dp_scalar(reinterpret_cast<Homme::Scalar*>(vtheta_dp_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEVI]> phi_int_scalar(reinterpret_cast<Homme::Scalar*>(phi_int_i.data()), NLEVI);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> qv_scalar(reinterpret_cast<Homme::Scalar*>(qv_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> pnh_scalar(reinterpret_cast<Homme::Scalar*>(pnh_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> exner_scalar(reinterpret_cast<Homme::Scalar*>(exner_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> rstar_scalar(reinterpret_cast<Homme::Scalar*>(rstar_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> temperature_scalar(reinterpret_cast<Homme::Scalar*>(temperature_i.data()), NLEV);

      // Compute exner from EOS
      if (params.theta_hydrostatic_mode) {
        auto hydro_p_int = ws.take("hydro_p_int");
        Homme::ExecViewUnmanaged<Homme::Scalar[NLEVI]> hydro_p_int_scalar(reinterpret_cast<Homme::Scalar*>(hydro_p_int.data()), NLEVI);
        elem_ops.compute_hydrostatic_p(kv, dp3d_scalar, hydro_p_int_scalar, pnh_scalar);
        eos.compute_exner(kv, pnh_scalar, exner_scalar);
        ws.release(hydro_p_int);
      } else {
        eos.compute_pnh_and_exner(kv, vtheta_dp_scalar, phi_int_scalar, pnh_scalar, exner_scalar);
      }

      // Get the temperature from dynamics states
      elem_ops.get_temperature(kv, eos, use_moisture, dp3d_scalar, exner_scalar, vtheta_dp_scalar, qv_scalar, rstar_scalar, temperature_scalar);
    });

    // Release WS views
    ws.release_macro_block(rstar_slot, NGP*NGP);
    ws.release_macro_block(pnh_slot, NGP*NGP);
  });
  Kokkos::fence();

  // Apply IOP forcing
  Kokkos::parallel_for("apply_iop_forcing", policy_eamxx, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank()/(NGP*NGP);
    const int igp = (team.league_rank()/NGP)%NGP;
    const int jgp =  team.league_rank()%NGP;

    // Get temp views from workspace
    auto ws = eamxx_wsm.get_workspace(team);
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

    // Release WS views
    ws.release_many_contiguous<3>({&pmid, &pint, &pdel});
  });
  Kokkos::fence();

  // Postprocess homme states Qdp and vtheta_dp
  Kokkos::parallel_for("compute_qdp_and_vtheta_dp", policy_homme, KOKKOS_LAMBDA (const KT::MemberType& team) {
    KV kv(team);
    const int ie  =  team.league_rank();

    // Get temp views from workspace
    auto ws = homme_wsm.get_workspace(team);
    auto rstar_slot = ws.take_macro_block("rstar", NGP*NGP);
    uview_2d<Pack>
      rstar(reinterpret_cast<Pack*>(rstar_slot.data()), NGP*NGP, NLEV);

    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NGP*NGP), [&] (const int idx) {
      const int igp = idx/NGP;
      const int jgp = idx%NGP;

      auto dp3d_i      = ekat::subview(dp3d_dyn, ie, igp, jgp);
      auto vtheta_dp_i = ekat::subview(vtheta_dp_dyn, ie, igp, jgp);
      auto qv_i        = ekat::subview(Q_dyn, ie, 0, igp, jgp);
      auto Q_i         = Kokkos::subview(Q_dyn, ie, Kokkos::ALL(), igp, jgp, Kokkos::ALL());
      auto Qdp_i       = Kokkos::subview(Qdp_dyn, ie, Kokkos::ALL(), igp, jgp, Kokkos::ALL());
      auto rstar_i = ekat::subview(rstar, idx);
      auto exner_i       = ekat::subview(exner, ie, igp, jgp);
      auto temperature_i = ekat::subview(temperature, ie, igp, jgp);

      // Reinterperate into views of Homme::Scalar for calling Hommexx function.
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> qv_scalar(reinterpret_cast<Homme::Scalar*>(qv_i.data()), NLEV);
      Homme::ExecViewUnmanaged<Homme::Scalar[NLEV]> rstar_scalar(reinterpret_cast<Homme::Scalar*>(rstar_i.data()), NLEV);

      // Compute Qdp from updated Q
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, NLEV*qsize), [&] (const int k) {
        const int ilev = k/qsize;
        const int q = k%qsize;

        Qdp_i(q, ilev) = Q_i(q, ilev)*dp3d_i(ilev);
        // For BFB on restarts, Q needs to be updated after we compute Qdp
        // TODO: Is this needed?
        Q_i(q, ilev) = Qdp_i(q, ilev)/dp3d_i(ilev);
      });

      // Convert updated temperature back to potential temperature
      elem_ops.get_R_star(kv, use_moisture, qv_scalar, rstar_scalar);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, NLEV), [&] (const int k) {
        vtheta_dp_i(k) = temperature_i(k)*rstar_i(k)*dp3d_i(k)/(Rair*exner_i(k));
      });
    });

    // Release WS views
    ws.release_macro_block(rstar_slot, NGP*NGP);
  });
}

} // namespace scream
