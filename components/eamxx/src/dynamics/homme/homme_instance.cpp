#include "dynamics/homme/homme_instance.hpp"

// HOMMEXX Includes
#include "Context.hpp"
#include "Elements.hpp"
#include "ElementsState.hpp"
#include "HommexxEnums.hpp"
#include "SimulationParams.hpp"
#include "ElementsForcing.hpp"
#include "EulerStepFunctor.hpp"
#include "ComposeTransport.hpp"
#include "Diagnostics.hpp"
#include "DirkFunctor.hpp"
#include "ForcingFunctor.hpp"
#include "CaarFunctor.hpp"
#include "VerticalRemapManager.hpp"
#include "HyperviscosityFunctor.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"
#include "mpi/ConnectivityHelpers.hpp"
#include "utilities/MathUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "ExecSpaceDefs.hpp"
#include "Types.hpp"

// Scream includes
#include "dynamics/homme/homme_dimensions.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_column_ops.hpp"

// Ekat includes
#include <ekat_assert.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_subview_utils.hpp>
#include <ekat_pack.hpp>
#include <ekat_math_utils.hpp>
#include <ekat_pack_utils.hpp>
#include <ekat_workspace.hpp>

namespace scream
{

std::shared_ptr<HommeInstance> HommeInstance::get ()
{
  static std::weak_ptr<HommeInstance> weak_instance;

  auto instance = weak_instance.lock();
  if (not instance) {
    // The constructor is private, so we cannot use make_shared
    instance.reset(new HommeInstance());
    weak_instance = instance;
  }
  return instance;
}

void HommeInstance::add_dyn_field (const std::string& name, const Field& f)
{
  EKAT_REQUIRE_MSG (f.get_header().get_identifier().get_grid_name()=="dynamics",
      "Error! HommeInstance can only store fields on the dynamics grid.\n"
      " - field name: " + name + "\n"
      " - grid name : " + f.get_header().get_identifier().get_grid_name() + "\n");
  m_dyn_fields[name] = f;
}

Field HommeInstance::get_dyn_field (const std::string& name) const
{
  EKAT_REQUIRE_MSG (has_dyn_field(name),
      "Error! The dyn field '" + name + "' was not added to the HommeInstance.\n");
  return m_dyn_fields.at(name);
}

void HommeInstance::init_views ()
{
  const auto& c = Homme::Context::singleton();
  auto& params  = c.get<Homme::SimulationParams>();
  auto& state   = c.get<Homme::ElementsState>();
  auto& derived = c.get<Homme::ElementsDerivedState>();
  auto& tracers = c.get<Homme::Tracers>();
  auto& forcing = c.get<Homme::ElementsForcing>();

  constexpr int NTL  = HOMMEXX_NUM_TIME_LEVELS;
  constexpr int QTL  = HOMMEXX_Q_NUM_TIME_LEVELS;
  constexpr int QSZ  = HOMMEXX_QSIZE_D;
  constexpr int NVL  = HOMMEXX_NUM_LEV;
  constexpr int NVLI = HOMMEXX_NUM_LEV_P;

  const int nelem = get_dyn_field("v_dyn").get_header().get_identifier().get_layout().dim(0);
  const int qsize = tracers.num_tracers();

  // Velocity
  auto v_in = get_dyn_field("v_dyn").get_view<Homme::Scalar*[NTL][2][NP][NP][NVL]>();
  using v_type = std::remove_reference<decltype(state.m_v)>::type;
  state.m_v = v_type (v_in.data(),nelem);

  // Virtual potential temperature
  auto vtheta_in = get_dyn_field("vtheta_dp_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using vtheta_type = std::remove_reference<decltype(state.m_vtheta_dp)>::type;
  state.m_vtheta_dp = vtheta_type(vtheta_in.data(),nelem);

  // Geopotential
  auto phi_in = get_dyn_field("phi_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using phi_type = std::remove_reference<decltype(state.m_phinh_i)>::type;
  state.m_phinh_i = phi_type(phi_in.data(),nelem);

  // Vertical velocity
  auto w_in = get_dyn_field("w_int_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVLI]>();
  using w_type = std::remove_reference<decltype(state.m_w_i)>::type;
  state.m_w_i = w_type(w_in.data(),nelem);

  // Pseudo-density
  auto dp3d_in = get_dyn_field("dp3d_dyn").get_view<Homme::Scalar*[NTL][NP][NP][NVL]>();
  using dp3d_type = std::remove_reference<decltype(state.m_dp3d)>::type;
  state.m_dp3d = dp3d_type(dp3d_in.data(),nelem);

  // Surface pressure
  auto ps_in = get_dyn_field("ps_dyn").get_view<Real*[NTL][NP][NP]>();
  using ps_type = std::remove_reference<decltype(state.m_ps_v)>::type;
  state.m_ps_v = ps_type(ps_in.data(),nelem);

  // Vertical pressure velocity
  auto omega_in = get_dyn_field("omega_dyn").get_view<Homme::Scalar*[NP][NP][NVL]>();
  using omega_type = std::remove_reference<decltype(derived.m_omega_p)>::type;
  derived.m_omega_p = omega_type(omega_in.data(),nelem);

  // Tracers mixing ratio
  auto q_in = get_dyn_field("Q_dyn").get_view<Homme::Scalar**[NP][NP][NVL]>();
  using q_type = std::remove_reference<decltype(tracers.Q)>::type;
  tracers.Q = q_type(q_in.data(),nelem,qsize);

  // Tracers mass
  auto qdp_in = get_dyn_field("Qdp_dyn").get_view<Homme::Scalar*[QTL][QSZ][NP][NP][NVL]>();
  using qdp_type = std::remove_reference<decltype(tracers.qdp)>::type;
  tracers.qdp = qdp_type(qdp_in.data(),nelem);

  // Tracers forcing
  auto fq_in = get_dyn_field("FQ_dyn").get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fq_type = std::remove_reference<decltype(tracers.fq)>::type;
  tracers.fq = fq_type(fq_in.data(),nelem,qsize);

  // Temperature forcing
  auto ft_in = get_dyn_field("FT_dyn").get_view<Homme::Scalar*[NP][NP][NVL]>();
  using ft_type = std::remove_reference<decltype(forcing.m_ft)>::type;
  forcing.m_ft = ft_type(ft_in.data(),nelem);

  // Momentum forcing
  auto fm_in = get_dyn_field("FM_dyn").get_view<Homme::Scalar**[NP][NP][NVL]>();
  using fm_type = std::remove_reference<decltype(forcing.m_fm)>::type;
  forcing.m_fm = fm_type(fm_in.data(),nelem);

  // Homme has 3 components for FM, but the 3d (the omega forcing) is not computed
  // by EAMxx, so we set FM(3)=0 right away
  get_dyn_field("FM_dyn").get_component(2).deep_copy(0);

  if (params.do_3d_turbulence) {
    // SGS Eddy diffusivity for momentum
    auto Km_in = get_dyn_field("Km_dyn").get_view<Homme::Scalar*[NP][NP][NVL]>();
    using turb_type_mom = std::remove_reference<decltype(derived.m_turb_diff_mom)>::type;
    derived.m_turb_diff_mom = turb_type_mom(Km_in.data(), nelem);

    // SGS Eddy diffusivity for heat
    auto Kh_in = get_dyn_field("Kh_dyn").get_view<Homme::Scalar*[NP][NP][NVL]>();
    using turb_type_heat = std::remove_reference<decltype(derived.m_turb_diff_heat)>::type;
    derived.m_turb_diff_heat = turb_type_heat(Kh_in.data(), nelem);
  }
}

void HommeInstance::copy_dyn_states_to_all_timelevels ()
{
  const auto& c = Homme::Context::singleton();

  // State time slices
  const int nm1     = c.get<Homme::TimeLevel>().nm1;
  const int n0      = c.get<Homme::TimeLevel>().n0;
  const int np1     = c.get<Homme::TimeLevel>().np1;
  const int n0_qdp  = c.get<Homme::TimeLevel>().n0_qdp;
  const int np1_qdp = c.get<Homme::TimeLevel>().np1_qdp;

  // States
  const auto dp3d      = get_dyn_field("dp3d_dyn");
  const auto ps        = get_dyn_field("ps_dyn");
  const auto v         = get_dyn_field("v_dyn");
  const auto w_i       = get_dyn_field("w_int_dyn");
  const auto phinh_i   = get_dyn_field("phi_int_dyn");
  const auto vtheta_dp = get_dyn_field("vtheta_dp_dyn");
  const auto qdp       = get_dyn_field("Qdp_dyn");

  // Note: it might be somewhat faster to do a single parallel region,
  //       rather than 13 deep copies. However, this is much clearer
  //       to read, and since this method is called only during
  //       initialization, we prefer a more readable method
  dp3d.subfield(1,nm1).deep_copy(dp3d.subfield(1,n0));
  dp3d.subfield(1,np1).deep_copy(dp3d.subfield(1,n0));
  ps.subfield(1,nm1).deep_copy(ps.subfield(1,n0));
  ps.subfield(1,np1).deep_copy(ps.subfield(1,n0));
  v.subfield(1,nm1).deep_copy(v.subfield(1,n0));
  v.subfield(1,np1).deep_copy(v.subfield(1,n0));
  w_i.subfield(1,nm1).deep_copy(w_i.subfield(1,n0));
  w_i.subfield(1,np1).deep_copy(w_i.subfield(1,n0));
  phinh_i.subfield(1,nm1).deep_copy(phinh_i.subfield(1,n0));
  phinh_i.subfield(1,np1).deep_copy(phinh_i.subfield(1,n0));
  vtheta_dp.subfield(1,nm1).deep_copy(vtheta_dp.subfield(1,n0));
  vtheta_dp.subfield(1,np1).deep_copy(vtheta_dp.subfield(1,n0));
  qdp.subfield(1,np1_qdp).deep_copy(qdp.subfield(1,n0_qdp));
}

void HommeInstance::
init_dyn_state_from_gll_ic (AbstractRemapper& gll2dyn,
                            const Field& horiz_winds_gll,
                            const Field& T_mid_gll,
                            const Field& ps_gll,
                            const Field& phis_gll,
                            const Field& tracers_gll,
                            Field& pseudo_density_gll)
{
  // Some types
  using Pack         = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using KT           = KokkosTypes<DefaultDevice>;
  using ColOps       = ColumnOps<DefaultDevice,Real>;
  using PF           = PhysicsFunctions<DefaultDevice>;
  using TPF          = ekat::TeamPolicyFactory<KT::ExeSpace>;
  using EOS          = Homme::EquationOfState;
  using WorkspaceMgr = ekat::WorkspaceManager<Pack, DefaultDevice>;

  // Some Homme structures
  const auto& c = Homme::Context::singleton();
  auto& params = c.get<Homme::SimulationParams>();
  const auto& hvcoord = c.get<Homme::HybridVCoord>();
  const auto& tl = c.get<Homme::TimeLevel>();

  // Some extents
  const auto dp_ref = pseudo_density_gll.get_view<Real**>();
  const auto ncols = dp_ref.extent_int(0);
  const auto nlevs = dp_ref.extent_int(1);
  constexpr int NGP = HOMMEXX_NP;
  const int nelem = get_dyn_field("dp3d_dyn").get_header().get_identifier().get_layout().dim(0);
  const int qsize = params.qsize;
  const int npacks_mid = ekat::npack<Pack>(nlevs);
  const int npacks_int = ekat::npack<Pack>(nlevs+1);

  // Bootstrap dp on phys grid, and let the ic remapper transfer dp on dyn grid
  // NOTE: HybridVCoord already stores hyai and hybi deltas as packed views,
  //       but it used KokkosKernels packs, which are incompatible with ekat::Pack.
  //       If Homme switched to ekat::Pack, you could do the loop below with packs.
  const auto ps0 = hvcoord.ps0;
  const auto ps_ref = ps_gll.get_view<const Real*>();
  const auto hyai = hvcoord.hybrid_ai;
  const auto hybi = hvcoord.hybrid_bi;
  const auto policy_dp = TPF::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy_dp, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int icol = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,nlevs),
                        [&](const int ilev) {
       dp_ref(icol,ilev) = (hyai(ilev+1)-hyai(ilev))*ps0
                         + (hybi(ilev+1)-hybi(ilev))*ps_ref(icol);
    });
    team.team_barrier();
  });

  // Import IC from the GLL grid to the dyn grid
  // NOTE: the remapper supports (dynamic) subfields, so we remap directly into
  //       the n0 time slice of the states.
  gll2dyn.register_field(horiz_winds_gll,get_dyn_field("v_dyn").subfield(1,tl.n0,true));
  gll2dyn.register_field(pseudo_density_gll,get_dyn_field("dp3d_dyn").subfield(1,tl.n0,true));
  gll2dyn.register_field(ps_gll,get_dyn_field("ps_dyn").subfield(1,tl.n0,true));
  gll2dyn.register_field(phis_gll,get_dyn_field("phis_dyn"));
  gll2dyn.register_field(T_mid_gll,get_dyn_field("vtheta_dp_dyn").subfield(1,tl.n0,true));
  gll2dyn.register_field(tracers_gll,get_dyn_field("Q_dyn"));
  gll2dyn.registration_ends();
  gll2dyn.remap_fwd();

  // Wheter w_int is computed or not, Homme still does some global reduction on w_int when
  // printing the state, so we need to make sure it doesn't contain NaNs
  get_dyn_field("w_int_dyn").deep_copy(0.0);

  // Homme states
  const auto dp_view  = get_dyn_field("dp3d_dyn").get_view<Pack*****>();
  const auto vth_view = get_dyn_field("vtheta_dp_dyn").get_view<Pack*****>();
  const auto Q_view   = get_dyn_field("Q_dyn").get_view<Pack*****>();

  // State time slices
  const int n0  = tl.n0;
  const int n0_qdp  = tl.n0_qdp;

  const auto phis_dyn_view = get_dyn_field("phis_dyn").get_view<const Real***>();
  const auto phi_int_view = get_dyn_field("phi_int_dyn").get_view<Pack*****>();
  const auto hyai0 = hvcoord.hybrid_ai0;
  // Need two temporaries, for pi_mid and pi_int
  const auto policy = TPF::get_thread_range_parallel_scan_team_policy(nelem*NGP*NGP,npacks_mid);
  WorkspaceMgr wsm(npacks_int,2,policy);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int ie  =  team.league_rank() / (NGP*NGP);
    const int igp = (team.league_rank() / NGP) % NGP;
    const int jgp =  team.league_rank() % NGP;

    // Compute dp assuming hydrostatic initial state
    auto dp = ekat::subview(dp_view,ie,n0,igp,jgp);

    // Compute p_mid
    auto ws = wsm.get_workspace(team);
    ekat::Unmanaged<WorkspaceMgr::view_1d<Pack> > p_int, p_mid;
    ws.take_many_and_reset<2>({"p_int", "p_mid"}, {&p_int, &p_mid});

    ColOps::column_scan<true>(team,nlevs,dp,p_int,ps0*hyai0);
    team.team_barrier();
    ColOps::compute_midpoint_values(team,nlevs,p_int,p_mid);
    team.team_barrier();

    // Convert T->Theta->VTheta->VTheta*dp in place
    auto T      = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto vTh_dp = ekat::subview(vth_view,ie,n0,igp,jgp);
    auto qv     = ekat::subview(Q_view,ie,0,igp,jgp);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team,npacks_mid),
                         [&](const int ilev) {
      const auto th = PF::calculate_theta_from_T(T(ilev),p_mid(ilev));
      vTh_dp(ilev) = PF::calculate_virtual_temperature(th,qv(ilev))*dp(ilev);
    });
    team.team_barrier();

    // Init geopotential
    auto dphi   = [&](const int ilev)->Pack {
      return EOS::compute_dphi(vTh_dp(ilev), p_mid(ilev));
    };
    auto phi_int = ekat::subview(phi_int_view,ie,n0,igp,jgp);
    ColOps::column_scan<false>(team,nlevs,dphi,phi_int,phis_dyn_view(ie,igp,jgp));
  });

  // For initial runs, it's easier to prescribe IC for Q, and compute Qdp = Q*dp
  auto& tracers = c.get<Homme::Tracers>();
  const auto qdp = tracers.qdp;
  const auto q   = tracers.Q;
  const auto dp  = c.get<Homme::ElementsState>().m_dp3d;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,nelem*qsize*NGP*NGP*npacks_mid),
                       KOKKOS_LAMBDA (const int idx) {
    const int ie =  idx / (qsize*NGP*NGP*npacks_mid);
    const int iq = (idx / (NGP*NGP*npacks_mid)) % qsize;
    const int ip = (idx / (NGP*npacks_mid)) % NGP;
    const int jp = (idx / npacks_mid) % NGP;
    const int k  = idx % npacks_mid;

    qdp(ie,n0_qdp,iq,ip,jp,k) = q(ie,iq,ip,jp,k) * dp(ie,n0,ip,jp,k);
  });

  // If "instant" averaging type is used for output,
  // an initial output is performed before AD processes
  // are run. If omega_dyn output is requested, it will
  // not have valid computed values for this initial
  // output. Set to zero avoid potential FPE.
  get_dyn_field("omega_dyn").deep_copy(0);

  // Copy IC states on all timelevel slices
  copy_dyn_states_to_all_timelevels ();
}

} // namespace scream
