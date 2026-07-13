/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "HyperviscosityFunctorImpl.hpp"

#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "profiling.hpp"
#include "ColumnOps.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "mpi/Comm.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/Connectivity.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Homme
{

namespace {

constexpr int max_dynamic_sgs_subcycles = 12;

constexpr Real get_lambda_vis ()
{
  switch (NP) {
  case 2: return 12.0;
  case 3: return 30.0;
  case 4: return 91.6742;
  case 5: return 190.1176;
  case 6: return 374.7788;
  case 7: return 652.3015;
  default: return 0.0;
  }
}

KOKKOS_INLINE_FUNCTION
Real get_local_laplace_metric (const Real a, const Real b, const Real c, const Real d,
                               const Real lambda_vis, const Real scale_factor_inv)
{
  const Real s11 = a*a + c*c;
  const Real s22 = b*b + d*d;
  const Real s12 = a*b + c*d;
  const Real disc = (s11 - s22)*(s11 - s22) + 4.0*s12*s12;
  const Real max_eig = 0.5 * (s11 + s22 + std::sqrt(disc));
  const Real norm_dinv = std::sqrt(max_eig);
  return lambda_vis * (scale_factor_inv * norm_dinv) * (scale_factor_inv * norm_dinv);
}

} // anonymous namespace

HyperviscosityFunctorImpl::
HyperviscosityFunctorImpl (const SimulationParams&     params,
                           const ElementsGeometry&     geometry,
                           const ElementsState&        state,
                           const ElementsDerivedState& derived)
 : m_num_elems(state.num_elems())
 , m_data (params.hypervis_subcycle,params.hypervis_subcycle_sgs,params.hypervis_subcycle_tom,
		       params.nu_ratio1,params.nu_ratio2,params.nu_top,params.nu,
		       params.nu_p,params.nu_s,params.hypervis_scaling,
                       params.do_3d_turbulence)
 , m_state   (state)
 , m_derived (derived)
 , m_geometry (geometry)
 , m_sphere_ops (Context::singleton().get<SphereOperators>())
 , m_hvcoord (Context::singleton().get<HybridVCoord>())
 , m_policy_update_states (Homme::get_default_team_policy<ExecSpace,TagUpdateStates>(m_num_elems))
 , m_policy_first_laplace (Homme::get_default_team_policy<ExecSpace,TagFirstLaplaceHV>(m_num_elems))
 , m_policy_pre_exchange (Homme::get_default_team_policy<ExecSpace, TagHyperPreExchange>(m_num_elems))
 , m_policy_nutop_laplace (Homme::get_default_team_policy<ExecSpace, TagNutopLaplace>(m_num_elems))
 , m_policy_nutop_update_states (Homme::get_default_team_policy<ExecSpace,TagNutopUpdateStates>(m_num_elems))
 , m_policy_sgsturb_laplace (Homme::get_default_team_policy<ExecSpace, TagSGSTurbLaplace>(m_num_elems))
 , m_policy_sgsturb_update_states (Homme::get_default_team_policy<ExecSpace,TagSGSTurbUpdateStates>(m_num_elems))
 , m_tu(m_policy_update_states)
{
  init_params(params);

  // Make sure the sphere operators have buffers large enough to accommodate this functor's needs
  m_sphere_ops.allocate_buffers(m_tu);
}

HyperviscosityFunctorImpl::
HyperviscosityFunctorImpl (const int num_elems, const SimulationParams &params)
  : m_num_elems(num_elems)
  , m_data (params.hypervis_subcycle,params.hypervis_subcycle_sgs,params.hypervis_subcycle_tom,
		        params.nu_ratio1,params.nu_ratio2,params.nu_top,params.nu,
		        params.nu_p,params.nu_s,params.hypervis_scaling,
                        params.do_3d_turbulence)
  , m_hvcoord (Context::singleton().get<HybridVCoord>())
  , m_policy_update_states (Homme::get_default_team_policy<ExecSpace,TagUpdateStates>(m_num_elems))
  , m_policy_first_laplace (Homme::get_default_team_policy<ExecSpace,TagFirstLaplaceHV>(m_num_elems))
  , m_policy_pre_exchange (Homme::get_default_team_policy<ExecSpace, TagHyperPreExchange>(m_num_elems))
  , m_policy_nutop_laplace (Homme::get_default_team_policy<ExecSpace, TagNutopLaplace>(m_num_elems))
  , m_policy_nutop_update_states (Homme::get_default_team_policy<ExecSpace,TagNutopUpdateStates>(m_num_elems))
  , m_policy_sgsturb_laplace (Homme::get_default_team_policy<ExecSpace, TagSGSTurbLaplace>(m_num_elems))
  , m_policy_sgsturb_update_states (Homme::get_default_team_policy<ExecSpace,TagSGSTurbUpdateStates>(m_num_elems))
  , m_tu(m_policy_update_states)
{
  init_params(params);
}

void HyperviscosityFunctorImpl::init_params(const SimulationParams& params)
{
  // Sanity check
  assert(params.params_set);

  if (m_data.nu_top>0) {

    m_nu_scale_top = ExecViewManaged<Scalar[NUM_LEV]>("nu_scale_top");
    ExecViewManaged<Scalar[NUM_LEV]>::host_mirror_type h_nu_scale_top;
    h_nu_scale_top = Kokkos::create_mirror_view(m_nu_scale_top);

    const auto etai_h = Kokkos::create_mirror_view(m_hvcoord.etai);
    const auto etam_h = Kokkos::create_mirror_view(m_hvcoord.etam);
    Kokkos::deep_copy(etai_h, m_hvcoord.etai);
    Kokkos::deep_copy(etam_h, m_hvcoord.etam);

    for (int phys_lev=0; phys_lev < NUM_LEV*VECTOR_SIZE; ++phys_lev) {
      const int ilev = phys_lev / VECTOR_SIZE;
      const int ivec = phys_lev % VECTOR_SIZE;

      Real ptop_over_press;

      //prevent padding of nu_scale to get nans to avoid last interface levels
      //of w, phi to be nans, too
      if( phys_lev < NUM_PHYSICAL_LEV ){
        //etai is num_interface_lev, that is, 129 or 73
        //etam is num_lev, so packs
        if ( etai_h(0) == 0.0) {
          ptop_over_press = etam_h(0)[0] / etam_h(ilev)[ivec];
        }else{
          ptop_over_press = etai_h(0) / etam_h(ilev)[ivec];
        }
      }else{
          ptop_over_press = 0.0;
      }

      auto val = 16.0*ptop_over_press*ptop_over_press / (ptop_over_press*ptop_over_press + 1);
      if ( val < 0.15 ) val = 0.0;
      h_nu_scale_top(ilev)[ivec] = val;

      // This is the equivalent of nlev_tom in the F90 code.
      if (val != 0) m_nu_scale_top_ilev_pack_lim = phys_lev + 1;
    }
    Kokkos::deep_copy(m_nu_scale_top, h_nu_scale_top);

    // Convert to pack index.
    m_nu_scale_top_ilev_pack_lim = ((m_nu_scale_top_ilev_pack_lim + VECTOR_SIZE - 1) /
                                    VECTOR_SIZE);
  }

  // Init ElementOps
  m_elem_ops.init(m_hvcoord);

  // Init Equation of state
  m_eos.init(params.theta_hydrostatic_mode,m_hvcoord);

#ifdef HOMMEXX_BFB_TESTING
  m_process_nh_vars = 1;
#else
  m_process_nh_vars = not params.theta_hydrostatic_mode;
#endif
}

void HyperviscosityFunctorImpl::setup(const ElementsGeometry&     geometry,
                                      const ElementsState&        state,
                                      const ElementsDerivedState& derived)
{
  m_state = state;
  assert(m_num_elems == m_state.num_elems()); //Sanity check
  m_derived = derived;
  m_geometry = geometry;
  m_sphere_ops = Context::singleton().get<SphereOperators>();

  // Make sure the sphere operators have buffers large enough to accommodate this functor's needs
  m_sphere_ops.allocate_buffers(m_tu);
}

int HyperviscosityFunctorImpl::requested_buffer_size () const {
  constexpr int size_mid_scalar =   NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar =   NP*NP*NUM_LEV_P*VECTOR_SIZE;

  // Number of scalar/vector int/mid buffers needed, with size nelems
  const int mid_vectors_nelems = 1;
  const int int_scalars_nelems = 0 + (m_process_nh_vars ? 2 : 0);
  const int mid_scalars_nelems = 2 + (m_process_nh_vars ? 2 : 0);

  const int size = m_num_elems*(mid_scalars_nelems*size_mid_scalar +
                                mid_vectors_nelems*size_mid_vector +
                                int_scalars_nelems*size_int_scalar);

  return size;
}

void HyperviscosityFunctorImpl::init_buffers (const FunctorsBuffersManager& fbm) {
  Errors::runtime_check(fbm.allocated_size()>=requested_buffer_size(), "Error! Buffers size not sufficient.\n");

  constexpr int size_mid_scalar =   NP*NP*NUM_LEV;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV;
  constexpr int size_int_scalar =   NP*NP*NUM_LEV_P;

  auto mem_in = fbm.get_memory();
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
  const int nelems = m_geometry.num_elems();

  // Tens quantities (persistent views => nelems)
  m_buffers.dptens = decltype(m_buffers.dptens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  m_buffers.ttens = decltype(m_buffers.ttens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  if (m_process_nh_vars) {
    m_buffers.wtens = decltype(m_buffers.wtens)(mem,nelems);
    mem += size_mid_scalar*nelems;

    m_buffers.phitens = decltype(m_buffers.phitens)(mem,nelems);
    mem += size_mid_scalar*nelems;

    m_buffers.turb_diff_heat_i = decltype(m_buffers.turb_diff_heat_i)(mem,nelems);
    mem += size_int_scalar*nelems;

    m_buffers.turb_diff_mom_i = decltype(m_buffers.turb_diff_mom_i)(mem,nelems);
    mem += size_int_scalar*nelems;
  }

  m_buffers.vtens = decltype(m_buffers.vtens)(mem,nelems);
  mem += size_mid_vector*nelems;

  const int used_mem = reinterpret_cast<Real*>(mem)-mem_in;
  if (used_mem < requested_buffer_size()) {
    printf("[HyperviscosityFunctorImpl] Warning! We used less memory than we said we would: %d instead of %d\n",
           used_mem, requested_buffer_size());
  } else if (used_mem > requested_buffer_size()) {
    std::string msg = "[HyperviscosityFunctorImpl] Error! We used more memory than we said we would: "
                    + std::to_string(used_mem) + " instead of "
                    + std::to_string(requested_buffer_size()) + "\n";
    Errors::runtime_abort(msg);
  }
}

void HyperviscosityFunctorImpl::init_boundary_exchanges () {
  auto bm_exchange = Context::singleton().get<MpiBuffersManagerMap>()[MPI_EXCHANGE];
  const auto& sp = Context::singleton().get<SimulationParams>();
  m_be = std::make_shared<BoundaryExchange>();
  m_be_tom = std::make_shared<BoundaryExchange>();
  m_be_sgs = std::make_shared<BoundaryExchange>();
  m_be->set_label("Hyperviscosity-std");
  m_be_tom->set_label("Hyperviscosity-TOM");
  m_be_sgs->set_label("Hyperviscosity-SGS");
  std::shared_ptr<BoundaryExchange> bes[] = {m_be, m_be_tom, m_be_sgs};
  const int nlevs[] = {NUM_LEV, m_nu_scale_top_ilev_pack_lim, NUM_LEV};
  for (int i = 0; i < 3; ++i) {
    if (i == 1 && m_data.nu_top <= 0) continue;
    auto be = bes[i];
    be->set_diagnostics_level(sp.internal_diagnostics_level);
    const auto nlev = nlevs[i];
    be->set_buffers_manager(bm_exchange);
    if (m_process_nh_vars) {
      be->set_num_fields(0, 0, 6);
    } else {
      be->set_num_fields(0, 0, 4);
    }
    be->register_field(m_buffers.dptens, nlev);
    be->register_field(m_buffers.ttens, nlev);
    if (m_process_nh_vars) {
      be->register_field(m_buffers.wtens, nlev);
      be->register_field(m_buffers.phitens, nlev);
    }
    be->register_field(m_buffers.vtens, 2, 0, nlev);
    be->registration_completed();
  }
}//initBE

void HyperviscosityFunctorImpl::run (const int np1, const Real dt, const Real eta_ave_w)
{
  m_data.np1 = np1;

  m_data.dt = dt;
  if (m_data.hypervis_subcycle > 0) { 
    m_data.dt_hvs = dt/m_data.hypervis_subcycle;
  }else{
    //won't be used
    m_data.dt_hvs = -1.0;
  }
  if (m_data.hypervis_subcycle_tom > 0) { 
    m_data.dt_hvs_tom = dt/m_data.hypervis_subcycle_tom;
  }else{
    //won't be used
    m_data.dt_hvs_tom = -1.0;
  }
  m_data.hypervis_subcycle_sgs_eff = compute_sgs_subcycle_count(dt);
  if (m_data.hypervis_subcycle_sgs_eff > 0) {
    m_data.dt_hvs_sgs = dt/m_data.hypervis_subcycle_sgs_eff;
  } else {
    m_data.dt_hvs_sgs = -1.0;
  }
  if (m_data.hypervis_subcycle_sgs > 0) {
    clip_sgs_diffusivities_for_fixed_subcycling();
  }
  if (m_data.hypervis_subcycle_sgs_eff > 1 && Context::singleton().get<Comm>().root()) {
    if (m_data.hypervis_subcycle_sgs > 0) {
      printf("Warning: SGS horizontal diffusion is using fixed subcycling of %d for this run.\n",
             m_data.hypervis_subcycle_sgs_eff);
    } else if (m_data.hypervis_subcycle_sgs_eff == max_dynamic_sgs_subcycles) {
      printf("Warning: SGS horizontal diffusion reached the maximum allowed dynamic subcycling of %d for this step.\n",
             m_data.hypervis_subcycle_sgs_eff);
    } else {
      printf("Warning: SGS horizontal diffusion increased dynamic subcycling to %d for this step.\n",
             m_data.hypervis_subcycle_sgs_eff);
    }
  }
  m_data.eta_ave_w = eta_ave_w;

  // Convert vtheta_dp -> theta
  auto state = m_state;
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(state.num_elems()),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    const int ie = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // theta->vtheta
      auto vtheta = Homme::subview(state.m_vtheta_dp,ie,np1,igp,jgp);
      auto dp = Homme::subview(state.m_dp3d,ie,np1,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) /= dp(ilev);
      });
    });
  });
  Kokkos::fence();

  for (int icycle = 0; icycle < m_data.hypervis_subcycle; ++icycle) {
    GPTLstart("hvf-bhwk");
    biharmonic_wk_theta ();
    GPTLstop("hvf-bhwk");

    Kokkos::parallel_for(m_policy_pre_exchange, *this);
    Kokkos::fence();

    // Exchange
    assert (m_be->is_registration_completed());
    GPTLstart("hvf-bexch");
    m_be->exchange();
    GPTLstop("hvf-bexch");

    // Update states
    Kokkos::parallel_for(m_policy_update_states, *this);
    Kokkos::fence();

  } //subcycle

  // SGS Horizontal turbulent diffusion
  if (m_data.do_3d_turbulence > 0) {
    for (int icycle = 0; icycle < m_data.hypervis_subcycle_sgs_eff; ++icycle) {
      // laplace(fields) --> ttens, etc.
      Kokkos::parallel_for(m_policy_sgsturb_laplace, *this);
      Kokkos::fence();

      // exchange is done on ttens, dptens, vtens, etc.
      assert (m_be_sgs->is_registration_completed());
      GPTLstart("sgsturb-bexch");
      m_be_sgs->exchange();
      GPTLstop("sgsturb-bexch");

      // update states
      Kokkos::parallel_for(m_policy_sgsturb_update_states, *this);
      Kokkos::fence();
    }
  } // SGS horizontal turbulent diffusion

  // Convert theta back to vtheta, and adjust w at surface
  auto geo = m_geometry;
  auto process_nh_vars = m_process_nh_vars;
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(state.num_elems()),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    const int ie = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // theta->vtheta
      auto vtheta = Homme::subview(state.m_vtheta_dp,ie,np1,igp,jgp);
      auto dp = Homme::subview(state.m_dp3d,ie,np1,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) *= dp(ilev);
      });

      // Fix w at surface:
      // Adjust w_i at the surface, since velocity has changed
      if (process_nh_vars) {
        Kokkos::single(Kokkos::PerThread(team),[&](){
          using InfoI = ColInfo<NUM_INTERFACE_LEV>;
          using InfoM = ColInfo<NUM_PHYSICAL_LEV>;
          constexpr int LAST_MID_PACK     = InfoM::LastPack;
          constexpr int LAST_MID_PACK_END = InfoM::LastPackEnd;
          constexpr int LAST_INT_PACK     = InfoI::LastPack;
          constexpr int LAST_INT_PACK_END = InfoI::LastPackEnd;
          constexpr Real g = PhysicalConstants::g;

          const auto& grad_x = geo.m_gradphis(ie,0,igp,jgp);
          const auto& grad_y = geo.m_gradphis(ie,1,igp,jgp);
          const auto& u = state.m_v(ie,np1,0,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END];
          const auto& v = state.m_v(ie,np1,1,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END];

          auto& w = state.m_w_i(ie,np1,igp,jgp,LAST_INT_PACK)[LAST_INT_PACK_END];

          w = (u*grad_x+v*grad_y) / g;
        });
      }
    });
  });//conversion back to vtheta

  Kokkos::fence();

  // sponge layer 
  if (m_data.nu_top > 0) {
    for (int icycle = 0; icycle < m_data.hypervis_subcycle_tom; ++icycle) {
      // laplace(fields) --> ttens, etc.
      Kokkos::parallel_for(m_policy_nutop_laplace, *this);
      Kokkos::fence();

      // exchange is done on ttens, dptens, vtens, etc.
      assert (m_be->is_registration_completed());
      GPTLstart("hvf-bexch");
      m_be_tom->exchange();
      GPTLstop("hvf-bexch");

      Kokkos::parallel_for(m_policy_nutop_update_states, *this);
      Kokkos::fence();
    }
  } // for sponge layer
} // run()

int HyperviscosityFunctorImpl::compute_sgs_subcycle_count (const Real dt) const
{
  // Semantics:
  //   -1: inherit hypervis_subcycle, with no dynamic SGS adaptation
  //    0: dynamically choose SGS subcycling, starting from a baseline of 1
  //   >0: fixed SGS subcycling for the full run, with no dynamic adaptation
  if (m_data.hypervis_subcycle_sgs > 0) {
    return m_data.hypervis_subcycle_sgs;
  }

  int nsub = 1;

  if (not m_data.do_3d_turbulence) {
    return nsub;
  }

  const Real lambda_vis = get_lambda_vis();
  if (lambda_vis <= 0) {
    return nsub;
  }

  const auto dinv = m_geometry.m_dinv;
  const auto Km = m_derived.m_turb_diff_mom;
  const auto Kh = m_derived.m_turb_diff_heat;
  const Real scale_factor_inv = 1.0 / m_geometry.m_scale_factor;

  // Estimate the largest explicit diffusive CFL over all local GLL points and
  // vertical levels using the local SGS diffusivities and the local element
  // metric. This mirrors HOMME's existing Laplacian stability estimate.
  Real local_max_cfl = 0.0;
  Kokkos::parallel_reduce(
      "compute_sgs_max_cfl",
      Kokkos::RangePolicy<ExecSpace>(0, m_num_elems * NP * NP),
      KOKKOS_LAMBDA (const int idx, Real& thread_max) {
        const int ie = idx / (NP * NP);
        const int rem = idx % (NP * NP);
        const int igp = rem / NP;
        const int jgp = rem % NP;

        const Real a = dinv(ie,0,0,igp,jgp);
        const Real b = dinv(ie,0,1,igp,jgp);
        const Real c = dinv(ie,1,0,igp,jgp);
        const Real d = dinv(ie,1,1,igp,jgp);

        // Recover a local 2-norm of D^{-1} from the largest eigenvalue of
        // D^{-1} (D^{-1})^T, then convert it into the Laplacian stability
        // factor used in HOMME's viscous CFL estimate.
        const Real laplace_metric = get_local_laplace_metric(a, b, c, d, lambda_vis, scale_factor_inv);

        // Use the largest of Km and Kh at this horizontal point as a scalar
        // bound for the explicit SGS diffusion step.
        Real point_max_k = 0.0;
        for (int k = 0; k < NUM_LEV; ++k) {
          const auto km = Km(ie,igp,jgp,k);
          const auto kh = Kh(ie,igp,jgp,k);
          for (int s = 0; s < VECTOR_SIZE; ++s) {
            const int phys_lev = k * VECTOR_SIZE + s;
            if (phys_lev < NUM_PHYSICAL_LEV) {
              if (km[s] > point_max_k) point_max_k = km[s];
              if (kh[s] > point_max_k) point_max_k = kh[s];
            }
          }
        }

        // Forward-Euler diffusion is stable for CFL <= 1 in this normalized
        // estimate, so values above 1 imply we need more SGS substeps.
        const Real point_cfl = 0.5 * dt * laplace_metric * point_max_k;
        if (point_cfl > thread_max) thread_max = point_cfl;
      },
      Kokkos::Max<Real>(local_max_cfl));

  // Promote the local maximum to a global one so every rank uses the same SGS
  // subcycle count for this dynamics step.
  Real global_max_cfl = local_max_cfl;
  const auto& comm = Context::singleton().get<Comm>();
  MPI_Allreduce(&local_max_cfl, &global_max_cfl, 1, MPI_DOUBLE, MPI_MAX, comm.mpi_comm());

  // If the estimated CFL is, e.g., 2.3, take 3 SGS substeps so the effective
  // per-substep CFL is brought back below 1. Cap the adaptive path to avoid
  // runaway cost in pathological cases.
  nsub = std::max(nsub, static_cast<int>(std::ceil(global_max_cfl)));
  const int nsub_uncapped = nsub;
  nsub = std::min(nsub, max_dynamic_sgs_subcycles);

  return nsub;
}

void HyperviscosityFunctorImpl::clip_sgs_diffusivities_for_fixed_subcycling () const
{
  const Real lambda_vis = get_lambda_vis();
  if (not m_data.do_3d_turbulence || lambda_vis <= 0 || m_data.dt_hvs_sgs <= 0) {
    return;
  }

  const auto dinv = m_geometry.m_dinv;
  const auto Km = m_derived.m_turb_diff_mom;
  const auto Kh = m_derived.m_turb_diff_heat;
  const Real scale_factor_inv = 1.0 / m_geometry.m_scale_factor;

  Kokkos::parallel_for(
      "clip_sgs_diffusivities_fixed_subcycling",
      Kokkos::RangePolicy<ExecSpace>(0, m_num_elems * NP * NP),
      KOKKOS_LAMBDA (const int idx) {
        const int ie = idx / (NP * NP);
        const int rem = idx % (NP * NP);
        const int igp = rem / NP;
        const int jgp = rem % NP;

        const Real a = dinv(ie,0,0,igp,jgp);
        const Real b = dinv(ie,0,1,igp,jgp);
        const Real c = dinv(ie,1,0,igp,jgp);
        const Real d = dinv(ie,1,1,igp,jgp);
        const Real laplace_metric = get_local_laplace_metric(a, b, c, d, lambda_vis, scale_factor_inv);

        if (laplace_metric <= 0) {
          return;
        }

        const Real max_diffusivity = 2.0 / (m_data.dt_hvs_sgs * laplace_metric);
        for (int k = 0; k < NUM_LEV; ++k) {
          auto km = Km(ie,igp,jgp,k);
          auto kh = Kh(ie,igp,jgp,k);
          for (int s = 0; s < VECTOR_SIZE; ++s) {
            const int phys_lev = k * VECTOR_SIZE + s;
            if (phys_lev < NUM_PHYSICAL_LEV) {
              if (km[s] > max_diffusivity) km[s] = max_diffusivity;
              if (kh[s] > max_diffusivity) kh[s] = max_diffusivity;
            }
          }
          Km(ie,igp,jgp,k) = km;
          Kh(ie,igp,jgp,k) = kh;
        }
      });
  Kokkos::fence();
}

void HyperviscosityFunctorImpl::biharmonic_wk_theta() const
{
  // For the first laplacian we use a differnt kernel, which uses directly the states
  // at timelevel np1 as inputs, and subtracts the reference states.
  // This way we avoid copying the states to *tens buffers.
  Kokkos::parallel_for(m_policy_first_laplace, *this);
  Kokkos::fence();

  // Exchange
  assert (m_be->is_registration_completed());
  GPTLstart("hvf-bexch");
  m_be->exchange(m_geometry.m_rspheremp);
  GPTLstop("hvf-bexch");

  // Compute second laplacian, tensor or const hv
  const int ne = m_geometry.num_elems();
  if ( m_data.consthv ) {
    auto policy = Homme::get_default_team_policy<ExecSpace,TagSecondLaplaceConstHV>(ne);
    Kokkos::parallel_for(policy, *this);
  }else{
    auto policy = Homme::get_default_team_policy<ExecSpace,TagSecondLaplaceTensorHV>(ne);
    Kokkos::parallel_for(policy, *this);
  }
  Kokkos::fence();
} //biharmonic

// Laplace for nu_top
KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagNutopLaplace&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));

  // Laplacian of layer thickness
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.dptens,kv.ie),
                              m_nu_scale_top_ilev_pack_lim);
  // Laplacian of theta
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.ttens,kv.ie),
                              m_nu_scale_top_ilev_pack_lim);

  if (m_process_nh_vars) {
    // Laplacian of vertical velocity (do not compute last interface)
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_w_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.wtens,kv.ie),
                                                   m_nu_scale_top_ilev_pack_lim);
    // Laplacian of geopotential (do not compute last interface)
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.phitens,kv.ie),
                                                   m_nu_scale_top_ilev_pack_lim);
  }

  // Laplacian of velocity
  m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
                                         Homme::subview(m_state.m_v,kv.ie,m_data.np1),
                                         Homme::subview(m_buffers.vtens,kv.ie),
                                         m_nu_scale_top_ilev_pack_lim);

  kv.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(kv.team,NP*NP),
    [&] (const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      const auto utens  = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
      const auto vtens  = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
      const auto ttens  = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
      const auto dptens = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
     
      MidColumn wtens, phitens;
      if (m_process_nh_vars) {
        wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
        phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
      }

      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(kv.team, m_nu_scale_top_ilev_pack_lim),
        [&] (const int ilev) {
          
          const auto xf = m_data.dt_hvs_tom  * m_nu_scale_top(ilev) * m_data.nu_top;
          utens(ilev)  *= xf;
          vtens(ilev)  *= xf;
          ttens(ilev)  *= xf;
          dptens(ilev) *= xf;

          if (m_process_nh_vars) {
            wtens(ilev)   *= xf;
            phitens(ilev) *= xf;
          }

        }); // threadvectorrange
    }); // teamthreadrange
} // TagNutopLaplace

KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagNutopUpdateStates&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));
  using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                       [&](const int idx) {
    const int igp = idx / NP;
    const int jgp = idx % NP;

    // Add Xtens quantities back to the states, except for vtheta
    auto u = Homme::subview(m_state.m_v,kv.ie,m_data.np1,0,igp,jgp);
    auto v = Homme::subview(m_state.m_v,kv.ie,m_data.np1,1,igp,jgp);
    auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
    auto dp     = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);

    auto utens   = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
    auto vtens   = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
    auto ttens   = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
    auto dptens  = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
    const auto& rspheremp = m_geometry.m_rspheremp(kv.ie,igp,jgp);

    MidColumn wtens, phitens;
    IntColumn w, phi_i;

    if (m_process_nh_vars) {
      wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
      phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
      w       = Homme::subview(m_state.m_w_i,kv.ie,m_data.np1,igp,jgp);
      phi_i   = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
    }

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, m_nu_scale_top_ilev_pack_lim),
                         [&](const int ilev) {
      utens(ilev)   *= rspheremp;
      vtens(ilev)   *= rspheremp;
      ttens(ilev)   *= rspheremp;
      dptens(ilev)  *= rspheremp;
      u(ilev)      += utens(ilev);
      v(ilev)      += vtens(ilev);
      vtheta(ilev) += ttens(ilev);
      dp(ilev)     += dptens(ilev);

      if (m_process_nh_vars) {
        wtens(ilev)   *= rspheremp;
        phitens(ilev) *= rspheremp;
        w(ilev)     += wtens(ilev);
        phi_i(ilev) += phitens(ilev);
      }
    }); // threadvectorrange
  }); // threadteamrange
} // tagUpdateStates2

// Laplace for horizontal SGS turbulent diffusion
KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagSGSTurbLaplace&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));
  using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                       [&](const int idx) {
    const int igp = idx / NP;
    const int jgp = idx % NP;

    auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
    auto dp = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);
    auto theta_ref = Homme::subview(m_state.m_ref_states.theta_ref,kv.ie,igp,jgp);
    auto dp_ref = Homme::subview(m_state.m_ref_states.dp_ref,kv.ie,igp,jgp);

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      vtheta(ilev) -= theta_ref(ilev);
      dp(ilev) -= dp_ref(ilev);
    });
  });

  kv.team_barrier();

  if (m_process_nh_vars) {
    // Diffuse only the perturbational geopotential, not the terrain-following
    // reference profile tied to phis.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto phi_i = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
      auto phi_i_ref = Homme::subview(m_state.m_ref_states.phi_i_ref,kv.ie,igp,jgp);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        phi_i(ilev) -= phi_i_ref(ilev);
      });

#ifndef XX_NONBFB_COMING
      if (NUM_LEV!=NUM_LEV_P) {
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          phi_i(NUM_LEV_P-1) -= phi_i_ref(NUM_LEV_P-1);
        });
      }
#endif
    });

    kv.team_barrier();
  }

  // Laplacian of layer thickness
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.dptens,kv.ie));
  // Laplacian of theta
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.ttens,kv.ie));
  if (m_process_nh_vars) {
    // Laplacian of vertical velocity
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_w_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.wtens,kv.ie));
    // Laplacian of geopotential
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.phitens,kv.ie));
  }

  // Laplacian of velocity
  m_sphere_ops.vlaplace_sphere_wk_contra(kv,
                                         1.0, // no nu_ratio here, we want plain Lap(v)
                                         Homme::subview(m_state.m_v,kv.ie,m_data.np1),
                                         Homme::subview(m_buffers.vtens,kv.ie));

  kv.team_barrier();

  if (m_process_nh_vars) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto phi_i = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
      auto phi_i_ref = Homme::subview(m_state.m_ref_states.phi_i_ref,kv.ie,igp,jgp);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        phi_i(ilev) += phi_i_ref(ilev);
      });

#ifndef XX_NONBFB_COMING
      if (NUM_LEV!=NUM_LEV_P) {
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          phi_i(NUM_LEV_P-1) += phi_i_ref(NUM_LEV_P-1);
        });
      }
#endif
    });

    kv.team_barrier();
  }

  Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                       [&](const int idx) {
    const int igp = idx / NP;
    const int jgp = idx % NP;

    auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
    auto dp = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);
    auto theta_ref = Homme::subview(m_state.m_ref_states.theta_ref,kv.ie,igp,jgp);
    auto dp_ref = Homme::subview(m_state.m_ref_states.dp_ref,kv.ie,igp,jgp);

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      vtheta(ilev) += theta_ref(ilev);
      dp(ilev) += dp_ref(ilev);
    });
  });

  kv.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(kv.team,NP*NP),
    [&] (const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      const auto utens  = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
      const auto vtens  = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
      const auto ttens  = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
      const auto dptens = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);

      const auto Km = Homme::subview(m_derived.m_turb_diff_mom,kv.ie,igp,jgp);
      const auto Kh = Homme::subview(m_derived.m_turb_diff_heat,kv.ie,igp,jgp);

      MidColumn wtens, phitens;
      IntColumn Km_i, Kh_i;

      if (m_process_nh_vars) {
        wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
        phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);

        // Diffusivities on the interface grid
        Km_i = Homme::subview(m_buffers.turb_diff_mom_i,kv.ie,igp,jgp);
        Kh_i = Homme::subview(m_buffers.turb_diff_heat_i,kv.ie,igp,jgp);

        // Get diffusivities on the interface vertical grid from those
        //  on the mid-point grid after any fixed-subcycle clipping
        ColumnOps::compute_interface_values(kv, Km, Km_i);
        ColumnOps::compute_interface_values(kv, Kh, Kh_i);
      }

      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
        [&] (const int k) {

          const auto xf_m = m_data.dt_hvs_sgs * Km(k); // Momentum diffusivity
          const auto xf_h = m_data.dt_hvs_sgs * Kh(k); // Heat diffusivity
          utens(k)  *= xf_m;
          vtens(k)  *= xf_m;
          ttens(k)  *= xf_h;
          dptens(k) *= xf_h;

          if (m_process_nh_vars) {
            const auto xf_mi = m_data.dt_hvs_sgs * Km_i(k); // Momentum diffusivity on interface
            const auto xf_hi = m_data.dt_hvs_sgs * Kh_i(k); // Heat diffusivity on interface
            wtens(k)   *= xf_mi;
            phitens(k) *= xf_hi;
          }

        }); // threadvectorrange
    }); // teamthreadrange
} // TagSGSTurbLaplace

// SGS Horizontal turbulent diffusion, update states
KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagSGSTurbUpdateStates&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));
  using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                       [&](const int idx) {
    const int igp = idx / NP;
    const int jgp = idx % NP;

    // Add Xtens quantities back to the states, except for vtheta
    auto u = Homme::subview(m_state.m_v,kv.ie,m_data.np1,0,igp,jgp);
    auto v = Homme::subview(m_state.m_v,kv.ie,m_data.np1,1,igp,jgp);
    auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
    auto dp     = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);

    auto utens   = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
    auto vtens   = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
    auto ttens   = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
    auto dptens  = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
    const auto& rspheremp = m_geometry.m_rspheremp(kv.ie,igp,jgp);

    MidColumn wtens, phitens;
    IntColumn w, phi_i;

    if (m_process_nh_vars) {
      wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
      phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
      w       = Homme::subview(m_state.m_w_i,kv.ie,m_data.np1,igp,jgp);
      phi_i   = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
    }

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                         [&](const int k) {
      utens(k)   *= rspheremp;
      vtens(k)   *= rspheremp;
      ttens(k)   *= rspheremp;
      dptens(k)  *= rspheremp;
      u(k)      += utens(k);
      v(k)      += vtens(k);
      vtheta(k) += ttens(k);
      dp(k)     += dptens(k);

      if (m_process_nh_vars) {
        wtens(k)   *= rspheremp;
        phitens(k) *= rspheremp;
        w(k)     += wtens(k);
        phi_i(k) += phitens(k);
      }

    }); // threadvectorrange
  }); // threadteamrange
} // tagSGSTurbUpdateStates

} // namespace Homme
