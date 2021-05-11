/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_REMAP_STATE_PROVIDER_HPP
#define HOMMEXX_REMAP_STATE_PROVIDER_HPP

#include "Elements.hpp"
#include "EquationOfState.hpp"
#include "ElementOps.hpp"
#include "ErrorDefs.hpp"
#include "ColumnOps.hpp"
#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "SimulationParams.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "Types.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme {
namespace Remap {

struct RemapStateProvider {

  EquationOfState   m_eos;
  ElementOps        m_elem_ops;
  ElementsState     m_state;
  ElementsGeometry  m_geometry;
  HybridVCoord      m_hvcoord;
  bool              m_process_nh_vars;

  // These two morally are d(w_i)/ds and d(phinh_i)/ds.
  // However, since in the remap we need to multiply by ds
  // (the layer thickness, aka dp), we simply compute
  // d(w_i) and d(phinh_i).
  ExecViewManaged<Scalar*  [NP][NP][NUM_LEV]> m_delta_w;
  ExecViewManaged<Scalar*  [NP][NP][NUM_LEV]> m_delta_phinh;

  ExecViewUnmanaged<Scalar*  [NP][NP][NUM_LEV  ]> m_temp;
  ExecViewUnmanaged<Scalar*  [NP][NP][NUM_LEV_P]> m_phi_ref;

  explicit RemapStateProvider(const Elements& elements)
   : m_state(elements.m_state)
   , m_geometry(elements.m_geometry)
  {
    // Fetch SimulationParams and HybridVCoord from the context
    const auto& params = Context::singleton().get<SimulationParams>();
    assert (params.params_set);

#ifdef HOMMEXX_BFB_TESTING
    m_process_nh_vars = true;
#else
    m_process_nh_vars = !params.theta_hydrostatic_mode;
#endif

    if (m_process_nh_vars) {
      m_delta_w     = decltype(m_delta_w) ("w_i increments",elements.num_elems());
      m_delta_phinh = decltype(m_delta_phinh) ("phinh_i increments",elements.num_elems());
    }

    m_hvcoord = Context::singleton().get<HybridVCoord>();
    assert (m_hvcoord.m_inited);

    m_eos.init(params.theta_hydrostatic_mode,m_hvcoord);
    m_elem_ops.init(m_hvcoord);
  }

  int requested_buffer_size (int num_teams) const {
    if (!m_process_nh_vars) {
      return 0;
    }

    using temp_type = decltype(m_temp);
    using phi_type = decltype(m_phi_ref);
    const int temp_size = temp_type::shmem_size(num_teams)/sizeof(Real);
    const int phi_size  = phi_type::shmem_size(num_teams)/sizeof(Real);
    return temp_size + phi_size;
  }

  void init_buffers(const FunctorsBuffersManager& fbm, int num_teams) {
    if (!m_process_nh_vars) {
      return;
    }

    Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

    m_temp = decltype(m_temp)(mem,num_teams);
    mem += m_temp.size();

    m_phi_ref = decltype(m_phi_ref)(mem,num_teams);
  }

  KOKKOS_INLINE_FUNCTION
  int num_states_remap() const {
    return (m_process_nh_vars ? 5 : 3);
  }

  KOKKOS_INLINE_FUNCTION
  int num_states_preprocess() const {
    return (m_process_nh_vars ? 2 : 0);
  }

  KOKKOS_INLINE_FUNCTION
  int num_states_postprocess() const {
    return (m_process_nh_vars ? 2 : 0);
  }

  KOKKOS_INLINE_FUNCTION
  bool is_intrinsic_state (const int istate) const {
    assert (istate>=0 && istate<num_states_remap());

    if (istate==0 || istate==1) {
      // Horizontal velocity needs to be rescaled by dp
      return true;
    }

    // Other quantities are already scaled by dp
    return false;
  }

  KOKKOS_INLINE_FUNCTION
  void preprocess_state (const KernelVariables& kv,
                         const int istate,
                         const int np1,
                         ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> dp) const {
    assert (m_process_nh_vars);
    // Note: w/phi are the state 3/4, but when it comes to pre/post processing they are 0/1
    if (istate==0) {
      // Compute delta_w
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;
        auto w_i = Homme::subview(m_state.m_w_i,kv.ie,np1,igp,jgp);
        auto delta_w = Homme::subview(m_delta_w,kv.ie,igp,jgp);

        ColumnOps::compute_midpoint_delta(kv,w_i,delta_w);
      });
    } else if (istate==1) {
      // Compute delta_phinh
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        // Subviews (note: recycle phi_ref for the temporary p_i)
        auto phinh_i     = Homme::subview(m_state.m_phinh_i,kv.ie,np1,igp,jgp);
        auto delta_phinh = Homme::subview(m_delta_phinh,kv.ie,igp,jgp);
        auto phi_ref     = Homme::subview(m_phi_ref,kv.team_idx,igp,jgp);
        auto dp_pt       = Homme::subview(dp,igp,jgp);
        auto p_i         = Homme::subview(m_phi_ref,kv.team_idx,igp,jgp);
        auto p           = Homme::subview(m_temp,kv.team_idx,igp,jgp);

        // Compute hydrostatic phi
        m_elem_ops.compute_hydrostatic_p(kv,dp_pt,p_i,p);
        m_eos.compute_phi_i(kv,m_geometry.m_phis(kv.ie,igp,jgp),
                            Homme::subview(m_state.m_vtheta_dp,kv.ie,np1,igp,jgp),
                            p,
                            phi_ref);

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV_P),
                             [&](const int ilev) {
          phinh_i(ilev) -= phi_ref(ilev);
        });

        // Now compute phinh_i increments
        ColumnOps::compute_midpoint_delta(kv,phinh_i,delta_phinh);
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void postprocess_state (const KernelVariables& kv,
                          const int istate,
                          const int np1,
                          ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> dp) const {
    assert (m_process_nh_vars);
    using InfoI = ColInfo<NUM_INTERFACE_LEV>;
    using InfoM = ColInfo<NUM_PHYSICAL_LEV>;
    // Note: w/phi are the state 3/4, but when it comes to pre/post processing they are 0/1
    if (istate==0) {
      // Update w_i
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;
        auto w_i = Homme::subview(m_state.m_w_i,kv.ie,np1,igp,jgp);
        auto delta_w = Homme::subview(m_delta_w,kv.ie,igp,jgp);

        auto minus_delta = [&](const int ilev)->Scalar { return -delta_w(ilev); };
        ColumnOps::column_scan_mid_to_int<false>(kv,minus_delta,w_i);

        // Since u changed, update w_i b.c. at the surface
        Kokkos::single(Kokkos::PerThread(kv.team),[&](){
          constexpr int LAST_MID_PACK     = InfoM::LastPack;
          constexpr int LAST_MID_PACK_END = InfoM::LastPackEnd;
          constexpr int LAST_INT_PACK     = InfoI::LastPack;
          constexpr int LAST_INT_PACK_END = InfoI::LastPackEnd;
          constexpr auto g = PhysicalConstants::g;

          const auto gradphis = Homme::subview(m_geometry.m_gradphis,kv.ie);
          const auto v        = Homme::subview(m_state.m_v,kv.ie,np1);

          w_i(LAST_INT_PACK)[LAST_INT_PACK_END] = 
                (v(0,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END]*gradphis(0,igp,jgp) +
                 v(1,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END]*gradphis(1,igp,jgp)) / g;
        });
      });
    } else if (istate==1) {
      // Update phinh_i
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        // Subviews (note: recycle phi_ref for the temporary p_i)
        auto phinh_i     = Homme::subview(m_state.m_phinh_i,kv.ie,np1,igp,jgp);
        auto delta_phinh = Homme::subview(m_delta_phinh,kv.ie,igp,jgp);
        auto phi_ref     = Homme::subview(m_phi_ref,kv.team_idx,igp,jgp);
        auto dp_pt       = Homme::subview(dp,igp,jgp);
        auto p_i         = Homme::subview(m_phi_ref,kv.team_idx,igp,jgp);
        auto p           = Homme::subview(m_temp,kv.team_idx,igp,jgp);

        // Recompute hydrostatic phi with new theta
        m_elem_ops.compute_hydrostatic_p(kv,dp_pt,p_i,p);
        m_eos.compute_phi_i(kv,m_geometry.m_phis(kv.ie,igp,jgp),
                            Homme::subview(m_state.m_vtheta_dp,kv.ie,np1,igp,jgp),
                            p,
                            phi_ref);

        // // phinh_i(k) = phinh_i(k+1) - delta_phinh(k), so do a backward scan sum of -delta_phinh
        auto minus_delta = [&](const int ilev)->Scalar { return -delta_phinh(ilev); };
        ColumnOps::column_scan_mid_to_int<false>(kv,minus_delta,phinh_i);

        // Add hydrostatic phi
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV_P),
                             [&](const int ilev) {
          phinh_i(ilev) += phi_ref(ilev);
        });
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>
  get_state(const KernelVariables &kv, int np1, int var) const {
    assert(var>=0 && var<=4);
    switch (var) {
    case 0:
      return Homme::subview(m_state.m_v, kv.ie, np1, 0);
    case 1:
      return Homme::subview(m_state.m_v, kv.ie, np1, 1);
    case 2:
      return Homme::subview(m_state.m_vtheta_dp, kv.ie, np1);
    case 3:
      return Homme::subview(m_delta_w, kv.ie);
    case 4:
      return Homme::subview(m_delta_phinh, kv.ie);
    default:
      Kokkos::abort("RemapStateProvider: invalid variable index.\n");
      return ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]>();
    }
  }
};

} // namespace Remap
} // namespace Homme

#endif // HOMMEXX_REMAP_STATE_PROVIDER_HPP
