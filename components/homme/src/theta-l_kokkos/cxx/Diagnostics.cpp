#include "Diagnostics.hpp"

#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HybridVCoord.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"

#include "utilities/SyncUtils.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme
{

void Diagnostics::init (const ElementsState& state, const ElementsGeometry& geometry,
                             const HybridVCoord& hvcoord, const bool theta_hydrostatic_mode,
                             F90Ptr& elem_state_q_ptr,
                             F90Ptr& elem_accum_qvar_ptr,  F90Ptr& elem_accum_qmass_ptr,
                             F90Ptr& elem_accum_q1mass_ptr,F90Ptr& elem_accum_iener_ptr,
                             F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr)
{
  // Check state/geometry/hvcoord were inited
  assert (state.num_elems()>0);
  assert (geometry.num_elems()>0);
  assert (hvcoord.m_inited);

  // Check initialization was done right
  assert (m_num_elems==state.num_elems());

  m_state    = state;
  m_geometry = geometry;
  m_hvcoord  = hvcoord;
  m_theta_hydrostatic_mode = theta_hydrostatic_mode;

  m_eos.init(m_theta_hydrostatic_mode,m_hvcoord);
  m_elem_ops.init(m_hvcoord);

  // F90 ptr to array (n1,n2,...,nK,nelemd) can be stuffed directly in an unmanaged view
  // with scalar type Real*[nK]...[n2][n1] (with runtime dimension nelemd)
  h_Q      = decltype(h_Q)(elem_state_q_ptr, m_num_elems);
  h_Qvar   = decltype(h_Qvar)(elem_accum_qvar_ptr, m_num_elems);
  h_Qmass  = decltype(h_Qmass)(elem_accum_qmass_ptr, m_num_elems);
  h_Q1mass = decltype(h_Q1mass)(elem_accum_q1mass_ptr, m_num_elems);

  h_IEner  = decltype(h_IEner)(elem_accum_iener_ptr, m_num_elems);
  h_KEner  = decltype(h_KEner)(elem_accum_kener_ptr, m_num_elems);
  h_PEner  = decltype(h_PEner)(elem_accum_pener_ptr, m_num_elems);

  m_IEner  = decltype(m_IEner)("Internal  Energy", m_num_elems);
  m_KEner  = decltype(m_KEner)("Kinetic   Energy", m_num_elems);
  m_PEner  = decltype(m_PEner)("Potential Energy", m_num_elems);
}

int Diagnostics::requested_buffer_size () const {
  const int nteams = m_tu.get_num_concurrent_teams();

  constexpr int size_mid_scalar = NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar = NP*NP*NUM_LEV_P*VECTOR_SIZE;
  return nteams * (Buffers::num_3d_scalar_mid_buf*size_mid_scalar +
                   Buffers::num_3d_scalar_mid_buf*size_int_scalar);
}

void Diagnostics::init_buffers (const FunctorsBuffersManager& fbm) {
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

  const int nteams = m_tu.get_num_concurrent_teams();

  // If nteams is 0, something is wrong
  assert (nteams>0);

  m_buffers.pnh    = decltype(m_buffers.pnh)(mem,nteams);
  mem += nteams*NP*NP*NUM_LEV;

  m_buffers.exner  = decltype(m_buffers.exner)(mem,nteams);
  mem += nteams*NP*NP*NUM_LEV;

  m_buffers.phi    = decltype(m_buffers.phi)(mem,nteams);
  mem += nteams*NP*NP*NUM_LEV;

  m_buffers.dp_ref = decltype(m_buffers.dp_ref)(mem,nteams);
  mem += nteams*NP*NP*NUM_LEV;

  m_buffers.dpnh_dp_i = decltype(m_buffers.dpnh_dp_i)(mem,nteams);
}

void Diagnostics::sync_diags_to_host () {
  Kokkos::deep_copy(h_IEner,m_IEner);
  Kokkos::deep_copy(h_PEner,m_PEner);
  Kokkos::deep_copy(h_KEner,m_KEner);
}

void Diagnostics::run_diagnostics (const bool before_advance, const int ivar)
{
  prim_diag_scalars(before_advance, ivar);
  prim_energy_halftimes(before_advance, ivar);
}

void Diagnostics::prim_diag_scalars (const bool before_advance, const int ivar)
{
  // Get simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  assert(params.params_set);

  // Get time info
  TimeLevel& tl = Context::singleton().get<TimeLevel>();

  // Make sure tracers timelevels are updated
  tl.update_tracers_levels(params.qsplit);

  // Pick tracers time-level, depending on when this routine was called
  int t2_qdp;
  if (before_advance) {
    t2_qdp = tl.n0_qdp;
  } else {
    t2_qdp = tl.np1_qdp;
  }

  const Tracers& tracers = Context::singleton().get<Tracers>();

  sync_to_host(tracers.Q,h_Q);

  // Copy back tracers concentration and mass
  auto qdp_h = Kokkos::create_mirror_view(tracers.qdp);
  Kokkos::deep_copy(qdp_h,tracers.qdp);
  for (int ie=0; ie<m_num_elems; ++ie) {
    for (int iq=0; iq<params.qsize; ++iq) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          Real accum_qdp_q = 0;
          Real accum_qdp = 0;

          HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> qdp(&qdp_h(ie, t2_qdp, iq, igp, jgp, 0)[0]);
          for (int level=0; level<NUM_PHYSICAL_LEV; ++level) {
            accum_qdp_q += qdp(level)*h_Q(ie, iq, level, igp, jgp);
            accum_qdp   += qdp(level);
          }
          h_Qvar(ie, ivar, iq, igp, jgp) = accum_qdp_q;

          h_Qmass(ie, ivar, iq, igp, jgp) = accum_qdp;
          h_Q1mass(ie, iq, igp, jgp) = accum_qdp;
        }
      }
    }
  }
}

void Diagnostics::prim_energy_halftimes (const bool before_advance, const int ivar)
{
  m_ivar = ivar;

  // Get simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  assert(params.params_set);

  // Get time info
  TimeLevel& tl = Context::singleton().get<TimeLevel>();

  // Make sure tracers timelevels are updated
  tl.update_tracers_levels(params.qsplit);

  // Pick tracers time-level, depending on when this routine was called
  if (before_advance) {
    t1 = tl.n0;
    t1_qdp = tl.n0_qdp;
  } else {
    t1 = tl.np1;
    t1_qdp = tl.np1_qdp;
  }

  Kokkos::parallel_for(m_policy, *this);

  Kokkos::deep_copy(h_KEner, m_KEner);
}

} // namespace Homme
