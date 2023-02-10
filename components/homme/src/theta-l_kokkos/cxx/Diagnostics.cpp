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

#include "profiling.hpp"

namespace Homme
{

void Diagnostics::init (const ElementsState& state, const ElementsGeometry& geometry,
                        const HybridVCoord& hvcoord,  const Tracers& tracers,
                        F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,
                        F90Ptr& elem_accum_qmass_ptr, F90Ptr& elem_accum_q1mass_ptr,
                        F90Ptr& elem_accum_iener_ptr, F90Ptr& elem_accum_kener_ptr,
                        F90Ptr& elem_accum_pener_ptr) {
  // Check state/geometry/hvcoord were inited
  assert (state.num_elems()>0);
  assert (geometry.num_elems()>0);
  assert (hvcoord.m_inited);
  assert (tracers.inited());

  // Check initialization was done right
  assert (m_num_elems==state.num_elems());
  assert (m_num_tracers==tracers.num_tracers());

  m_state    = state;
  m_geometry = geometry;
  m_hvcoord  = hvcoord;
  m_tracers  = tracers;

  m_eos.init(m_theta_hydrostatic_mode,m_hvcoord);
  m_elem_ops.init(m_hvcoord);

  // F90 ptr to array (n1,n2,...,nK,nelemd) can be stuffed directly in an unmanaged view
  // with scalar type Real*[nK]...[n2][n1] (with runtime dimension nelemd)
  h_Qvar   = decltype(h_Qvar)  (elem_accum_qvar_ptr,   m_num_elems);
  h_Qmass  = decltype(h_Qmass) (elem_accum_qmass_ptr,  m_num_elems);
  h_Q1mass = decltype(h_Q1mass)(elem_accum_q1mass_ptr, m_num_elems);

  d_Qvar   = decltype(d_Qvar)  ("Qvar",   m_num_elems);
  d_Qmass  = decltype(d_Qmass) ("Qmass",  m_num_elems);
  d_Q1mass = decltype(d_Q1mass)("Q1mass", m_num_elems);

  h_IEner  = decltype(h_IEner) (elem_accum_iener_ptr, m_num_elems);
  h_KEner  = decltype(h_KEner) (elem_accum_kener_ptr, m_num_elems);
  h_PEner  = decltype(h_PEner) (elem_accum_pener_ptr, m_num_elems);

  d_IEner  = decltype(d_IEner) ("Internal  Energy", m_num_elems);
  d_KEner  = decltype(d_KEner) ("Kinetic   Energy", m_num_elems);
  d_PEner  = decltype(d_PEner) ("Potential Energy", m_num_elems);
}

int Diagnostics::requested_buffer_size () const {
  const int nslots = m_tu.get_num_ws_slots();

  constexpr int size_mid_scalar = NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar = NP*NP*NUM_LEV_P*VECTOR_SIZE;
  const int num_3d_int = m_theta_hydrostatic_mode ? Buffers::num_3d_scalar_int_buf_hy
                                                  : Buffers::num_3d_scalar_int_buf_nh;
  return nslots * (Buffers::num_3d_scalar_mid_buf*size_mid_scalar +
                   num_3d_int*size_int_scalar);
}

void Diagnostics::init_buffers (const FunctorsBuffersManager& fbm) {
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

  const int nslots = m_tu.get_num_ws_slots();

  // If nslots is 0, something is wrong
  assert (nslots>0);

  m_buffers.pnh    = decltype(m_buffers.pnh)(mem,nslots);
  mem += nslots*NP*NP*NUM_LEV;

  m_buffers.exner  = decltype(m_buffers.exner)(mem,nslots);
  mem += nslots*NP*NP*NUM_LEV;

  m_buffers.phi    = decltype(m_buffers.phi)(mem,nslots);
  mem += nslots*NP*NP*NUM_LEV;

  m_buffers.dp_ref = decltype(m_buffers.dp_ref)(mem,nslots);
  mem += nslots*NP*NP*NUM_LEV;

  m_buffers.dpnh_dp_i = decltype(m_buffers.dpnh_dp_i)(mem,nslots);
  mem += nslots*NP*NP*NUM_LEV_P;

  if (m_theta_hydrostatic_mode) {
    m_buffers.phi_i = decltype(m_buffers.phi_i)(mem,nslots);
  }
}

void Diagnostics::sync_diagnostics_to_host () {
  Kokkos::deep_copy(h_IEner,  d_IEner);
  Kokkos::deep_copy(h_PEner,  d_PEner);
  Kokkos::deep_copy(h_KEner,  d_KEner);
  Kokkos::deep_copy(h_Qvar,   d_Qvar);
  Kokkos::deep_copy(h_Qmass,  d_Qmass);
  Kokkos::deep_copy(h_Q1mass, d_Q1mass);
}

void Diagnostics::run_diagnostics (const bool before_advance, const int ivar)
{
  GPTLstart("prim_diag");
  prim_diag_scalars(before_advance, ivar);
  prim_energy_halftimes(before_advance, ivar);
  GPTLstop("prim_diag");
}

void Diagnostics::prim_diag_scalars (const bool before_advance, const int ivar)
{
  m_ivar = ivar;

  // Get simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  assert(params.params_set);

  // Get time info
  TimeLevel& tl = Context::singleton().get<TimeLevel>();

  // Make sure tracers timelevels are updated
  tl.update_tracers_levels(params.dt_tracer_factor);

  // Pick tracers time-level, depending on when this routine was called
  if (before_advance) {
    t2_qdp = tl.n0_qdp;
  } else {
    t2_qdp = tl.np1_qdp;
  }

  Kokkos::parallel_for(m_policy_diag_scalars, *this);
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
  tl.update_tracers_levels(params.dt_tracer_factor);

  // Pick tracers time-level, depending on when this routine was called
  if (before_advance) {
    t1 = tl.n0;
    t1_qdp = tl.n0_qdp;
  } else {
    t1 = tl.np1;
    t1_qdp = tl.np1_qdp;
  }

  Kokkos::parallel_for(m_policy_energy_halftimes, *this);
}

} // namespace Homme
