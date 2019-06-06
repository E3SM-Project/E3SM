#include "DiagnosticsTheta.hpp"

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

void DiagnosticsTheta::init (const ElementsState& state, const ElementsGeometry& geometry,
                             const HybridVCoord& hvcoord, const bool theta_hydrostatic_mode,
                             F90Ptr& elem_state_q_ptr,
                             F90Ptr& elem_accum_qvar_ptr,  F90Ptr& elem_accum_qmass_ptr,
                             F90Ptr& elem_accum_q1mass_ptr,F90Ptr& elem_accum_iener_ptr,
                             F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr)
{
  DiagnosticsBase::init(state.num_elems(), elem_state_q_ptr, elem_accum_qvar_ptr, elem_accum_qmass_ptr, elem_accum_q1mass_ptr);

  m_state    = state;
  m_hvcoord  = hvcoord;
  m_theta_hydrostatic_mode = theta_hydrostatic_mode;

  h_IEner  = HostViewUnmanaged<Real*[4][NP][NP]>(elem_accum_iener_ptr, m_num_elems);
  h_KEner  = HostViewUnmanaged<Real*[4][NP][NP]>(elem_accum_kener_ptr, m_num_elems);
  h_PEner  = HostViewUnmanaged<Real*[4][NP][NP]>(elem_accum_pener_ptr, m_num_elems);

  m_IEner  = ExecViewManaged<Real*[4][NP][NP]>("Internal  Energy", m_num_elems);
  m_KEner  = ExecViewManaged<Real*[4][NP][NP]>("Kinetic   Energy", m_num_elems);
  m_PEner  = ExecViewManaged<Real*[4][NP][NP]>("Potential Energy", m_num_elems);
}

int DiagnosticsTheta::requested_buffer_size () const {
  const int nteams = get_num_concurrent_teams(Homme::get_default_team_policy<ExecSpace,EnergyHalfTimesTag>(m_num_elems));

  constexpr int size_mid_scalar = NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar = NP*NP*NUM_LEV_P*VECTOR_SIZE;
  return nteams * (Buffers::num_3d_scalar_mid_buf*size_mid_scalar +
                   Buffers::num_3d_scalar_mid_buf*size_int_scalar);
}

void DiagnosticsTheta::init_buffers (const FunctorsBuffersManager& fbm) {
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());

  const int nteams = get_num_concurrent_teams(Homme::get_default_team_policy<ExecSpace,EnergyHalfTimesTag>(m_num_elems));

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

void DiagnosticsTheta::prim_energy_halftimes (const bool before_advance, const int ivar)
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

  auto policy = Homme::get_default_team_policy<ExecSpace,EnergyHalfTimesTag>(m_num_elems);
  Kokkos::parallel_for(policy, *this);

  Kokkos::deep_copy(h_KEner, m_KEner);
}

} // namespace Homme
