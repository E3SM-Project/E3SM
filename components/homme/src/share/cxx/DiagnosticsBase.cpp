#include "DiagnosticsBase.hpp"

#include "Context.hpp"
#include "Elements.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"

#include "utilities/SyncUtils.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme
{

void DiagnosticsBase::init (const int num_elems,
                            F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,
                            F90Ptr& elem_accum_qmass_ptr, F90Ptr& elem_accum_q1mass_ptr)
{
  assert (num_elems>0);
  m_num_elems = num_elems;

  // F90 ptr to array (n1,n2,...,nK,nelemd) can be stuffed directly in an unmanaged view
  // with scalar type Real*[nK]...[n2][n1] (with runtime dimension nelemd)
  h_Q      = HostViewUnmanaged<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>(elem_state_q_ptr, num_elems);
  h_Qvar   = HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>(elem_accum_qvar_ptr, num_elems);
  h_Qmass  = HostViewUnmanaged<Real*[4][QSIZE_D][NP][NP]>(elem_accum_qmass_ptr, num_elems);
  h_Q1mass = HostViewUnmanaged<Real*   [QSIZE_D][NP][NP]>(elem_accum_q1mass_ptr, num_elems);
}

void DiagnosticsBase::prim_diag_scalars (const bool before_advance, const int ivar)
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

  if (params.time_step_type>0) {
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
}

} // namespace Homme
