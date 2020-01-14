#include "Diagnostics.hpp"

#include "Context.hpp"
#include "Elements.hpp"
#include "HybridVCoord.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "TimeLevel.hpp"
#include "Tracers.hpp"

#include "utilities/SyncUtils.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme
{

void Diagnostics::init (const int num_elems,
                        F90Ptr& elem_state_q_ptr, F90Ptr& elem_accum_qvar_ptr,
                        F90Ptr& elem_accum_qmass_ptr, F90Ptr& elem_accum_q1mass_ptr,
                        F90Ptr& elem_accum_iener_ptr, F90Ptr& elem_accum_iener_wet_ptr,
                        F90Ptr& elem_accum_kener_ptr, F90Ptr& elem_accum_pener_ptr)
{
  assert (num_elems>0);
  m_num_elems = num_elems;

  // F90 ptr to array (n1,n2,...,nK,nelemd) can be stuffed directly in an unmanaged view
  // with scalar type Real*[nK]...[n2][n1] (with runtime dimension nelemd)
  h_Q      = decltype(h_Q)(elem_state_q_ptr, m_num_elems);
  h_Qvar   = decltype(h_Qvar)(elem_accum_qvar_ptr, m_num_elems);
  h_Qmass  = decltype(h_Qmass)(elem_accum_qmass_ptr, m_num_elems);
  h_Q1mass = decltype(h_Q1mass)(elem_accum_q1mass_ptr, m_num_elems);

  h_IEner     = decltype(h_IEner)(elem_accum_iener_ptr, m_num_elems);
  h_IEner_wet = decltype(h_IEner_wet)(elem_accum_iener_wet_ptr, m_num_elems);
  h_KEner     = decltype(h_KEner)(elem_accum_kener_ptr, m_num_elems);
  h_PEner     = decltype(h_PEner)(elem_accum_pener_ptr, m_num_elems);
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
  // Get simulation params
  SimulationParams& params = Context::singleton().get<SimulationParams>();
  assert(params.params_set);

  // Get time info
  TimeLevel& tl = Context::singleton().get<TimeLevel>();

  // Make sure tracers timelevels are updated
  tl.update_tracers_levels(params.qsplit);

  // Pick tracers time-level, depending on when this routine was called
  int t1,t1_qdp;
  if (before_advance) {
    t1 = tl.n0;
    t1_qdp = tl.n0_qdp;
  } else {
    t1 = tl.np1;
    t1_qdp = tl.np1_qdp;
  }

  // Getting stuff we need on host
  const HybridVCoord& hvcoord = Context::singleton().get<HybridVCoord>();
  auto dhyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai_delta);
  auto dhybi = Kokkos::create_mirror_view(hvcoord.hybrid_bi_delta);
  Kokkos::deep_copy(dhyai,hvcoord.hybrid_ai_delta);
  Kokkos::deep_copy(dhybi,hvcoord.hybrid_bi_delta);

  const Elements& elements = Context::singleton().get<Elements>();
  const Tracers& tracers = Context::singleton().get<Tracers>();

  auto h_ps_v = Kokkos::create_mirror_view(elements.m_state.m_ps_v);
  auto h_phis = Kokkos::create_mirror_view(elements.m_geometry.m_phis);
  auto h_T = Kokkos::create_mirror_view(elements.m_state.m_t);
  auto h_v = Kokkos::create_mirror_view(elements.m_state.m_v);
  auto qdp_h = Kokkos::create_mirror_view(tracers.qdp);
  Kokkos::deep_copy(h_ps_v, elements.m_state.m_ps_v);
  Kokkos::deep_copy(h_phis, elements.m_geometry.m_phis);
  Kokkos::deep_copy(h_T, elements.m_state.m_t);
  Kokkos::deep_copy(h_v, elements.m_state.m_v);

  if (params.use_cpstar) {
    Kokkos::deep_copy(qdp_h, tracers.qdp);
  }

  auto virtual_specific_heat = [] (const Real rin)->Real{
    using namespace PhysicalConstants;
    return cp*(1.0 + (Cpwater_vapor/cp - 1.0)*rin);
  };

  // !   IE   Cp*dpdn*T  + (Cpv-Cp) Qdpdn*T
  // !        Cp*dpdn(n)*T(n+1) + (Cpv-Cp) Qdpdn(n)*T(n+1)
  // !        [Cp + (Cpv-Cp) Q(n)] *dpdn(n)*T(n+1)

  Real cp_star1;
  constexpr Real cp = PhysicalConstants::cp;
  for (int ie=0; ie<m_num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++jgp) {
        // Accumulation values
        Real IEner = 0.0;
        Real IEner_wet = 0.0;
        Real KEner = 0.0;
        Real PEner = 0.0;

        auto u = Homme::subview(h_v, ie, t1, 0,igp, jgp);
        auto v = Homme::subview(h_v, ie, t1, 1,igp, jgp);
        auto T = Homme::subview(h_T, ie, t1, igp, jgp);

        for (int ilev=0; ilev<NUM_LEV; ++ilev) {
          Scalar dpt1 = dhyai(ilev)*hvcoord.ps0 + dhybi(ilev)*h_ps_v(ie, t1, igp, jgp);

          const int vector_end = (ilev==NUM_LEV-1 ? (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE : VECTOR_SIZE - 1);
          for (int ivec=0; ivec<=vector_end; ++ivec) {
            if (params.use_cpstar) {
              Real qval_t1 = qdp_h(ie, t1_qdp, 0, igp, jgp, ilev)[ivec]/dpt1[ivec];
              cp_star1 = virtual_specific_heat(qval_t1);
            } else {
              cp_star1 = cp;
            }

            IEner     += cp_star1*T(ilev)[ivec]*dpt1[ivec];
            IEner_wet += (cp_star1 - cp)*T(ilev)[ivec]*dpt1[ivec];
            KEner     += (std::pow(u(ilev)[ivec],2) + std::pow(v(ilev)[ivec],2))*0.5*dpt1[ivec];
            PEner     += h_phis(ie, igp, jgp)*dpt1[ivec];
          }
        } // ilev
        h_IEner(ie, ivar, igp, jgp) = IEner;
        h_IEner_wet(ie, ivar, igp, jgp)  = IEner_wet;
        h_KEner(ie, ivar, igp, jgp) = KEner;
        h_PEner(ie, ivar, igp, jgp) = PEner;
      } // jgp
    } // igp
  } // ie
}

} // namespace Homme
