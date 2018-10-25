/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Context.hpp"
#include "CamForcing.hpp"
#include "Tracers.hpp"
#include "Elements.hpp"
#include "TimeLevel.hpp"
#include "HybridVCoord.hpp"
#include "SimulationParams.hpp"
#include "KernelVariables.hpp"
#include "vector/vector_pragmas.hpp"
#include "profiling.hpp"

namespace Homme {

void state_forcing(
    const ExecViewUnmanaged<const Scalar * [NP][NP][NUM_LEV]> &f_t,
    const ExecViewUnmanaged<const Scalar * [2][NP][NP][NUM_LEV]> &f_m,
    const int &np1, const Real &dt,
    const ExecViewUnmanaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> &t,
    const ExecViewUnmanaged<Scalar * [NUM_TIME_LEVELS][2][NP][NP][NUM_LEV]> &
        v) {
  const int num_e = f_t.extent_int(0);
  Kokkos::parallel_for(
      "state temperature forcing",
      Kokkos::RangePolicy<ExecSpace>(0, num_e * NP * NP * NUM_LEV),
      KOKKOS_LAMBDA(const int & idx) {
        const int ie = ((idx / NUM_LEV) / NP) / NP;
        const int igp = ((idx / NUM_LEV) / NP) % NP;
        const int jgp = (idx / NUM_LEV) % NP;
        const int k = idx % NUM_LEV;
        t(ie, np1, igp, jgp, k) += dt * f_t(ie, igp, jgp, k);
      });
  Kokkos::parallel_for(
      "state velocity forcing",
      Kokkos::RangePolicy<ExecSpace>(0, num_e * 2 * NP * NP * NUM_LEV),
      KOKKOS_LAMBDA(const int & idx) {
        const int ie = (((idx / NUM_LEV) / NP) / NP) / 2;
        const int dim = (((idx / NUM_LEV) / NP) / NP) % 2;
        const int igp = ((idx / NUM_LEV) / NP) % NP;
        const int jgp = (idx / NUM_LEV) % NP;
        const int k = idx % NUM_LEV;
        v(ie, np1, dim, igp, jgp, k) += dt * f_m(ie, dim, igp, jgp, k);
      });
}

void tracer_forcing(
    const ExecViewUnmanaged<const Scalar * [QSIZE_D][NP][NP][NUM_LEV]> &f_q,
    const HybridVCoord &hvcoord, const TimeLevel &tl, const int &num_q,
    const MoistDry &moisture, const double &dt,
    const ExecViewManaged<Real * [NUM_TIME_LEVELS][NP][NP]> &ps_v,
    const ExecViewManaged<
        Scalar * [Q_NUM_TIME_LEVELS][QSIZE_D][NP][NP][NUM_LEV]> &qdp,
    const ExecViewManaged<Scalar * [QSIZE_D][NP][NP][NUM_LEV]> &Q) {

  const int num_e = ps_v.extent_int(0);
  const int np1 = tl.n0;
  const int np1_qdp = tl.n0_qdp;

  if (moisture == MoistDry::MOIST) {
    // Remove the m_fq_ps_v buffer since it's not actually needed.
    // Instead apply the forcing to m_ps_v directly
    // Bonus - one less parallel reduce in dry cases!

    // This block must go force so that qdp_s, which is read but not written,
    // has its original value.

    // This conserves the dry mass in the process of forcing tracer 0
    const auto policy = Homme::get_default_team_policy<ExecSpace>(num_e);
    Kokkos::parallel_for("tracer FQ_PS forcing", policy,
                         KOKKOS_LAMBDA(const TeamMember & team) {
      KernelVariables kv(team);
      const int &ie = kv.ie;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                           [&](const int pt_idx) {
        const int igp = pt_idx / NP;
        const int jgp = pt_idx % NP;

        Real ps_v_forcing = 0.0;

        Dispatch<ExecSpace>::parallel_reduce(
            kv.team, Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
            [&](const int &k, Real &accumulator) {
              const int ilev = k / VECTOR_SIZE;
              const int vlev = k % VECTOR_SIZE;
              Real v1 = dt * f_q(ie, 0, igp, jgp, ilev)[vlev];
              const Real &qdp_s =
                  qdp(ie, np1_qdp, 0, igp, jgp, ilev)[vlev];
              if (qdp_s + v1 < 0.0 && v1 < 0.0) {
                if (qdp_s < 0.0) {
                  v1 = 0.0;
                } else {
                  v1 = -qdp_s;
                }
              }
              accumulator += v1;
            },
            ps_v_forcing);
        ps_v(ie, np1, igp, jgp) += ps_v_forcing;
      });
    });
  }

  Kokkos::parallel_for(
      "tracer qdp forcing",
      Kokkos::RangePolicy<ExecSpace>(0, num_e * num_q * NP * NP * NUM_LEV),
      KOKKOS_LAMBDA(const int & idx) {
        const int ie = (((idx / NUM_LEV) / NP) / NP) / num_q;
        const int iq = (((idx / NUM_LEV) / NP) / NP) % num_q;
        const int igp = ((idx / NUM_LEV) / NP) % NP;
        const int jgp = (idx / NUM_LEV) % NP;
        const int k = idx % NUM_LEV;
        Scalar v1 = dt * f_q(ie, iq, igp, jgp, k);
        Scalar &qdp_s = qdp(ie, np1_qdp, iq, igp, jgp, k);
        VECTOR_SIMD_LOOP
        for (int vlev = 0; vlev < VECTOR_SIZE; ++vlev) {
          if (qdp_s[vlev] + v1[vlev] < 0.0 && v1[vlev] < 0.0) {
            if (qdp_s[vlev] < 0.0) {
              v1[vlev] = 0.0;
            } else {
              v1[vlev] = -qdp_s[vlev];
            }
          }
          qdp_s[vlev] += v1[vlev];
        }
      });

  Kokkos::parallel_for("tracer forcing ps_v",
                       Kokkos::RangePolicy<ExecSpace>(0, num_e * num_q * NP *
                                                             NP * NUM_LEV),
                       KOKKOS_LAMBDA(const int & idx) {
    const int ie = (((idx / NUM_LEV) / NP) / NP) / num_q;
    const int iq = (((idx / NUM_LEV) / NP) / NP) % num_q;
    const int igp = ((idx / NUM_LEV) / NP) % NP;
    const int jgp = (idx / NUM_LEV) % NP;
    const int k = idx % NUM_LEV;

    const Scalar dp = hvcoord.hybrid_ai_delta(k) * hvcoord.ps0 +
                      hvcoord.hybrid_bi_delta(k) * ps_v(ie, np1, igp, jgp);
    Q(ie, iq, igp, jgp, k) = qdp(ie, np1_qdp, iq, igp, jgp, k) / dp;
  });
}

void apply_cam_forcing(const Real &dt) {
  GPTLstart("ApplyCAMForcing");
  const Elements &elems = Context::singleton().get_elements();
  const TimeLevel &tl = Context::singleton().get_time_level();

  state_forcing(elems.m_ft, elems.m_fm, tl.n0, dt, elems.m_t, elems.m_v);

  const SimulationParams &sim_params =
      Context::singleton().get_simulation_params();
  const HybridVCoord &hvcoord = Context::singleton().get_hvcoord();
  Tracers &tracers = Context::singleton().get_tracers();
  if(tracers.fq.data() == nullptr) {
    tracers.fq = decltype(tracers.fq)("fq", elems.num_elems());
  }
  tracer_forcing(tracers.fq, hvcoord, tl, tracers.num_tracers(),
                 sim_params.moisture, dt, elems.m_ps_v, tracers.qdp, tracers.Q);
  GPTLstop("ApplyCAMForcing");
}

void apply_cam_forcing_dynamics(const Real &dt) {
  GPTLstart("ApplyCAMForcing_dynamics");
  const Elements &elems = Context::singleton().get_elements();
  const TimeLevel &tl = Context::singleton().get_time_level();
  state_forcing(elems.m_ft, elems.m_fm, tl.n0, dt, elems.m_t, elems.m_v);
  GPTLstop("ApplyCAMForcing_dynamics");
}

// ---------------------------------------- //
}
