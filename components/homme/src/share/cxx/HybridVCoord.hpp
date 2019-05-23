/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYBRID_V_COORD_HPP
#define HOMMEXX_HYBRID_V_COORD_HPP

#include "KernelVariables.hpp"
#include "Types.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme
{

struct HybridVCoord
{
  HybridVCoord () : m_inited(false) {}

  // This method should only be called from the host
  void init(const Real ps0_in,
            CRCPtr hybrid_am_ptr,
            CRCPtr hybrid_ai_ptr,
            CRCPtr hybrid_bm_ptr,
            CRCPtr hybrid_bi_ptr);

  void random_init(int seed);
  void compute_deltas ();
  void compute_eta ();

  Real ps0;
  Real hybrid_ai0;

  // hybrid ai
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> hybrid_ai;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_am;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_ai_delta;

  // hybrid bi
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> hybrid_bi;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_bm;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_bi_delta;

  // Eta levels (at midpoints and interfaces)
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> etai;
  ExecViewManaged<Scalar[NUM_LEV]>         etam;

  // So far these seem to never be used
  // ExecViewManaged<Real[NUM_INTERFACE_LEV]> delta_etai;
  // ExecViewManaged<Scalar[NUM_LEV]>         delta_etam;

  ExecViewManaged<Scalar[NUM_LEV]> dp0;

  bool m_inited;

  // This reference dp is computed several times in the code, so we decide
  KOKKOS_INLINE_FUNCTION
  void compute_dp_ref (KernelVariables& kv,
                       ExecViewUnmanaged<const Real[NP][NP]> ps,
                       ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV]> dp) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int& point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;
      auto point_dp = Homme::subview(dp,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int& ilev) {
        point_dp(ilev) = hybrid_ai_delta(ilev)*ps0 +
                         hybrid_bi_delta(ilev)*ps(igp,jgp);
      });
    });

    // Should we remove this and let the user put it from outside if needed?
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void compute_ps_ref (KernelVariables& kv,
                       ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]> dp,
                       ExecViewUnmanaged<      Real   [NP][NP]> ps) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      Scalar tmp;
      Dispatch<>::parallel_reduce(kv.team,Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                  [&](const int ilev, Scalar& sum) {
        sum += dp(igp,jgp,ilev);
      },tmp);
      if (OnGpu<ExecSpace>::value) {
        ps(igp,jgp) = hybrid_ai(0)*ps0 + tmp[0];
      } else {
        ps(igp,jgp) = hybrid_ai(0)*ps0 + tmp.reduce_add();
      }
    });
    kv.team_barrier();

    // Should we remove this and let the user put it from outside if needed?
    kv.team_barrier();
  }
};

} // namespace Homme

#endif // HOMMEXX_HYBRID_V_COORD_HPP
