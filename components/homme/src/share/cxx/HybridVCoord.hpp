/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_HYBRID_V_COORD_HPP
#define HOMMEXX_HYBRID_V_COORD_HPP

#include "KernelVariables.hpp"
#include "PhysicalConstants.hpp"
#include "ColumnOps.hpp"
#include "Types.hpp"
#include "utilities/SubviewUtils.hpp"

namespace Homme
{

class HybridVCoord
{
public:
  KOKKOS_INLINE_FUNCTION
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
  ExecViewUnmanaged<Real[NUM_INTERFACE_LEV]> hybrid_ai;
  ExecViewManaged<Scalar[NUM_LEV_P]> hybrid_ai_packed;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_am;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_ai_delta;

  // hybrid bi
  ExecViewUnmanaged<Real[NUM_INTERFACE_LEV]> hybrid_bi;
  ExecViewManaged<Scalar[NUM_LEV_P]> hybrid_bi_packed;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_bm;
  ExecViewManaged<Scalar[NUM_LEV]> hybrid_bi_delta;

  // Eta levels (at midpoints and interfaces)
  ExecViewManaged<Real[NUM_INTERFACE_LEV]> etai;
  ExecViewManaged<Scalar[NUM_LEV]>         etam;

  // So far these seem to never be used
  // ExecViewManaged<Real[NUM_INTERFACE_LEV]> delta_etai;
  // ExecViewManaged<Scalar[NUM_LEV]>         delta_etam;

  ExecViewManaged<Scalar[NUM_LEV]> dp0;
  ExecViewManaged<Scalar[NUM_LEV]> exner0;

  bool m_inited;

  // This reference p is computed several times in the code, so we decide
  KOKKOS_INLINE_FUNCTION
  void compute_dp_ref (const KernelVariables& kv,
                       const ExecViewUnmanaged<const Real[NP][NP]>& ps,
                       const ExecViewUnmanaged<Scalar [NP][NP][NUM_LEV]>& dp) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int& point_idx) {
      const int igp = point_idx / NP;
      const int jgp = point_idx % NP;
      auto point_dp = Homme::subview(dp,igp,jgp);
      const Real point_ps = ps(igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int& ilev) {
        point_dp(ilev) = hybrid_ai_delta(ilev)*ps0 +
                         hybrid_bi_delta(ilev)*point_ps;
      });
    });

    // Should we remove this and let the user put it from outside if needed?
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void compute_dp_ref (const KernelVariables& kv, const Real ps,
                       const ExecViewUnmanaged<Scalar [NUM_LEV]>& dp) const
  {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int& ilev) {
      dp(ilev) = hybrid_ai_delta(ilev)*ps0 +
                 hybrid_bi_delta(ilev)*ps;
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_ps_ref_from_dp (const KernelVariables& kv,
                       const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]>& dp,
                       const ExecViewUnmanaged<      Real   [NP][NP]>& ps) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto dp_ij = Homme::subview(dp,igp,jgp);

      ColumnOps::column_reduction<NUM_PHYSICAL_LEV>(kv,dp_ij,ps(igp,jgp));
      kv.team_barrier();

      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        ps(igp,jgp) += hybrid_ai0*ps0;
      });
    });
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void compute_ps_ref_from_phis (const KernelVariables& kv,
                       const ExecViewUnmanaged<const Real [NP][NP]>& phis,
                       const ExecViewUnmanaged<      Real [NP][NP]>& ps) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      ps(igp,jgp) = ps0*exp( - phis(igp,jgp)/ (PhysicalConstants::Rgas*300) );
    });
    kv.team_barrier();
  }
};

} // namespace Homme

#endif // HOMMEXX_HYBRID_V_COORD_HPP
