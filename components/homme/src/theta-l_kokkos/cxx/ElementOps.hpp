#ifndef HOMMEXX_ELEMENT_OPS_HPP
#define HOMMEXX_ELEMENT_OPS_HPP

#include "Types.hpp"
#include "KernelVariables.hpp"
#include "HybridVCoord.hpp"
#include "ColumnOps.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/BfbUtils.hpp"

namespace Homme {

class ElementOps {
public:
  KOKKOS_INLINE_FUNCTION
  ElementOps () = default;

  KOKKOS_INLINE_FUNCTION
  ~ElementOps () = default;

  void init (const HybridVCoord& hvcoord) {
    m_hvcoord = hvcoord;
  }

  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  void get_R_star (const KernelVariables& kv,
                   const bool use_moisture,
                   const InputProvider& Q,
                   const ExecViewUnmanaged<Scalar[NUM_LEV]>& R) const {
    using namespace PhysicalConstants;
    if (use_moisture) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        R(ilev) = (Rgas + (Rwater_vapor-Rgas)*Q(ilev));
      });
    } else {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        R(ilev) = Rgas;
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_hydrostatic_p (const KernelVariables& kv,
                              const ExecViewUnmanaged<const Scalar[NUM_LEV  ]>& dp,
                              const ExecViewUnmanaged<      Scalar[NUM_LEV_P]>& p_i,
                              const ExecViewUnmanaged<      Scalar[NUM_LEV  ]>& pi) const
  {
    p_i(0)[0] = m_hvcoord.hybrid_ai0*m_hvcoord.ps0;
    ColumnOps::column_scan_mid_to_int<true>(kv,dp,p_i);
    ColumnOps::compute_midpoint_values(kv,p_i,pi);
  }

  template<typename InputProvider>
  KOKKOS_INLINE_FUNCTION
  void compute_theta_ref (const KernelVariables& kv,
                          const InputProvider& p,
                          const ExecViewUnmanaged<Scalar[NUM_LEV]>& theta_ref) const {
    assert (m_hvcoord.m_inited);
    // theta_ref = T0/exner + T1, with T0,T1 fixed
    // exner = (p/p0)^k
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      using namespace PhysicalConstants;
      // Compute exner, store in theta_ref
      // TODO: F90 does p(k) = (p_i(k)+p_i(k+1)) / (2*p0).
      //       If this is a non BFB source, incorporate p0 scaling
      //       in the calculation of p
#ifdef HOMMEXX_BFB_TESTING
      theta_ref(ilev) = bfb_pow(p(ilev)/p0,kappa);
#else
      theta_ref(ilev) = pow(p(ilev)/p0,kappa);
#endif

      // Compute theta_ref
      theta_ref(ilev) = T0/theta_ref(ilev) + T1;
    });
  }

private:

  static constexpr Real TREF = 288.0;
  static constexpr Real T1 = 0.0065f*TREF*PhysicalConstants::cp/PhysicalConstants::g;
  static constexpr Real T0 = TREF-T1;

  HybridVCoord    m_hvcoord;
};

} // namespace Homme

#endif // HOMMEXX_ELEMENT_OPS_HPP
