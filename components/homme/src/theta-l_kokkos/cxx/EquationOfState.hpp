#ifndef HOMMEXX_EQUATION_OF_STATE_HPP
#define HOMMEXX_EQUATION_OF_STATE_HPP

#include "Types.hpp"
#include "ElementOps.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/VectorUtils.hpp"
#include "utilities/ViewUtils.hpp"

namespace Homme {

class EquationOfState {
public:

  EquationOfState () = default;

  void init (const bool theta_hydrostatic_mode,
             const HybridVCoord& hvcoord) {
    m_theta_hydrostatic_mode = theta_hydrostatic_mode;
    m_hvcoord = hvcoord;
  }

  KOKKOS_INLINE_FUNCTION
  void compute_pnh_and_exner (KernelVariables& kv,
                              ExecViewUnmanaged<const Scalar[NUM_LEV  ]> vtheta_dp,
                              ExecViewUnmanaged<const Scalar[NUM_LEV_P]> phi_i,
                              ExecViewUnmanaged<const Scalar[NUM_LEV  ]> pi,
                              ExecViewUnmanaged<      Scalar[NUM_LEV  ]> pnh,
                              ExecViewUnmanaged<      Scalar[NUM_LEV  ]> exner) const
  {
    if (m_theta_hydrostatic_mode) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        // Copy pi into pnh, and set exner = (pi/p0)^kappa
        pnh(ilev) = pi(ilev);

        // Avoid temporaries
        exner(ilev) = pi(ilev);
        exner(ilev) /= PhysicalConstants::p0;
        pow_update(exner(ilev),PhysicalConstants::kappa);
      });
    } else {
      // Compute p_over_exner = Rgas*vtheta_dp/delta(phi_i)
      // Then, pnh = p0 (p_over_exner/p0)^(1/(1-kappa))
      //       exner = pnh/p_over_exner

      // To avoid temporaries, store p_over_exner inside exner
      m_elem_ops.compute_midpoint_delta(kv,phi_i,exner);

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {

        // TODO: could do *= Rgas/p0, but may lose BFB with F90.
        pnh(ilev) = PhysicalConstants::Rgas*vtheta_dp(ilev) / (-exner(ilev));
        pnh(ilev) /= PhysicalConstants::p0;
        pow_update(pnh(ilev),1.0/(1.0-PhysicalConstants::kappa));
        pnh(ilev) *= PhysicalConstants::p0;

        exner(ilev) = pnh(ilev)/exner(ilev);
      });
    }
  }

  KOKKOS_INLINE_FUNCTION
  void compute_dpnh_dp_i (KernelVariables& kv,
                          ExecViewUnmanaged<const Scalar[NUM_LEV  ]> pnh,
                          ExecViewUnmanaged<const Scalar[NUM_LEV_P]> dp_i,
                          ExecViewUnmanaged<      Scalar[NUM_LEV_P]> dpnh_dp_i) const
  {
    if (m_theta_hydrostatic_mode) {
      // Set dpnh_dp_i to 1.0
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV_P),
                           [&](const int ilev) {
        dpnh_dp_i(ilev) = 1.0;
      });
    } else {
      // Start with dpnh_dp_i = delta(pnh)/dp_i
      m_elem_ops.compute_interface_delta(kv.team,pnh,dpnh_dp_i);

      // Note: top and bottom need special treatment, so we may as well stop at NUM_LEV here
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        dpnh_dp_i(ilev) /= dp_i(ilev);
      });

      // Note: the last midpoint is in the pack NUM_LEV, at position (NUM_PHYSICAL_LEV-1)%VECTOR_SIZE;
      //       the last interface is in the pack NUM_LEV_P, one position after the last midpoint (mod VECTOR_SIZE).
      //       This is true regardless of whether NUM_LEV_P>NUM_LEV or not, so one formula works for all cases
      constexpr int last_midpoint_idx = (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE; 
      constexpr int last_interface_idx = (last_midpoint_idx+1) % VECTOR_SIZE;

      // Boundaries: delta(x) = 2*(x_m(last)-x_i(last)).
      // Top: pnh_i = pi_i.
      // Bottom: approximate with hydrostatic, so that dpnh_dp_i=1
      dpnh_dp_i(0)[0] = 2*(pnh(0)[0] - m_hvcoord.hybrid_ai(0)*m_hvcoord.ps0)/dp_i(0)[0];
      const Real pnh_last = pnh(NUM_LEV-1)[last_midpoint_idx];
      const Real dp_last = dp_i(NUM_LEV_P-1)[last_interface_idx];
      const Real pnh_i_last = pnh_last + dp_last;
      dpnh_dp_i(NUM_LEV_P-1)[last_interface_idx] = 2*(pnh_i_last - pnh_last)/dp_last;
    }
  }

private:

  bool            m_theta_hydrostatic_mode;
  ElementOps      m_elem_ops;
  HybridVCoord    m_hvcoord;
};

} // namespace Homme

#endif // HOMMEXX_EQUATION_OF_STATE_HPP
