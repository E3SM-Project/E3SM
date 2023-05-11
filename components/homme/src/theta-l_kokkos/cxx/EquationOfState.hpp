#ifndef HOMMEXX_EQUATION_OF_STATE_HPP
#define HOMMEXX_EQUATION_OF_STATE_HPP

#include "Types.hpp"
#include "ColumnOps.hpp"
#include "HybridVCoord.hpp"
#include "KernelVariables.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/VectorUtils.hpp"
#include "utilities/ViewUtils.hpp"
#include "utilities/BfbUtils.hpp"

namespace Homme {

class EquationOfState {
public:

  enum : int { LAST_MID_PACK     = ColInfo<NUM_PHYSICAL_LEV>::LastPack     };
  enum : int { LAST_MID_PACK_END = ColInfo<NUM_PHYSICAL_LEV>::LastPackEnd  };
  enum : int { LAST_INT_PACK     = ColInfo<NUM_INTERFACE_LEV>::LastPack    };
  enum : int { LAST_INT_PACK_END = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd };

  EquationOfState () = default;

  void init (const bool theta_hydrostatic_mode,
             const HybridVCoord& hvcoord) {
    m_theta_hydrostatic_mode = theta_hydrostatic_mode;
    m_hvcoord = hvcoord;
    assert (m_hvcoord.m_inited);
  }

  // On input, pe is pressure; on output, the Exner function.
  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION
  static void pressure_to_exner (Scalar& pe) {
    pe /= PhysicalConstants::p0;
#ifdef HOMMEXX_BFB_TESTING
    pe = bfb_pow(pe,PhysicalConstants::kappa);
#else
    pe = pow(pe,PhysicalConstants::kappa);
#endif
  }

  // On input, pe is pressure; on output, the reciprocal of the Exner function.
  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION
  static void pressure_to_recip_exner (Scalar& pe) {
    pe = PhysicalConstants::p0 / pe;
#ifdef HOMMEXX_BFB_TESTING
    pe = bfb_pow(pe,PhysicalConstants::kappa);
#else
    pe = pow(pe,PhysicalConstants::kappa);
#endif
  }

  KOKKOS_INLINE_FUNCTION
  void compute_exner (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Scalar[NUM_LEV]>& pi,
                      const ExecViewUnmanaged<      Scalar[NUM_LEV]>& exner) const
  {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                         [&](const int ilev) {
      // Avoid temporaries
      exner(ilev) = pi(ilev);
      pressure_to_exner(exner(ilev));
    });
  }

  template<typename VThetaProvider, typename PhiProvider>
  KOKKOS_INLINE_FUNCTION
  bool compute_pnh_and_exner (const KernelVariables& kv,
                              const VThetaProvider& vtheta_dp,
                              const PhiProvider&    phi_i,
                              const ExecViewUnmanaged<Scalar[NUM_LEV]>& pnh,
                              const ExecViewUnmanaged<Scalar[NUM_LEV]>& exner) const
  {
    // If you're hydrostatic, check outside the function
    assert (!m_theta_hydrostatic_mode);

    // To avoid temporaries, use exner to store some temporaries
    ColumnOps::compute_midpoint_delta(kv,phi_i,exner);

    int nerr;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                            [&](const int ilev, int& nerr) {
      // check inputs
      const int vec_len = (ilev==(NUM_LEV-1)) ? ColInfo<NUM_PHYSICAL_LEV>::LastPackLen : VECTOR_SIZE;
      for (int iv=0; iv<vec_len; ++iv) {
        if (vtheta_dp(ilev)[iv] < 0.0 || exner(ilev)[iv] > 0.0)
          ++nerr;
      }
      if (nerr) return;
      compute_pnh_and_exner(vtheta_dp(ilev), exner(ilev), pnh(ilev), exner(ilev));
    }, nerr);
    return nerr == 0;
  }

  // Compute:
  //  1) p_over_exner = -Rgas*vtheta_dp/delta(phi_i)
  //  2) pnh = p0 (p_over_exner/p0)^(1/(1-kappa))
  //  3) exner = pnh/p_over_exner
  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION
  static void compute_pnh_and_exner (const Scalar& vtheta_dp, const Scalar& dphi,
                                     Scalar& pnh, Scalar& exner) {
    exner = (-PhysicalConstants::Rgas)*vtheta_dp / dphi;
    pnh = exner/PhysicalConstants::p0;
#ifndef HOMMEXX_BFB_TESTING
    pnh = pow(pnh,1.0/(1.0-PhysicalConstants::kappa));
#else
    pnh = bfb_pow(pnh,1.0/(1.0-PhysicalConstants::kappa));
#endif
    pnh *= PhysicalConstants::p0;
    exner = pnh/exner;    
  }

  KOKKOS_INLINE_FUNCTION
  void compute_dpnh_dp_i (const KernelVariables& kv,
                          const ExecViewUnmanaged<const Scalar[NUM_LEV  ]>& pnh,
                          const ExecViewUnmanaged<const Scalar[NUM_LEV_P]>& dp_i,
                          const ExecViewUnmanaged<      Scalar[NUM_LEV_P]>& dpnh_dp_i) const
  {
    if (m_theta_hydrostatic_mode) {
      // Set dpnh_dp_i to 1.0
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV_P),
                           [&](const int ilev) {
        dpnh_dp_i(ilev) = 1.0;
      });
    } else {
      // dp_i CANNOT alias dpnh_dp_i
      assert(dpnh_dp_i.data()!=dp_i.data());

      // Start with dpnh_dp_i = delta(pnh)/dp_i.
      ColumnOps::compute_interface_delta<CombineMode::Replace>(kv.team,pnh,dpnh_dp_i);

      // Note: top and bottom need special treatment, so we may as well stop at NUM_LEV here (rather than NUM_LEV_P)
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                           [&](const int ilev) {
        dpnh_dp_i(ilev) /= dp_i(ilev);
      });

      // Boundaries: delta(x) = 2*(x_m(last)-x_i(last)).
      // Top: pnh_i = pi_i = hyai(0)*ps0.
      // Bottom: approximate with hydrostatic, so that dpnh_dp_i=1
      dpnh_dp_i(0)[0] = 2*(pnh(0)[0] - m_hvcoord.hybrid_ai(0)*m_hvcoord.ps0)/dp_i(0)[0];
      const Real pnh_last = pnh(LAST_MID_PACK)[LAST_MID_PACK_END];
      const Real dp_last = dp_i(LAST_INT_PACK)[LAST_INT_PACK_END];
      const Real pnh_i_last = pnh_last + dp_last/2;

      dpnh_dp_i(LAST_INT_PACK)[LAST_INT_PACK_END] = 2*(pnh_i_last - pnh_last)/dp_last;
    }
  }

  // Note: if p is hydrostatic, this will compute the hydrostatic geopotential,
  //       otherwise it will be the non-hydrostatic. In particular, if the pressure
  //       p is computed using dp from pnh, this will be the discrete inverse of
  //       the compute_pnh_and_exner method.
  KOKKOS_INLINE_FUNCTION static
  void compute_phi_i (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Real   [NP][NP]           >& phis,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& p,
                      const ExecViewUnmanaged<      Scalar [NP][NP][NUM_LEV_P]>& phi_i) {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      compute_phi_i(kv, phis(igp,jgp),
                    Homme::subview(vtheta_dp,igp,jgp),
                    Homme::subview(p,igp,jgp),
                    Homme::subview(phi_i,igp,jgp));
    });
  }

  // VThetaProvider can be either a 1d view or a lambda,
  // as long as vtheta_dp(ilev) returns vtheta_dp at pack ilev
  template<typename VThetaProvider>
  KOKKOS_INLINE_FUNCTION static
  void compute_phi_i (const KernelVariables& kv, const Real phis,
                      const VThetaProvider& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& p,
                      const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i)
  {
    // Init phi on surface with phis
    phi_i(LAST_INT_PACK)[LAST_INT_PACK_END] = phis;

    // Use ColumnOps to do the scan sum
    auto integrand_provider = [&](const int ilev)->Scalar {
      return compute_dphi(vtheta_dp(ilev), p(ilev));
    };

    ColumnOps::column_scan_mid_to_int<false>(kv,integrand_provider,phi_i);
  }

  template<typename Scalar>
  KOKKOS_INLINE_FUNCTION static
  Scalar compute_dphi (const Scalar& vtheta_dp, const Scalar& p) {
    constexpr Real p0    = PhysicalConstants::p0;
    constexpr Real kappa = PhysicalConstants::kappa;
    constexpr Real Rgas  = PhysicalConstants::Rgas;
#ifdef HOMMEXX_BFB_TESTING
    return (Rgas*vtheta_dp * bfb_pow(p/p0,kappa-1)) / p0;
#else
    // TODO: remove temporaries
    return (Rgas*vtheta_dp * pow(p/p0,kappa-1)) / p0;
#endif    
  }

  // If exner is available, then use exner/p instead of (p/p0)^(k-1)/p0, to avoid dealing with exponentials
  // VThetaProvider can be either a 1d view or a lambda,
  // as long as vtheta_dp(ilev) returns vtheta_dp at pack ilev
  KOKKOS_INLINE_FUNCTION
  void compute_phi_i (const KernelVariables& kv,
                      const ExecViewUnmanaged<const Real   [NP][NP]           >& phis,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NP][NP][NUM_LEV]  >& p,
                      const ExecViewUnmanaged<      Scalar [NP][NP][NUM_LEV_P]>& phi_i) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;
      compute_phi_i(kv, phis(igp,jgp),
                    Homme::subview(vtheta_dp,igp,jgp),
                    Homme::subview(exner,igp,jgp),
                    Homme::subview(p,igp,jgp),
                    Homme::subview(phi_i,igp,jgp));
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_phi_i (const KernelVariables& kv, const Real phis,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& vtheta_dp,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& exner,
                      const ExecViewUnmanaged<const Scalar [NUM_LEV]  >& p,
                      const ExecViewUnmanaged<      Scalar [NUM_LEV_P]>& phi_i) const
  {
    // Init phi on surface with phis
    phi_i(LAST_INT_PACK)[LAST_INT_PACK_END] = phis;

    // Use ColumnOps to do the scan sum
    auto integrand_provider = [&](const int ilev)->Scalar {
      constexpr Real Rgas  = PhysicalConstants::Rgas;
      return Rgas*vtheta_dp(ilev)*exner(ilev)/p(ilev);
    };

    ColumnOps::column_scan_mid_to_int<false>(kv,integrand_provider,phi_i);
  }

public:

  bool            m_theta_hydrostatic_mode;
  HybridVCoord    m_hvcoord;
};

} // namespace Homme

#endif // HOMMEXX_EQUATION_OF_STATE_HPP
