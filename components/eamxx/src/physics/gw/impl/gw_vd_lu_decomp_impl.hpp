#ifndef GW_VD_LU_DECOMP_IMPL_HPP
#define GW_VD_LU_DECOMP_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw vd_lu_decomp. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_lu_decomp(
  // Inputs
  const MemberType& team,
  const Int& pver,
  const Real& ksrf,
  const uview_1d<const Real>& kv,
  const uview_1d<const Real>& tmpi,
  const uview_1d<const Real>& rpdel,
  const Real& ztodt,
  const Real& cc_top,
  const Int& ntop,
  const Int& nbot,
  // Outputs
  const uview_1d<Real>& decomp_ca,
  const uview_1d<Real>& decomp_cc,
  const uview_1d<Real>& decomp_dnom,
  const uview_1d<Real>& decomp_ze)
{
  // Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the
  // tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
  // a combination of ca and cc; they are not required by the solver.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, ntop, nbot), [&] (const int k) {
      decomp_ca(k  ) = kv(k+1) * tmpi(k+1) * rpdel(k  );
      decomp_cc(k+1) = kv(k+1) * tmpi(k+1) * rpdel(k+1);
      // The bottom element of the upper diagonal (ca) is zero (not used).
      // The subdiagonal (cc) is not needed in the solver.
      if (k == nbot - 1) {
        decomp_ca(nbot) = 0;
      }
    });

  // Calculate e(k). This term is required in the solution of the
  // tridiagonal matrix as defined by the implicit diffusion equation.
  team.team_barrier();
  const Real dnom = 1 / (1 + decomp_cc(nbot) + ksrf*ztodt*C::gravit*rpdel(nbot));
  decomp_dnom(nbot) = dnom;
  decomp_ze(nbot)   = decomp_cc(nbot) * dnom;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, ntop+1, nbot), [&] (const int k) {
      decomp_dnom(k) = 1 / (1 + decomp_ca(k) + decomp_cc(k) -
                            decomp_ca(k)*decomp_ze(k+1));
      decomp_ze(k)   = decomp_cc(k) * decomp_dnom(k);
    });

  team.team_barrier();
  decomp_dnom(ntop) = 1 / (1 + decomp_ca(ntop) + cc_top - decomp_ca(ntop)*decomp_ze(ntop+1));
}

} // namespace gw
} // namespace scream

#endif
