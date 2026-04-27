#ifndef GW_VD_LU_SOLVE_IMPL_HPP
#define GW_VD_LU_SOLVE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw vd_lu_solve. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_lu_solve(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver,
  const uview_1d<const Real>& decomp_ca,
  const uview_1d<const Real>& decomp_cc,
  const uview_1d<const Real>& decomp_dnom,
  const uview_1d<const Real>& decomp_ze,
  const Int& ntop,
  const Int& nbot,
  const Real& cd_top,
  // Inputs/Outputs
  const uview_1d<Real>& q)
{
  auto zf = workspace.take("zf");





  int nan_count = 0, inf_count = 0;

  






  // All the loops have dependencies on prior iterations of the loop, no parallelism possible here
  Kokkos::single(Kokkos::PerTeam(team), [&] {



    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(q(k))) nan_count++;
      if (Kokkos::isinf(q(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] q (in) NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(decomp_dnom(k))) nan_count++;
      if (Kokkos::isinf(decomp_dnom(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] decomp_dnom NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(decomp_ca(k))) nan_count++;
      if (Kokkos::isinf(decomp_ca(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] decomp_ca NaN=%d Inf=%d\n", nan_count, inf_count);
    }

    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(decomp_ze(k))) nan_count++;
      if (Kokkos::isinf(decomp_ze(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] decomp_ze NaN=%d Inf=%d\n", nan_count, inf_count);
    }



    for (int k = 0; k < pver; k++) {
      zf(k) = 0;
    }


    // Calculate zf(k). Terms zf(k) and ze(k) are required in solution of
    // tridiagonal matrix defined by implicit diffusion equation.
    // Note that only levels ntop through nbot need be solved for.
    zf(nbot) = q(nbot)*decomp_dnom(nbot);
    for (int k = nbot - 1; k >= ntop + 1; --k) {
      zf(k) = (q(k) + decomp_ca(k)*zf(k+1)) * decomp_dnom(k);
    }

    // Include boundary condition on top element
    zf(ntop) = (q(ntop) + cd_top + decomp_ca(ntop)*zf(ntop+1)) * decomp_dnom(ntop);



    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(zf(k))) nan_count++;
      if (Kokkos::isinf(zf(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] zf NaN=%d Inf=%d\n", nan_count, inf_count);
    }





    // Perform back substitution
    q(ntop) = zf(ntop);
    for (int k = ntop + 1; k <= nbot; ++k) {
      q(k) = zf(k) + decomp_ze(k)*q(k-1);
    }


    nan_count = 0; inf_count = 0;
    for (int k = 0; k < pver; k++) {
      if (Kokkos::isnan(q(k))) nan_count++;
      if (Kokkos::isinf(q(k))) inf_count++;
    }
    if (nan_count > 0 || inf_count > 0) {
      printf("[vd_lu_solve] q (out) NaN=%d Inf=%d\n", nan_count, inf_count);
    }



  });







  workspace.release(zf);
}

} // namespace gw
} // namespace scream

#endif
