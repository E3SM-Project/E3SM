#ifndef GW_GWD_PRECALC_RHOI_IMPL_HPP
#define GW_GWD_PRECALC_RHOI_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_precalc_rhoi. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_precalc_rhoi(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver,
  const Int& pgwv,
  const Real& dt,
  const Int& tend_level,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& t,
  const uview_2d<const Real>& gwut,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& nm,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& c,
  const uview_2d<const Real>& q,
  const uview_1d<const Real>& dse,
  // Outputs
  const uview_1d<Real>& egwdffi,
  const uview_2d<Real>& qtgw,
  const uview_1d<Real>& dttdf,
  const uview_1d<Real>& dttke,
  const uview_1d<Real>& ttgw)
{
#if 0
  // rhoi_kludge: Recalculated rhoi to preserve answers.
  uview_1d<Real> rhoi_kludge, decomp_ca_1d, decomp_cc_1d, decomp_dnom_1d, decomp_ze_1d;
  workspace2.template take_many_contiguous_unsafe<5>(
    {"rhoi_kludge", "decomp_ca", "decomp_cc", "decomp_dnom", "decomp_ze"},
    {&rhoi_kludge, &decomp_ca_1d, &decomp_cc_1d, &decomp_dnom_1d, &decomp_ze_1d});

  // Note: pint is from 0 to pver instead of 1 to pver+1.
  rhoi_kludge(:,1) = pint(:,0) / (rair * t(:,1))
  do k = 2, pver
     rhoi_kludge(:,k) = pint(:,k-1) * 2._r8 / &
          (rair * (t(:,k) + t(:,k-1)))
  end do
  rhoi_kludge(:,pver+1) = pint(:,pver) / (rair * t(:,pver))

  // Calculate effective diffusivity and LU decomposition for the
  // vertical diffusion solver.
  call gw_ediff (ncol, pver, ngwv, kbotbg, ktop, tend_level, &
       gwut, ubm, nm, rhoi_kludge, dt, gravit, pmid, rdpm, c, &
       egwdffi, decomp)

  // Calculate tendency on each constituent.
  do m = 1, size(q,3)
     call gw_diff_tend(ncol, pver, kbotbg, ktop, q(:,:,m), dt, &
          decomp, qtgw(:,:,m))
  enddo

  // Calculate tendency from diffusing dry static energy (dttdf).
  call gw_diff_tend(ncol, pver, kbotbg, ktop, dse, dt, decomp, dttdf)

  // Evaluate second temperature tendency term: Conversion of kinetic
  // energy into thermal.
  do l = -ngwv, ngwv
     do k = ktop+1, kbotbg
        dttke(:,k) = dttke(:,k) + c(:,l) * gwut(:,k,l)
     end do
  end do

  ttgw = dttke + dttdf

  workspace.template release_many_contiguous<5>(
    {&rhoi_kludge, &decomp_ca_1d, &decomp_cc_1d, &decomp_dnom_1d, &decomp_ze_1d});
#endif
}

} // namespace gw
} // namespace scream

#endif
