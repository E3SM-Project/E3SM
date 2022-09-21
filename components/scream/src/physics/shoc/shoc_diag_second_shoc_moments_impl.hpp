#ifndef SHOC_DIAG_SECOND_SHOC_MOMENTS_IMPL_HPP
#define SHOC_DIAG_SECOND_SHOC_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_second_shoc_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_second_shoc_moments(const MemberType& team, const Int& nlev, const Int& nlevi,
       const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind,
       const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
       const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi,
       const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix,
       const Scalar& wthl_sfc, const Scalar& wqw_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc, Scalar& ustar2, Scalar& wstar,
       const Workspace& workspace, const uview_1d<Spack>& thl_sec,
       const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec, const uview_1d<Spack>& qwthl_sec,
       const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec, const uview_1d<Spack>& w_sec)
{
  // This is the main routine to compute the second
  // order moments in SHOC.

  // Define temporary variables
  uview_1d<Spack> isotropy_zi, tkh_zi, tk_zi;
  workspace.template take_many_contiguous_unsafe<3>(
    {"isotropy_zi", "tkh_zi", "tk_zi"},
    {&isotropy_zi, &tkh_zi, &tk_zi});

  // Calculate surface properties needed for lower
  //  boundary conditions
  shoc_diag_second_moments_srf(wthl_sfc, uw_sfc, vw_sfc, ustar2, wstar);
  team.team_barrier();

  // Scalarize views for single entry access
  const auto s_wthl_sec  = ekat::scalarize(wthl_sec);
  const auto s_wqw_sec   = ekat::scalarize(wqw_sec);
  const auto s_uw_sec    = ekat::scalarize(uw_sec);
  const auto s_vw_sec    = ekat::scalarize(vw_sec);
  const auto s_wtke_sec  = ekat::scalarize(wtke_sec);
  const auto s_thl_sec   = ekat::scalarize(thl_sec);
  const auto s_qw_sec    = ekat::scalarize(qw_sec);
  const auto s_qwthl_sec = ekat::scalarize(qwthl_sec);

  // Diagnose the second order moments flux, for the lower boundary
  shoc_diag_second_moments_lbycond(wthl_sfc, wqw_sfc, uw_sfc, vw_sfc, ustar2, wstar,
                                   s_wthl_sec(nlevi-1), s_wqw_sec(nlevi-1),
                                   s_uw_sec(nlevi-1), s_vw_sec(nlevi-1),
                                   s_wtke_sec(nlevi-1), s_thl_sec(nlevi-1),
                                   s_qw_sec(nlevi-1), s_qwthl_sec(nlevi-1));
  team.team_barrier();

  // Diagnose the second order moments, for points away from boundaries.  this is
  //  the main computation for the second moments
  diag_second_moments(team, nlev, nlevi,
                     thetal, qw, u_wind,v_wind, tke, isotropy,tkh, tk, dz_zi, zt_grid, zi_grid, shoc_mix,
                     isotropy_zi, tkh_zi, tk_zi, thl_sec, qw_sec, wthl_sec, wqw_sec,
                     qwthl_sec, uw_sec, vw_sec, wtke_sec, w_sec);
  team.team_barrier();

  // Diagnose the second order moments, calculate the upper boundary conditions
  shoc_diag_second_moments_ubycond(s_thl_sec(0), s_qw_sec(0), s_wthl_sec(0),
                                   s_wqw_sec(0), s_qwthl_sec(0), s_uw_sec(0),
                                   s_vw_sec(0), s_wtke_sec(0));

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<3>(
    {&isotropy_zi, &tkh_zi, &tk_zi});
}

} // namespace shoc
} // namespace scream

#endif
