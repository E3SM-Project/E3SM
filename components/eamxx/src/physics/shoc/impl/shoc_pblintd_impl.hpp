#ifndef SHOC_PBLINTD_IMPL_HPP
#define SHOC_PBLINTD_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc pblintd. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

//-----------------------------------------------------------------------
//
// Purpose:
// Diagnose standard PBL variables
//
// Method:
// Diagnosis of PBL depth.
// The PBL depth follows:
//    Holtslag, A.A.M., and B.A. Boville, 1993:
//    Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
//    Model. J. Clim., vol. 6., p. 1825--1842.
//
// Updated by Holtslag and Hack to exclude the surface layer from the
// definition of the boundary layer Richardson number. Ri is now defined
// across the outer layer of the pbl (between the top of the surface
// layer and the pbl top) instead of the full pbl (between the surface and
// the pbl top). For simiplicity, the surface layer is assumed to be the
// region below the first model level (otherwise the boundary layer depth
// determination would require iteration).
//
// Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
// (Use ricr = 0.3 in this formulation)
//
// Author: B. Stevens (extracted from pbldiff, August 2000)
//

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::pblintd(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   npbl,
  const uview_1d<const Spack>& z,
  const uview_1d<const Spack>& zi,
  const uview_1d<const Spack>& thl,
  const uview_1d<const Spack>& ql,
  const uview_1d<const Spack>& q,
  const uview_1d<const Spack>& u,
  const uview_1d<const Spack>& v,
  const Scalar&                ustar,
  const Scalar&                obklen,
  const Scalar&                kbfs,
  const uview_1d<const Spack>& cldn,
  const Workspace&             workspace,
  Scalar&                      pblh)
{
  // Define temporary variables
  uview_1d<Spack> rino, thv;
  workspace.template take_many_contiguous_unsafe<2>(
    {"rino", "thv"},
    {&rino, &thv});

  // Scalarize views for single entry access
  const auto s_z    = ekat::scalarize(z);
  const auto s_thv  = ekat::scalarize(thv);
  const auto s_rino = ekat::scalarize(rino);
  const auto s_zi   = ekat::scalarize(zi);
  const auto s_cldn = ekat::scalarize(cldn);

  // Compute Obukhov length virtual temperature flux and various arrays for use later:

  //Compute virtual potential temperature
  shoc_pblintd_init_pot(team,nlev,thl,ql,q,thv);

  // Initialize
  bool check = true;
  // The loop below fixes valgrind uninitialized mem errs. As in other
  // places in eamxx, we use SCREAM_SHORT_TESTS as a proxy for knowing
  // mem checking is on.
#if !defined(NDEBUG) || defined(SCREAM_SHORT_TESTS)
  for (size_t i=0; i<rino.size(); ++i) {
    rino(i)=0;
  }
#else
  s_rino(nlev-1) = 0;
#endif
  pblh = s_z(nlev-1);

  // PBL height calculation
  team.team_barrier();
  pblintd_height(team,nlev,npbl,z,u,v,ustar,
                 thv,s_thv(nlev-1),
                 pblh,rino,check);

  // Estimate an effective surface temperature to account for surface fluctuations
  Scalar tlv;
  pblintd_surf_temp(nlev,nlevi,npbl,z,ustar,obklen,kbfs,thv,tlv,pblh,check,rino);

  // Improve pblh estimate for unstable conditions using the convective
  // temperature excess as reference temperature:
  team.team_barrier();
  pblintd_height(team,nlev,npbl,z,u,v,ustar,
                 thv,tlv,
                 pblh,rino,check);

  // Check PBL height
  pblintd_check_pblh(nlevi,npbl,z,ustar,check,pblh);

  // PBL check over ocean
  shoc_pblintd_cldcheck(s_zi(nlev-1),s_cldn(nlev-1),pblh);

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<2>(
    {&rino, &thv});
}

} // namespace shoc
} // namespace scream

#endif
