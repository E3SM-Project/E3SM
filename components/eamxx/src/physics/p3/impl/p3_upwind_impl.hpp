#ifndef P3_UPWIND_IMPL_HPP
#define P3_UPWIND_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of upwind functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
template <Int kdir, int nfield>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_first_order_upwind_step (
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& inv_dz,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Scalar& dt_sub,
  const view_1d_ptr_array<Spack, nfield>& flux,
  const view_1d_ptr_array<Spack, nfield>& V,
  const view_1d_ptr_array<Spack, nfield>& r)
{
  const Int kmin_scalar = ( kdir == 1 ? k_bot : k_top);
  const Int kmax_scalar = ( kdir == 1 ? k_top : k_bot);
  Int
    kmin = kmin_scalar / Spack::n,
    // Add 1 to make [kmin, kmax). But then the extra term (Spack::n -
    // 1) to determine pack index cancels the +1.
    kmax = (kmax_scalar + Spack::n) / Spack::n;

  // The following is morally a const var, but there are issues with
  // gnu and std=c++14. The macro ConstExceptGnu is defined in ekat_kokkos_types.hpp.
  ConstExceptGnu Int k_top_pack = k_top / Spack::n;

  // calculate fluxes
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, kmax - kmin), [&] (Int k_) {
      const Int k = kmin + k_;
      for (int f = 0; f < nfield; ++f)
        (*flux[f])(k) = (*V[f])(k) * (*r[f])(k) * rho(k);
    });
  team.team_barrier();

  Kokkos::single(
    Kokkos::PerTeam(team), [&] () {
      const Int k = k_top_pack;
      {
        const auto range_pack = ekat::range<IntSmallPack>(k_top_pack*Spack::n);
        const auto mask = range_pack > kmax_scalar || range_pack < kmin_scalar;
        if (mask.any()) {
          for (int f = 0; f < nfield; ++f) {
            (*flux[f])(k_top_pack).set(mask, 0);
          }
        }
      }
      for (int f = 0; f < nfield; ++f) {
        // compute flux divergence
        const auto flux_pkdir = (kdir == -1) ?
          shift_right(0, (*flux[f])(k)) :
          shift_left (0, (*flux[f])(k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(k)) * inv_dz(k);

        // update prognostic variables
        (*r[f])(k) += fluxdiv * dt_sub * inv_rho(k);
      }
    });

  if (kdir == 1)
    --kmax;
  else
    ++kmin;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, kmax - kmin), [&] (Int k_) {
      const Int k = kmin + k_;
      for (int f = 0; f < nfield; ++f) {
        // compute flux divergence
        const auto flux_pkdir = (kdir == -1) ?
          shift_right((*flux[f])(k+kdir), (*flux[f])(k)) :
          shift_left ((*flux[f])(k+kdir), (*flux[f])(k));
        const auto fluxdiv = (flux_pkdir - (*flux[f])(k)) * inv_dz(k);
        // update prognostic variables
        (*r[f])(k) += fluxdiv * dt_sub * inv_rho(k);
      }
    });
}

template <typename S, typename D>
template <int nfield>
KOKKOS_FUNCTION
void Functions<S,D>
::generalized_sedimentation (
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& inv_dz,
  const MemberType& team,
  const Int& nk, const Int& k_qxtop, Int& k_qxbot, const Int& kbot, const Int& kdir, const Scalar& Co_max, Scalar& dt_left, Scalar& prt_accum,
  const view_1d_ptr_array<Spack, nfield>& fluxes,
  const view_1d_ptr_array<Spack, nfield>& Vs, // (behaviorally const)
  const view_1d_ptr_array<Spack, nfield>& rs)
{
  // compute dt_sub
  EKAT_KERNEL_ASSERT(Co_max >= 0);
  const Int tmpint1 = static_cast<int>(Co_max + 1);
  const Scalar dt_sub = dt_left/tmpint1;

  // Move bottom cell down by 1 if not at ground already
  const Int k_temp = (k_qxbot == kbot) ? k_qxbot : k_qxbot - kdir;

  calc_first_order_upwind_step<nfield>(rho, inv_rho, inv_dz, team, nk, k_temp, k_qxtop, kdir, dt_sub, fluxes, Vs, rs);
  team.team_barrier();

  // accumulated precip during time step
  if (k_qxbot == kbot) {
    const auto sflux0 = scalarize(*fluxes[0]);
    prt_accum += sflux0(kbot) * dt_sub;
  }
  else {
    k_qxbot -= kdir;
  }

  // update time remaining for sedimentation
  dt_left -= dt_sub;
}

template <typename S, typename D>
template <int nfield>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_first_order_upwind_step (
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& inv_dz,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
  const view_1d_ptr_array<Spack, nfield>& flux,
  const view_1d_ptr_array<Spack, nfield>& V,
  const view_1d_ptr_array<Spack, nfield>& r)
{
  if (kdir == 1)
    calc_first_order_upwind_step< 1, nfield>(
      rho, inv_rho, inv_dz, team, nk, k_bot, k_top, dt_sub, flux, V, r);
  else
    calc_first_order_upwind_step<-1, nfield>(
      rho, inv_rho, inv_dz, team, nk, k_bot, k_top, dt_sub, flux, V, r);
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_first_order_upwind_step (
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& inv_dz,
  const MemberType& team,
  const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, const Scalar& dt_sub,
  const uview_1d<Spack>& flux,
  const uview_1d<const Spack>& V,
  const uview_1d<Spack>& r)
{
  // B/c automatic casting to const does not work in the nested data
  // view_1d_ptr_array (C++ does not provide all legal const casts automatically
  // in nested data structures), we are not enforcing const in the array
  // versions of this function, as doing so would be a burden to the caller. But
  // in this version, we can. Thus, to call through to the impl, we explicitly
  // cast here.
  const auto V_nonconst = uview_1d<Spack>(const_cast<Spack*>(V.data()),
                                          V.extent_int(0));
  if (kdir == 1)
    calc_first_order_upwind_step< 1, 1>(
      rho, inv_rho, inv_dz, team, nk, k_bot, k_top, dt_sub,
      {&flux}, {&V_nonconst}, {&r});
  else
    calc_first_order_upwind_step<-1, 1>(
      rho, inv_rho, inv_dz, team, nk, k_bot, k_top, dt_sub,
      {&flux}, {&V_nonconst}, {&r});
}

} // namespace p3
} // namespace scream

#endif
