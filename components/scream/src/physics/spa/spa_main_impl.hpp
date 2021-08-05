#ifndef SPA_MAIN_IMPL_HPP
#define SPA_MAIN_IMPL_HPP

#include "physics/spa/spa_functions.hpp"
#include "physics/share/physics_constants.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
:: reconstruct_pressure_profile(
  const Int ncols,
  const Int nlevs,
  const view_1d<const Spack>& hya,
  const view_1d<const Spack>& hyb,
  const view_1d<const Real>&  PS,
  const view_2d<Spack>&       pres)
{
  using C = scream::physics::Constants<Real>;
  static constexpr auto P0 = C::P0;
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nlevs);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols, nk_pack);
  Kokkos::parallel_for(
    "reconstruct SPA pressure profile",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    const auto pres_sub = ekat::subview(pres, i);
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        pres_sub(k) = PS(i) * hyb(k) + P0 * hya(k);
      });
    });

} // reconstruct_pressure_profile
/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::aero_vertical_remap(
    const Int ncols,
    const Int nlevs_src,
    const Int nlevs_tgt,
    const view_2d<const Spack>& pres_src,
    const view_2d<const Spack>& pres_tgt,
    const view_2d<const Spack>& aero_src,
    const view_2d<Spack>& aero_tgt) 
{
  using LIV = ekat::LinInterp<Real,Spack::n>;
  Real minthreshold = 0.0;  // Hard-code a minimum value for CCN to zero.
  LIV VertInterp(ncols,nlevs_src,nlevs_tgt,minthreshold); 
  {
    Kokkos::parallel_for("vertical-interp-spa",
      VertInterp.m_policy,
      KOKKOS_LAMBDA(typename LIV::MemberType const& team_member) {
        const int i = team_member.league_rank();
        VertInterp.setup(team_member,
                         ekat::subview(pres_src,i),
                         ekat::subview(pres_tgt,i));
        team_member.team_barrier();
        VertInterp.lin_interp(team_member,
                              ekat::subview(pres_src,i),
                              ekat::subview(pres_tgt,i),
                              ekat::subview(aero_src,i),
                              ekat::subview(aero_tgt,i));
      });
  }

} // aero_vertical_remap
/*-----------------------------------------------------------------*/
template <typename S, typename D>
void SPAFunctions<S,D>
::aero_time_interp(
    const Real& t0,
    const Real& ts,
    const Real& tlen,
    const Real& y0,
    const Real& y1,
          Real& y_out)
{
 // Simple linear interpolation: y_out = b + m*x, b = y0, m = (y1-y0)/tlen
 y_out = y0 + (ts-t0) * (y1-y0)/tlen;
}
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_MAIN_IMPL_HPP
