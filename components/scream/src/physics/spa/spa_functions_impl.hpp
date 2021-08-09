#ifndef SPA_FUNCTIONS_IMPL_HPP
#define SPA_FUNCTIONS_IMPL_HPP

#include "physics/share/physics_constants.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

namespace scream {
namespace spa {

/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
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
KOKKOS_INLINE_FUNCTION
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
template <typename ScalarY>
KOKKOS_INLINE_FUNCTION
ScalarY SPAFunctions<S,D>::aero_time_interp(
  const Real& t0,
  const Real& ts,
  const Real& tlen,
  const ScalarY& y0,
  const ScalarY& y1)
{
  // Simple linear interpolation: y_out = b + m*x, b = y0, m = (y1-y0)/tlen
  return y0 + (ts-t0) * (y1-y0)/tlen;
} // aero_time_interp
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void SPAFunctions<S,D>::
aero_time_interp(
  const MemberType& team,
  const Real& t0,
  const Real& ts,
  const Real& tlen,
  const uview_1d<const Spack>& y0,
  const uview_1d<const Spack>& y1,
  const uview_1d<Spack>&       y_out)
{
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,y_out.extent(0)),
                       [&] (const int k) {
    y_out(k) = aero_time_interp(t0, ts, tlen, y0(k), y1(k));
  });
} // aero_time_interp
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void SPAFunctions<S,D>::
calculate_current_data_ps(
  const Int&  ncols,
  const Real& t0,
  const Real& ts,
  const Real& tlen,
  const view_1d<const Real>& ps_0,
  const view_1d<const Real>& ps_1,
  const view_1d<Real>&       ps_out)
{
  using ExeSpace = typename KT::ExeSpace;
  Kokkos::parallel_for("time-interp-PS",
    ncols,
    KOKKOS_LAMBDA(const int& i) {
      ps_out(i) = aero_time_interp(t0, ts, tlen, ps_0(i), ps_1(i));
    });
} // calculate_current_data_ps
/*-----------------------------------------------------------------*/
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void SPAFunctions<S,D>::
calculate_current_data_concentrations(
  const Int& ncols,
  const Int& nlevs,
  const Real& t0,
  const Real& ts,
  const Real& tlen,
  const view_2d<const Spack>& y0,
  const view_2d<const Spack>& y1,
  const view_2d<Spack>&       yout)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nlevs);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncols, nk_pack);
  Kokkos::parallel_for(
    "time-interp-aero",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

      const Int i = team.league_rank();
      const auto& y0_sub   = ekat::subview(y0,i);
      const auto& y1_sub   = ekat::subview(y1,i);
      const auto& yout_sub = ekat::subview(yout,i);
      aero_time_interp(team, t0, ts, tlen, y0_sub, y1_sub, yout_sub);
    });
} // calculate_current_data_concentrations
/*-----------------------------------------------------------------*/

} // namespace spa
} // namespace scream

#endif // SPA_FUNCTIONS_IMPL_HPP
