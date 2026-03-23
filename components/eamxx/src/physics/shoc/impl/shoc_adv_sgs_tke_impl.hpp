#ifndef SHOC_ADV_SGS_TKE_IMPL_HPP
#define SHOC_ADV_SGS_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc adv_sgs_tke. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::adv_sgs_tke(
  const MemberType&            team,
  const Int&                   nlev,
  const Real&                  dtime,
  const bool&                  shoc_1p5tke,
  const bool&                  do_3d_turb,
  const uview_1d<const Pack>& shoc_mix,
  const uview_1d<const Pack>& wthv_sec,
  const uview_1d<const Pack>& sterm_zt,
  const uview_1d<const Pack>& tk,
  const uview_1d<const Pack>& brunt,
  const uview_1d<const Pack>& strain2,
  const uview_1d<Pack>&       tke,
  const uview_1d<Pack>&       a_diss)
{

  //Shared constants
  static constexpr Scalar ggr      = C::gravit.value;
  static constexpr Scalar basetemp = C::basetemp;
  static constexpr Scalar mintke   = scream::shoc::Constants<Real>::mintke;
  static constexpr Scalar maxtke   = scream::shoc::Constants<Real>::maxtke;
  Pack a_prod_bu;

  //declare some constants
  static constexpr Scalar Cs  = 0.15;
  static constexpr Scalar Ck  = 0.1;
  static constexpr Scalar Ce  = (Ck*Ck*Ck)/((Cs*Cs)*(Cs*Cs));
  static constexpr Scalar Ce1 = Ce/sp(0.7)*sp(0.19);
  static constexpr Scalar Ce2 = Ce/sp(0.7)*sp(0.51);
  static constexpr Scalar Cee = Ce1 + Ce2;

  const Int nlev_pack = ekat::npack<Pack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {

    // Compute buoyant production term
    if (shoc_1p5tke){
       // If there is no SGS variability the buoyant production term is closed
       //   as a function of the local moist brunt vaisalla frequency.
       a_prod_bu = -tke(k)*brunt(k);
    }
    else{
       a_prod_bu = (ggr/basetemp)*wthv_sec(k);
    }

    tke(k) = ekat::max(0,tke(k));

    // Shear production term, use diffusivity from previous timestep
    if (do_3d_turb){
      const Pack a_prod_sh = tk(k)*sterm_zt(k);
    }
    else{
      const Pack a_prod_sh = tk(k)*strain2(k);
    }

    // Dissipation term
    a_diss(k)=Cee/shoc_mix(k)*ekat::pow(tke(k),sp(1.5));

    // March equation forward one timestep
    tke(k)=ekat::max(mintke,tke(k)+dtime*(ekat::max(0,a_prod_sh+a_prod_bu)-a_diss(k)));

    tke(k)=ekat::min(tke(k),maxtke);
  });

}

} // namespace shoc
} // namespace scream

#endif
