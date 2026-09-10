#ifndef SHOC_ASSEMBLE_SHOC_SHEAR_STRAIN3D_IMPL_HPP
#define SHOC_ASSEMBLE_SHOC_SHEAR_STRAIN3D_IMPL_HPP

#include "shoc_functions.hpp"

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::assemble_shoc_shear_strain3d(
  const MemberType&              team,
  const Int&                     nlev,
  const uview_2d<const Pack>&    shear_strain3d_components,
  const uview_1d<const Pack>&    du_dz_m,
  const uview_1d<const Pack>&    dv_dz_m,
  const uview_1d<const Pack>&    dw_dz_m,
  const uview_1d<Pack>&          shear_strain3d)
{
  const Int nlev_pack = ekat::npack<Pack>(nlev);

  // Assemble the full local velocity-gradient tensor from dycore horizontal
  // components and SHOC-computed vertical components, then form the symmetric strain invariant.
  team.team_barrier();
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    constexpr Scalar one_half = Scalar(0.5);
    constexpr Scalar two = Scalar(2.0);

    const Pack A00 = shear_strain3d_components(0,k);
    const Pack A01 = shear_strain3d_components(1,k);
    const Pack A10 = shear_strain3d_components(2,k);
    const Pack A11 = shear_strain3d_components(3,k);
    const Pack A20 = shear_strain3d_components(4,k);
    const Pack A21 = shear_strain3d_components(5,k);

    const Pack A02 = du_dz_m(k);
    const Pack A12 = dv_dz_m(k);
    const Pack A22 = dw_dz_m(k);

    const Pack S00 = A00;
    const Pack S11 = A11;
    const Pack S22 = A22;
    const Pack S01 = one_half * (A01 + A10);
    const Pack S02 = one_half * (A02 + A20);
    const Pack S12 = one_half * (A12 + A21);

    shear_strain3d(k) =
      two * (S00*S00 + S11*S11 + S22*S22
           + two*S01*S01 + two*S02*S02 + two*S12*S12);
  });
}

} // namespace shoc
} // namespace scream

#endif
