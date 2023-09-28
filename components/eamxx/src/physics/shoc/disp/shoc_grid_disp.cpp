#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_grid_disp(
  const Int&                  shcol,
  const Int&                  nlev,
  const Int&                  nlevi,
  const view_2d<const Spack>& zt_grid,
  const view_2d<const Spack>& zi_grid,
  const view_2d<const Spack>& pdel,
  const view_2d<Spack>&       dz_zt,
  const view_2d<Spack>&       dz_zi,
  const view_2d<Spack>&       rho_zt)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    shoc_grid(team, nlev, nlevi,
              ekat::subview(zt_grid, i),
              ekat::subview(zi_grid, i),
              ekat::subview(pdel, i),
              ekat::subview(dz_zt, i),
              ekat::subview(dz_zi, i),
              ekat::subview(rho_zt, i));
  });
}

} // namespace shoc
} // namespace scream
