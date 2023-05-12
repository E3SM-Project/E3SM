#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::shoc_diag_obklen_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const view_1d<const Scalar>& uw_sfc,
  const view_1d<const Scalar>& vw_sfc,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const view_2d<const Scalar>& thl_sfc,
  const view_2d<const Scalar>& cldliq_sfc,
  const view_2d<const Scalar>& qv_sfc,
  const view_1d<Scalar>&       ustar,
  const view_1d<Scalar>&       kbfs,
  const view_1d<Scalar>&       obklen)
{
  Kokkos::parallel_for(shcol, KOKKOS_LAMBDA(const Int& i) {
    shoc_diag_obklen(uw_sfc(i), vw_sfc(i), wthl_sfc(i), wqw_sfc(i),
                     ekat::subview(thl_sfc, i)(nlev-1),
                     ekat::subview(cldliq_sfc, i)(nlev-1),
                     ekat::subview(qv_sfc, i)(nlev-1),
                     ustar(i), kbfs(i), obklen(i));
  });
}

} // namespace shoc
} // namespace scream
