#include "p3_upwind_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing upwind functions on Reals using the
 * default device.
 */

#define ETI_UPWIND(nfield)                                              \
  template void Functions<Real,DefaultDevice>                           \
  ::calc_first_order_upwind_step<nfield>(                               \
    const uview_1d<const Spack>& rho,                                   \
    const uview_1d<const Spack>& inv_rho,                               \
    const uview_1d<const Spack>& inv_dz,                               \
    const MemberType& team,                                             \
    const Int& nk, const Int& k_bot, const Int& k_top, const Int& kdir, \
    const Scalar& dt_sub,                                               \
    const view_1d_ptr_array<Spack, nfield>& flux,                       \
    const view_1d_ptr_array<Spack, nfield>& V,                          \
    const view_1d_ptr_array<Spack, nfield>& r);
ETI_UPWIND(1)
ETI_UPWIND(2)
ETI_UPWIND(4)
#undef ETI_UPWIND

#define ETI_GENSED(nfield)                                              \
  template void Functions<Real,DefaultDevice>                           \
  ::generalized_sedimentation<nfield>(                                  \
    const uview_1d<const Spack>& rho,                                   \
    const uview_1d<const Spack>& inv_rho,                               \
    const uview_1d<const Spack>& inv_dz,                               \
    const MemberType& team,                                             \
    const Int& nk, const Int& k_qxtop, Int& k_qxbot, const Int& kbot, const Int& kdir, \
    const Scalar& Co_max, Scalar& dt_left, Scalar& prt_accum,           \
    const view_1d_ptr_array<Spack, nfield>& flux,                       \
    const view_1d_ptr_array<Spack, nfield>& V,                          \
    const view_1d_ptr_array<Spack, nfield>& r);
ETI_GENSED(1)
ETI_GENSED(2)
ETI_GENSED(4)
#undef ETI_GENSED

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
