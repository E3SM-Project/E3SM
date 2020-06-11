#include "p3_functions_dsd2_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instatiation for doing p3 dsd2 functions on Reals using the
 * default device.
 */

#define ETI_CLOUD_DSD2(zero_out)                                        \
template void Functions<Real,DefaultDevice>                             \
  ::get_cloud_dsd2<zero_out>(                                           \
    const Smask& qc_gt_small,const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu, \
    const view_dnu_table& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, const Spack& lcldm);
ETI_CLOUD_DSD2(true)
ETI_CLOUD_DSD2(false)
#undef ETI_CLOUD_DSD2

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
