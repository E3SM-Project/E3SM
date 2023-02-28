#include "shoc_adv_sgs_tke_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing adv_sgs_tke on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
