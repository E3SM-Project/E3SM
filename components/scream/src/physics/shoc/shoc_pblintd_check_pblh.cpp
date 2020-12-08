#include "shoc_pblintd_check_pblh_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing pblintd_check_pblh on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
