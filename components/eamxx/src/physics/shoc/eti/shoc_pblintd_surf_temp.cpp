#include "shoc_pblintd_surf_temp_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing pblintd_surf_temp on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
