#include "shoc_check_tke_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
