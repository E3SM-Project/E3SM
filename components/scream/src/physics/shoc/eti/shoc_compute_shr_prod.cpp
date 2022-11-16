#include "shoc_compute_shr_prod_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing compute_shr_prod on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
