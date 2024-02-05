#include "impl/dp_iop_intht_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing iop_intht on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
