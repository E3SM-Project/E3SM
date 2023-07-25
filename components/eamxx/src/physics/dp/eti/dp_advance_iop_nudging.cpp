#include "impl/dp_advance_iop_nudging_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing advance_iop_nudging on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
