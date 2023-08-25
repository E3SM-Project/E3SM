#include "impl/dp_advance_iop_forcing_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing advance_iop_forcing on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
