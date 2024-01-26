#include "impl/dp_iop_default_opts_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing iop_default_opts on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
