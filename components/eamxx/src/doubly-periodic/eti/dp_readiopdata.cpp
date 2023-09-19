#include "impl/dp_readiopdata_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing readiopdata on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
