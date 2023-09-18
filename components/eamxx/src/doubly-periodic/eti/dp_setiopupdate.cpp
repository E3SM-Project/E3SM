#include "impl/dp_setiopupdate_impl.hpp"

namespace scream {
namespace dp {

/*
 * Explicit instantiation for doing setiopupdate on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace dp
} // namespace scream
