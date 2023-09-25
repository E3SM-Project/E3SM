#include "p3_rain_self_collection_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing rain self collection function on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
