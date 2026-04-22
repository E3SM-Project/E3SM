#include "impl/zm_compute_dilute_parcel_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing compute_dilute_parcel on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
