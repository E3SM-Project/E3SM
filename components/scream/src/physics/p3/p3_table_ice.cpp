#include "p3_table_ice_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for doing ice table functions on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
