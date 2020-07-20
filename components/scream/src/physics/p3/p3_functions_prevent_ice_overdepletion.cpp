#include "p3_functions_prevent_ice_overdepletion_impl.hpp"
#include "ekat/scream_types.hpp"

namespace scream {
namespace p3 {

/*
 * Explicit instantiation for preventing ice overdepletion.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace p3
} // namespace scream
