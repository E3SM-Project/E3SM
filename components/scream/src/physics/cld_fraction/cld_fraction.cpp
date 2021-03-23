#include "physics/cld_fraction/cldfraction_main_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace cldfrac {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace cldfrac
} // namespace scream
