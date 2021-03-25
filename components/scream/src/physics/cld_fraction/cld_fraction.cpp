#include "physics/cld_fraction/cld_fraction_main_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace cldfrac {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct CldFractionFunctions<Real,DefaultDevice>;

} // namespace cldfrac
} // namespace scream
