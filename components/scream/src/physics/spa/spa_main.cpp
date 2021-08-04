#include "physics/spa/spa_main_impl.hpp"
#include "share/scream_types.hpp"

namespace scream {
namespace spa {

/*
 * Explicit instantiation for doing main on Reals using the
 * default device.
 */

template struct SPAFunctions<Real,DefaultDevice>;

} // namespace spa 
} // namespace scream
