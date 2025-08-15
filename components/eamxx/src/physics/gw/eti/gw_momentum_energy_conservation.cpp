#include "impl/gw_momentum_energy_conservation_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing momentum_energy_conservation on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
