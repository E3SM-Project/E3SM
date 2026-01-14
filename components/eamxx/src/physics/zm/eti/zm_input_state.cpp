#include "impl/zm_input_state_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_input_state on Reals using the
 * default device.
 */

template void Functions<Real,DefaultDevice>::zm_input_state::transpose<ekat::TransposeDirection::f2c>(int ncol, int nlev_mid);
template void Functions<Real,DefaultDevice>::zm_input_state::transpose<ekat::TransposeDirection::c2f>(int ncol, int nlev_mid);

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
