#include "impl/zm_input_state_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for ZmInputState struct on Reals using the
 * default device.
 */

template void Functions<Real,DefaultDevice>::ZmInputState::transpose<ekat::TransposeDirection::f2c>(int ncol, int nlev_mid);
template void Functions<Real,DefaultDevice>::ZmInputState::transpose<ekat::TransposeDirection::c2f>(int ncol, int nlev_mid);

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
