#include "impl/zm_output_tend_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for ZmOutputTent struct with Reals using the
 * default device.
 */

template void Functions<Real,DefaultDevice>::ZmOutputTend::transpose<ekat::TransposeDirection::f2c>(int ncol, int nlev_mid);
template void Functions<Real,DefaultDevice>::ZmOutputTend::transpose<ekat::TransposeDirection::c2f>(int ncol, int nlev_mid);

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
