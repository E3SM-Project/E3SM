#include "impl/gw_vd_lu_decomp_impl.hpp"

namespace scream {
namespace gw {

/*
 * Explicit instantiation for doing vd_lu_decomp on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace gw
} // namespace scream
