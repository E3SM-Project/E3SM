#include "impl/zm_calc_output_tend_impl.hpp"

namespace scream {
namespace zm {

/*
 * Explicit instantiation for doing zm_calc_output_tend on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace zm
} // namespace scream
