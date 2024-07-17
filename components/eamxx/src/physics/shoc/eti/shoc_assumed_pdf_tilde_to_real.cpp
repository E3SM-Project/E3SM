#include "shoc_assumed_pdf_tilde_to_real_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing shoc_assumed_pdf_tilde_to_real on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
