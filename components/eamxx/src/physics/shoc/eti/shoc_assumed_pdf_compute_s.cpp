#include "shoc_assumed_pdf_compute_s_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing shoc_assumed_pdf_compute_s on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
