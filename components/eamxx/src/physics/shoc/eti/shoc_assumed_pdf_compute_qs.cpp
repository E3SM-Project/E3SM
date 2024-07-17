#include "shoc_assumed_pdf_compute_qs_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing shoc_assumed_pdf_compute_qs on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
