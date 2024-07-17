#include "shoc_assumed_pdf_compute_temperature_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for doing shoc_assumed_pdf_compute_temperature on Reals using the
 * default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
