#include "shoc_assumed_pdf_impl.hpp"

namespace scream {
namespace shoc {

/*
 * Explicit instantiation for the default device.
 */

template struct Functions<Real,DefaultDevice>;

} // namespace shoc
} // namespace scream
