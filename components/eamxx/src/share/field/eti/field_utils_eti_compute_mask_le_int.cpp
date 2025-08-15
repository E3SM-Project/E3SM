#include "share/field/field_utils_impl.hpp"

namespace scream {
namespace impl {

template void compute_mask<Comparison::LE,int>(const Field&, int, Field&);

} // namespace impl
} // namespace scream
