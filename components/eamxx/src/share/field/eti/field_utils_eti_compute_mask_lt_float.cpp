#include "share/field/field_utils_impl.hpp"

namespace scream {
namespace impl {

template void compute_mask<Comparison::LT,float>(const Field&, float, Field&);

} // namespace impl
} // namespace scream
