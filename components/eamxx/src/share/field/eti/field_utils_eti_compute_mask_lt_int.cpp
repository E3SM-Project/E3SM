#include "share/field/field_utils_impl.hpp"

namespace scream {
namespace impl {

template void compute_mask<Comparison::LT,int>(const Field&, int, Field&);

} // namespace impl
} // namespace scream
