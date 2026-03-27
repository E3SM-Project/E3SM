#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<float>(const float, const Field&);

} // namespace scream
