#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<int>(const int, const Field&);

} // namespace scream
