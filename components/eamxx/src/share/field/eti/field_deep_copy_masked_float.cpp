#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<false,float>(const float, const Field&) const;
template void Field::deep_copy_masked<true, float>(const float, const Field&) const;

} // namespace scream
