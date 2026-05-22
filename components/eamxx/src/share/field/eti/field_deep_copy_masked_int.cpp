#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<false,int>(const int, const Field&) const;
template void Field::deep_copy_masked<true, int>(const int, const Field&) const;

} // namespace scream
