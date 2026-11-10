#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::deep_copy_masked<false,double>(const double, const Field&) const;
template void Field::deep_copy_masked<true, double>(const double, const Field&) const;

} // namespace scream
