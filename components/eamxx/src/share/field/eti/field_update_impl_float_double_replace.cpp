#include "share/field/field.hpp"
#include "share/field/field_update_impl.hpp"

namespace scream {

template void Field::update_impl<CombineMode::Replace, float, double>(const Field&, const float, const float, const float) const;

} // namespace scream
