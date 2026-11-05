#include "share/field/field.hpp"
#include "share/field/field_update_masked.hpp"

namespace scream {

template void Field::update_masked<CombineMode::Replace, float, float>(const Field&, const float, const float, const float, const Field&) const;

} // namespace scream
