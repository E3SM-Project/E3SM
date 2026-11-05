#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Multiply, float, double>(const Field&, const float, const float, const float) const;

} // namespace scream
