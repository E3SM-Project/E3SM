#include "share/field/field.hpp"
#include "share/field/field_update_fill_aware.hpp"

namespace scream {

template void Field::update_fill_aware<CombineMode::Replace,  float, double>(const Field&, const float, const float, const float) const;
template void Field::update_fill_aware<CombineMode::Update,   float, double>(const Field&, const float, const float, const float) const;
template void Field::update_fill_aware<CombineMode::Multiply, float, double>(const Field&, const float, const float, const float) const;
template void Field::update_fill_aware<CombineMode::Divide,   float, double>(const Field&, const float, const float, const float) const;
template void Field::update_fill_aware<CombineMode::Max,      float, double>(const Field&, const float, const float, const float) const;
template void Field::update_fill_aware<CombineMode::Min,      float, double>(const Field&, const float, const float, const float) const;

} // namespace scream
